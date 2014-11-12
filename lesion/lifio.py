import numpy as np
import re
import collections as coll
import itertools as it
import javabridge as jv
import bioformats as bf
from xml import etree as et


VM_STARTED = False
VM_KILLED = False
DEFAULT_DIM_ORDER = 'tzyxc'


BF2NP_DTYPE = {
    0: np.int8,
    1: np.uint8,
    2: np.int16,
    3: np.uint16,
    4: np.int32,
    5: np.uint32,
    6: np.float32,
    7: np.double
}


def start(max_heap_size='8G'):
    """Start the Java Virtual Machine, enabling bioformats IO.

    Parameters
    ----------
    max_heap_size : string, optional
        The maximum memory usage by the virtual machine. Valid strings
        include '256M', '64k', and '2G'. Expect to need a lot.
    """
    jv.start_vm(class_path=bf.JARS, max_heap_size=max_heap_size)
    global VM_STARTED
    VM_STARTED = True

def done():
    """Kill the JVM. Once killed, it cannot be restarted.

    Notes
    -----
    See the python-javabridge documentation for more information.
    """
    jv.kill_vm()
    global VM_KILLED
    VM_KILLED = True


def lif_metadata_string_size(filename):
    """Get the length in bytes of the metadata string of a LIF file.

    Parameters
    ----------
    filename : string
        Path to the LIF file.

    Returns
    -------
    length : int
        The length in bytes of the metadata string.

    Notes
    -----
    This is based on code by Lee Kamentsky. [1]

    References
    ----------
    [1] https://github.com/CellProfiler/python-bioformats/issues/8
    """
    with open(filename, 'rb') as fd:
        fd.read(9)
        length = np.frombuffer(fd.read(4), "<i4")[0]
        return length


def parse_xml_metadata(xml_string, array_order=DEFAULT_DIM_ORDER):
    """Get interesting metadata from the LIF file XML string.

    Parameters
    ----------
    xml_string : string
        The string containing the XML data.
    array_order : string
        The order of the dimensions in the multidimensional array.
        Valid orders are a permutation of "tzyxc" for time, the three
        spatial dimensions, and channels.

    Returns
    -------
    names : list of string
        The name of each image series.
    sizes : list of tuple of int
        The pixel size in the specified order of each series.
    resolutions : list of tuple of float
        The resolution of each series in the order given by
        `array_order`. Time and channel dimensions are ignored.
    """
    array_order = array_order.upper()
    names, sizes, resolutions = [], [], []
    spatial_array_order = [c for c in array_order if c in 'XYZ']
    size_tags = ['Size' + c for c in array_order]
    res_tags = ['PhysicalSize' + c for c in spatial_array_order]
    metadata_root = et.ElementTree.fromstring(xml_string)
    for child in metadata_root:
        if child.tag.endswith('Image'):
            names.append(child.attrib['Name'])
            for grandchild in child:
                if grandchild.tag.endswith('Pixels'):
                    att = grandchild.attrib
                    sizes.append(tuple([int(att[t]) for t in size_tags]))
                    resolutions.append(tuple([float(att[t])
                                              for t in res_tags]))
    return names, sizes, resolutions


def metadata(filename, array_order=DEFAULT_DIM_ORDER):
    """Get metadata from a BioFormats file.

    Parameters
    ----------
    filename : string
        The path to a BioFormats File.
    array_order : string
        The order of the dimensions in the multidimensional array.
        Valid orders are a permutation of "tzyxc" for time, the three
        spatial dimensions, and channels.

    Returns
    -------
    names : list of string
        The name of each image series.
    sizes : list of tuple of int
        The pixel size in the specified order of each series.
    resolutions : list of tuple of float
        The resolution of each series in the order given by
        `array_order`. Time and channel dimensions are ignored.
    """
    if not VM_STARTED:
        start()
    if VM_KILLED:
        raise RuntimeError("The Java Virtual Machine has already been "
                           "killed, and cannot be restarted. See the "
                           "python-javabridge documentation for more "
                           "information. You must restart your program "
                           "and try again.")
    md_string = bf.get_omexml_metadata(filename)
    return parse_xml_metadata(md_string)


def parse_series_name(name, interval=0.5):
    """Get series information (time points, embryo) from the name.

    Parameters
    ----------
    name : string
        The name of the image series.
    interval : float, optional
        The time interval between timepoints.

    Returns
    -------
    embryo : int
        The embryo ID / position number on the slide.
    times : array of float
        The timepoint corresponding to each image.

    Notes
    -----
    A "Pre" tag in the name is interpreted as ``t = -1``. A "Mark" tag
    is interpreted as ``t = -2``.

    Examples
    --------
    >>> name0 = "Pre lesion 2x/Pos008_S001"
    >>> e0, t0 = parse_series_name(name0)
    >>> e0, t0[:5]
    (8, array([-1.]))
    >>> name1 = "Mark_and_Find_001/Pos019_S001"
    >>> e1, t1 = parse_series_name(name1)
    >>> e1, t1[:5]
    (19, array([-2.]))
    >>> name2 = "22.5h to 41h pSCI/Pos013_S001"
    >>> e2, t2 = parse_series_name(name2)
    >>> e2, t2[:5]
    (13, array([ 22.5,  23. ,  23.5,  24. ,  24.5]))
    >>> name3 = "63 to 72.5hpSCI/Pos004_S001"
    >>> e3, t3 = parse_series_name(name3, interval=1)
    >>> e3, t3[:5]
    (4, array([ 63.,  64.,  65.,  66.,  67.]))
    """
    re_pre = r"Pre.*/Pos(\d+)_.*"
    re_mark = r"Mark.*/Pos(\d+)_.*"
    re_time = r"(\d+(\.\d+)?)h? to (\d+(\.\d+)?)h.*/Pos(\d+)_.*"
    m_pre = re.match(re_pre, name)
    m_mark = re.match(re_mark, name)
    m_time = re.match(re_time, name)
    if m_pre:
        return int(m_pre.groups()[0]), np.array([-1.])
    elif m_mark:
        return int(m_mark.groups()[0]), np.array([-2.])
    elif m_time:
        g = m_time.groups()
        t0 = float(g[0])
        t1 = float(g[2])
        embryo = int(g[4])
        times = np.arange(t0, t1 + float(interval)/2, interval)
        return embryo, times
    else:
        raise ValueError("Could not parse name string: %s" % name)


def image_reader(filelike):
    """Return a BioFormats ``ImageReader`` object from filelike.

    Parameters
    ----------
    filelike : string or bf.ImageReader
        If string, open the corresponding ImageReader. If ImageReader,
        this function is a no-op.

    Returns
    -------
    rdr : bf.ImageReader
        The relevant reader.

    Notes
    -----
    The purpose of this function is to provide a *robust* way to open
    a BioFormats file --- without having to start the JVM manually.
    """
    if not VM_STARTED:
        start()
    if VM_KILLED:
        raise RuntimeError("The Java Virtual Machine has already been "
                           "killed, and cannot be restarted. See the "
                           "python-javabridge documentation for more "
                           "information. You must restart your program "
                           "and try again.")
    if isinstance(filelike, bf.ImageReader):
        rdr = filelike
    else:
        rdr = bf.ImageReader(filelike)
    return rdr


def read_image_series(filelike, series_id=0, t=None, z=None, c=None,
                      desired_order=None):
    """Read an image volume from a file.

    Parameters
    ----------
    filelike : string or bf.ImageReader
        Either a filename containing a BioFormats image, or a
        `bioformats.ImageReader`.
    series_id : int, optional
        Load this series from the image file.
    t : int or list of int, optional
        Load this/these timepoint/s only from the image file. If 
        `None`, load all timepoints.
    z : int or list of int, optional
        Load this/these z-plane/s only from the image file. If `None`,
        load all planes.
    c : int or list of int, optional
        Load this/these channel/s only from the image file. If `None`,
        load all channels.
    desired_order : string, optional
        Store the image dimensions in this order, such as "TZCYX". By
        default, the order will be the exact inverse of the image's
        native order, since (Leica) files use Fortran order, while
        NumPy uses C order.

    Returns
    -------
    image : numpy ndarray, 5 dimensions
        The read image.
    """
    rdr = image_reader(filelike)
    reader = rdr.rdr
    total_series = reader.getSeriesCount()
    if not 0 <= series_id < total_series:
        raise ValueError("Series ID %i is not between 0 and the total "
                         "number of series, %i." % (series_id, total_series))
    reader.setSeries(series_id)
    order = reader.getDimensionOrder()
    # we invert the shape because numpy uses C order and (most?) BF images use
    # Fortran order
    old_shape = [getattr(reader, "getSize" + s)() for s in order]
    czt_list, old_shape = _sanitize_czt(c, z, t, old_shape, order)
    if desired_order is not None:
        desired_order = desired_order.upper()
        transposition = _get_ordering(order, desired_order)
        new_shape = [old_shape[i] for i in transposition]
    else:
        new_shape = old_shape[::-1]
        desired_order = order[::-1]
    image = np.empty(new_shape, dtype=BF2NP_DTYPE[reader.getPixelType()])
    for c, z, t in czt_list:
        indices = []
        for d in desired_order:
            if d == 'C':
                indices.append(c)
            elif d == 'Z':
                indices.append(z)
            elif d == 'T':
                indices.append(t)
            else:
                indices.append(slice(None))
        image[indices] = rdr.read(z=z, t=t, c=c, series=series_id,
                                  rescale=False)
    return image


def series_iterator(filelike, **kwargs):
    """Iterate over all the series in a file.

    Parameters
    ----------
    filelike : string or bf.ImageReader
        The input file.
    **kwargs : keyword arguments, optional
        Keyword arguments to be passed on to `read_image_series`.

    Returns
    -------
    seit : iterator
        Iterator over all series in `filelike`.
    """
    rdr = image_reader(filelike)
    total_series = rdr.rdr.getSeriesCount()
    for series_id in range(total_series):
        yield read_image_series(rdr, series_id, **kwargs)


def _get_ordering(actual, desired):
    """Find an ordering of indices so that desired[i] == actual[ordering[i]].

    Parameters
    ----------
    actual, desired : string, same length
        Two strings differring only in the permutation of their characters.

    Returns
    -------
    ordering : list of int
        A list of indices into `actual` such that
        ``desired[i] == actual[ordering[i]]``.

    Examples
    --------
    >>> actual = "XYCZT"
    >>> desired = "TZYXC"
    >>> _get_ordering(actual, desired)
    [4, 3, 1, 0, 2]
    """
    ordering = []
    for elem in desired:
        ordering.append(actual.find(elem))
    return ordering


def _sanitize_czt(c, z, t, shape, order):
    """Ensure cs, zs, and ts are lists, and the correct shape is returned.

    Parameters
    ----------
    c, z, t : int, list of int, or None
        The desired c, z, and t values for the image. `None` means all values
        are desired.
    shape : list of int
        The shape of the final image. *Will be edited in-place.*
    order : string
        The order of the dimensions in the file being read.

    Returns
    -------
    czt_list : list of (c, z, t) tuple
        A list of tuples to read in sequentially from the file.
    shape : list of int
        The shape of the final image, taking into account the subset of c, z,
        and t that will be read in.
    """
    out = {}
    for dim, label in zip((c, z, t), ('C', 'Z', 'T')):
        dim_idx = order.find(label)
        if dim is None:
            out[label] = range(shape[dim_idx])
        else:
            if not isinstance(dim, coll.Iterable):
                out[label] = [dim]
            shape[dim_idx] = len(out[label])
    czt_order = filter(lambda x: x in 'CZT', order)
    c, z, t = _get_ordering(czt_order, 'CZT')
    czt_list = [(tup[c], tup[z], tup[t])
                for tup in it.product(*[out[char] for char in czt_order])]
    return czt_list, shape
