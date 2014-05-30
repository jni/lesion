import numpy as np
import javabridge as jv
import bioformats as bf
from xml import etree as et

def start(max_heap_size='8G'):
    """Start the Java Virtual Machine, enabling bioformats IO.

    Parameters
    ----------
    max_heap_size : string, optional
        The maximum memory usage by the virtual machine. Valid strings
        include '256M', '64k', and '2G'. Expect to need a lot.
    """
    jv.start_vm(class_path=bf.JARS, max_heap_size=max_heap_size)

def done():
    """Kill the JVM. Once killed, it cannot be restarted.

    Notes
    -----
    See the python-javabridge documentation for more information.
    """
    jv.kill_vm()


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


def parse_xml_metadata(xml_string, array_order='tzyxc'):
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
