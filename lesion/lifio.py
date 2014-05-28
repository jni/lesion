import numpy as np
import javabridge as jv
import bioformats as bf

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
