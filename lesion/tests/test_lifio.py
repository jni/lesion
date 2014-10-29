import os
from lesion import lifio

from numpy.testing.decorators import skipif
from numpy.testing import assert_equal, assert_allclose

currdir = os.path.abspath(os.path.dirname(__file__))

test_lif = os.path.join(currdir, 'mouse-kidney.lif')
test_lif_unavailable = not os.path.isfile(test_lif)

@skipif(test_lif_unavailable)
def test_metadata_size():
    assert_equal(lifio.lif_metadata_string_size(test_lif), 107192)

@skipif(test_lif_unavailable)
def test_metadata():
    names, sizes, reso = lifio.metadata(test_lif)
    assert_equal(names, ['Series016', 'Series019'])
    assert_equal(sizes, [(1, 25, 512, 512, 4), (1, 46, 512, 512, 4)])
    assert_allclose(reso, [(0.9999, 1.5137, 1.5137),
                           (0.1395, 0.2539, 0.2539)], atol=1e-4)

def test_done():
    lifio.done()
    assert lifio.VM_KILLED
