import os
from lesion import lifio

import numpy as np

from numpy.testing.decorators import skipif
from numpy.testing import assert_equal, assert_allclose, assert_raises

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


@skipif(test_lif_unavailable)
def test_read_image_series():
    rdr0 = lifio.image_reader(test_lif)
    im0 = lifio.read_image_series(rdr0, desired_order=None)
    im0 = im0.transpose((0, 1, 3, 4, 2)) # native order is tzcyx
    im1 = lifio.read_image_series(test_lif, desired_order='ctzyx')
    im1 = im1.transpose((1, 2, 3, 4, 0))
    assert_equal(im0, im1)


@skipif(test_lif_unavailable)
def test_read_one_channel():
    im0 = lifio.read_image_series(test_lif, z=0, c=0)
    assert_equal(np.squeeze(im0).shape, (512, 512))


@skipif(test_lif_unavailable)
def test_series_iterator():
    import collections as coll
    series_iter = lifio.series_iterator(test_lif)
    assert isinstance(series_iter, coll.Iterator)
    for series in series_iter:
        assert series.shape[-3:] == (4, 512, 512)


@skipif(test_lif_unavailable)
def test_invalid_series_id():
    # test LIF has only 2 series (0 and 1).
    assert_raises(ValueError, lifio.read_image_series, test_lif, series_id=5)


def test_done():
    lifio.done()
    assert lifio.VM_KILLED


def test_metadata_raise_error():
    assert_raises(RuntimeError, lifio.metadata, test_lif)


def test_image_reader_killed_error():
    assert_raises(RuntimeError, lifio.image_reader, test_lif)


def test_bad_series_name():
    sname = 'Totally wrong string'
    assert_raises(ValueError, lifio.parse_series_name, sname)
