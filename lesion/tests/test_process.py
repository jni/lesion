import numpy as np
from lesion import process


def test_bad_image():
    im = np.zeros((64, 64), np.uint16)
    assert process.bad_image(im)
    im[32, :] = 6554
    assert not process.bad_image(im)
