"""
Compute statistics from linear traces of embryos. The most basic
statistic is min/max, but others will be compiled here.
"""

import numpy as np
from scipy import ndimage as nd


def min_max(tr):
    """Return the ratio of minimum to maximum of a trace.

    Parameters
    ----------
    tr : 1D array of float
        The input profile.

    Returns
    -------
    mm : float
        The ratio of the minimum value in `tr` over the maximum.

    Examples
    --------
    >>> tr = np.array([0.8, 0.9, 1.4, 2.0, 1.1])
    >>> min_max(tr) # 0.8 / 2.0
    0.4
    """
    tr = tr.astype(float)
    mm = tr.min() / tr.max()
    return mm


def slope(tr, sigma=None):
    """Compute the absolute slope between the max and min positions.

    Parameters
    ----------
    tr : 1D array of float
        The input profile.
    sigma : float, optional
        Smooth `tr` by a Gaussian filter with this sigma.

    Returns
    -------
    a : float
        The slope from max to min, in absolute value.

    Examples
    --------
    >>> tr = np.array([5, 5, 5, 0, 2, 1, 0, 5, 5, 5])
    >>> slope(tr)
    1.6666666666666667
    >>> slope(tr, sigma=1)
    0.75565422533672888
    """
    tr = tr.astype(float)
    if sigma is not None:
        tr = nd.gaussian_filter1d(tr, sigma=sigma)
    m, M = np.argmin(tr), np.argmax(tr)
    a = np.abs((tr[m] - tr[M]) / (m - M))
    return a
