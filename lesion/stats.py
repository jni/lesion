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
    >>> min_max(tr) # doctest: +ELLIPSIS
    0.4...
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
    >>> slope(tr) # doctest: +ELLIPSIS
    1.6666...
    >>> slope(tr, sigma=1) # doctest: +ELLIPSIS
    0.7556...
    """
    tr = tr.astype(float)
    if sigma is not None:
        tr = nd.gaussian_filter1d(tr, sigma=sigma)
    m, M = np.argmin(tr), np.argmax(tr)
    a = np.abs((tr[m] - tr[M]) / (m - M))
    return a


def missing_fluorescence(tr, sigma=None, height=None, margins=50):
    """Compute amount of fluorescence lost compared with no injury.

    Parameters
    ----------
    tr : 1D array of float
        The input profile.
    sigma : float, optional
        Smooth `tr` by a Gaussian filter with this sigma.
    height : float, optional
        The base level from which to estimate missed fluorescence. If
        this is not provided, it is estimated from the values near the
        edge of the profile.
    margins : int, optional
        How much of the edge of the profile is used to compute the
        normal fluorescence intensity of the profile.

        This is ignored if `height` is given.

    Returns
    -------
    m : float
        The total missing fluorescence under the trace.

    Examples
    --------
    >>> tr = np.array([3, 5, 4, 0, 2, 2, 0, 3, 4, 5])
    >>> missing_fluorescence(tr, margins=3) # 1 + 4 + 2 + 2 + 4 + 1
    14.0
    >>> missing_fluorescence(tr, height=3) # 3 + 1 + 1 + 3
    8.0
    """
    tr = tr.astype(float)
    if height is None:
        height = np.mean(tr[:margins] + tr[-margins:]) / 2
    if sigma is not None:
        tr = nd.gaussian_filter1d(tr, sigma=sigma)
    tr = np.clip(tr, 0, height)
    m = (height - tr).sum()
    return m
