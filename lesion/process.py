"""
Process one or more image files into a set of statistical time series.
"""
import itertools as it
import collections

import numpy as np
import pandas as pd

from . import lifio
from . import trace
from . import stats

# constants
all_stats = [stats.min_max, stats.slope, stats.missing_fluorescence]
all_stat_names = ['min_max', 'slope', 'missing']


def bad_image(im):
    """Heuristic to determine if the image was a bad acquisition.

    Sometimes, due to microscope misalignment or the embryo dying, the
    input image contains little or no fluorescence, making downstream
    statistics noisy or misleading. This function returns ``True`` when
    its heuristics (see source) determine that the image is not good.

    Parameters
    ----------
    im : array, shape (M, N)
        The input image.

    Returns
    -------
    is_bad : bool
        ``True`` when `im` is bad.
    """
    is_bad = im.max() < 100
    return is_bad


def traces_dict(fin, series=None, chan=0, return_images=False):
    """From a LIF file, produce image series, traces, stats.

    Parameters
    ----------
    fin : string
        The input filename.
    series : list of int, optional
        Which series to process. ``None`` is interpreted as all
        available series.
    chan : int, optional
        The channel containing the image to be traced.
    return_images : bool, optional
        If ``True``, the images underlying the statistics are returned as
        part of the traces dictionary.

    Returns
    -------
    traces : OrderedDict
        A dictionary mapping positions to timepoints and a list of
        traces.
    statistics : pandas DataFrame
        The statistics, with rows for each timepoint and columns for
        each position and statistic.
    """
    names, sizes, resolutions = lifio.metadata(fin)
    if series is None:
        series = range(len(names))
    names, sizes, resolutions = map(lambda x: [x[i] for i in series],
                                    (names, sizes, resolutions))
    rdr = lifio.image_reader(fin)
    all_times = _times(names)
    positions, times = zip(*map(lifio.parse_series_name, names))
    traces = collections.OrderedDict()
    ntimes, npositions, nstats = map(len, [all_times, positions, all_stats])
    statistics = (np.empty((ntimes, npositions * nstats), dtype=np.float32) *
                  np.float32(np.nan))
    statistics = pd.DataFrame(statistics, index=all_times,
                              columns=it.product(positions, all_stat_names))
    image_series = lifio.series_iterator(rdr, series,
                                         desired_order='tzcyx', c=chan)

    for i, (name, images) in enumerate(zip(names, image_series)):
        images2d = np.squeeze(images.sum(axis=1, dtype=np.uint16)) # squash z
        if images2d.ndim == 2:
            images2d = images2d[np.newaxis, ...]
        position, times = lifio.parse_series_name(name)
        if not traces.has_key(position):
            traces[position] = {'times': [], 'traces': [], 'images': []}
        current_traces = map(trace.trace_profile, images2d)
        traces[position]['times'].extend(times)
        traces[position]['traces'].extend(current_traces)
        if return_images:
            traces[position]['images'].extend(images2d)
        for tr, time in zip(current_traces, times):
            for stat, stat_name in zip(all_stats, all_stat_names):
                statistics.loc[time][(position, stat_name)] = stat(tr)

    return traces, statistics


def _times(names):
    """Get the set of all unique timepoints represented in `names`.

    Parameters
    ----------
    names : list of string
        The names of a set of image series. The timepoints should be
        encoded in the names according to the rules encoded in the
        `lifio.parse_series_name` function.

    Returns
    -------
    times : array of float
        array of all image times (in hours).

    Notes
    -----
    This function assumes there is a one-to-one correspondence between
    start times and unique time series.
    """
    positions, times = zip(*map(lifio.parse_series_name, names))
    start_times = [time[0] for time in times]
    _, indices = np.unique(start_times, return_index=True)
    return np.unique(np.concatenate([times[i] for i in indices]))
