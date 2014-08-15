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


def traces_dict(fin):
    """From a LIF file, produce image series, traces, stats.

    Parameters
    ----------
    fin : string
        The input filename.

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
    rdr = lifio.image_reader(fin)
    all_times = _times(names)
    positions, times = zip(*map(lifio.parse_series_name, names))
    traces = collections.OrderedDict()
    ntimes, npositions, nstats = map(len, [all_times, positions, all_stats])
    statistics = (np.empty((ntimes, npositions * nstats), dtype=np.float32) *
                  np.float32(np.nan))
    statistics = pd.DataFrame(statistics, index=all_times,
                              columns=it.product(positions, all_stat_names))

    for i, name in enumerate(names):
        images = lifio.read_image_series(rdr, i, desired_order='cztyx')[0]
        images2d = images.sum(axis=0) # squish z dimension
        position, times = lifio.parse_series_name(name)
        traces[position]['times'] = times
        current_traces = map(trace.trace_profile, images2d)
        traces[position]['traces'] = current_traces
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
    return np.concatenate([times[i] for i in indices])
