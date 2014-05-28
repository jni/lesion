import argparse


def get_series_properties(title):
    """Parse a series title to extract the relevant info.

    Parameters
    ----------
    title : string
        The title of the series.

    Returns
    -------
    props : string
        The relevant part of the title.

    Examples
    --------
    >>> title = 'MacAb 3dpf S001--.5dpSCI.lif - MacAb 3dpf S001-_C-E004'
    >>> get_series_properties(title)
    'S001-_C-E004'
    """
    props = title.split(' ')[-1]
    return props


def convert_files(files, nseries=18, channel=0, virtual=False, verbose=False):
    """Extract series from lif files, flatten, and save as png.

    Parameters
    ----------
    files : list of string
        A list of filenames to be converted.
    nseries : int, optional
        The number of image series within each file.
    channel : int, optional
        The channel within each series containing the image of
        interest.
    virtual : bool, optional
        Open a virtual stack instead of a stack. Saves precious memory.
    verbose : bool, optional
        Print out diagnostic information to stdout.
    """
    from jython_imports import (IJ, BF, ImporterOptions,
                                ImageConverter, ZProjector)
    files = sorted(filter(lambda f: f.endswith('.lif'), files))
    projector = ZProjector()
    projector.setMethod(ZProjector.SUM_METHOD)
    for fin in files:
        if verbose:
            print(fin)
        fout_base = fin[:-4] + '-%02i-%s.png'
        opts = ImporterOptions()
        opts.setId(fin)
        opts.setUngroupFiles(True)
        if virtual:
            opts.setVirtual(True)
        else:
            opts.setOpenAllSeries(True)
        for i in range(nseries):
            # keep only the required channel for each series
            opts.setCBegin(i, channel)
            opts.setCEnd(i, channel) 
        imps = BF.openImagePlus(opts)
        for i, imp in enumerate(imps):
            projector.setImage(imp)
            projector.doProjection()
            impout = projector.getProjection()
            converter = ImageConverter(impout)
            converter.convertToGray16()
            fout = fout_base % (i, get_series_properties(imp.getTitle()))
            if verbose:
                print "creating", fout
            IJ.saveAs(impout, 'png', fout)
            imp.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Extract relevant info from lif series files.')
    parser.add_argument('files', help='The files to be converted.', nargs='+')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Print out runtime information.')
    parser.add_argument('-V', '--virtual', action='store_true',
                        help='Open virtual stack.')
    parser.add_argument('-n', '--nseries', type=int, metavar='INT', default=18,
                        help='Number of series in each file.')
    parser.add_argument('-c', '--channel', type=int, metavar='INT', default=0,
                        help='Channel of interest.')

    args = parser.parse_args()
    convert_files(args.files, args.nseries, args.channel,
                  args.virtual, args.verbose)

