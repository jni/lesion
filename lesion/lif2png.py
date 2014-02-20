import argparse


def convert_files(files, nseries=18, channel=0, verbose=False):
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
    verbose : bool, optional
        Print out diagnostic information to stdout.
    """
    from jython_imports import IJ, BF, ImporterOptions, ZProjector
    files = sorted(filter(lambda f: f.endswith('.lif'), files))
    for fin in files:
        if verbose:
            print(fin)
        fout_base = fin[:-4] + '%02i-%s.png'
        opts = ImporterOptions()
        opts.setId(fin)
        opts.setUngroupFiles(True)
        opts.setOpenAllSeries(True)
        for i in range(nseries):
            # keep only the required channel for each series
            opts.setCBegin(i, channel)
            opts.setCEnd(i, channel) 
        imps = BF.openImagePlus(opts)
        projector = ZProjector()
        projector.setMethod(ZProjector.SUM_METHOD)
        for i, imp in enumerate(imps):
            projector.setImage(imp)
            projector.doProjection()
            impout = projector.getProjection()
            fout = fout_base % (i, imp.getTitle())
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
    parser.add_argument('-n', '--nseries', type=int, metavar='INT', default=18,
                        help='Number of series in each file.')
    parser.add_argument('-c', '--channel', type=int, metavar='INT', default=0,
                        help='Channel of interest.')

    args = parser.parse_args()
    convert_files(args.files, args.nseries, args.channel, args.verbose)

