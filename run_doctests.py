
if __name__ == '__main__':
    import doctest
    from lesion import lifio, stats, trace
    map(doctest.testmod, [lifio, stats, trace])

