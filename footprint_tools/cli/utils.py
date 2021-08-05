def chunkify(l, n):
    """Splits an iterable list into n chunks"""
	return [l[i::nchunks] for i in range(nchunks)]

def add_log_filehandler(logger, outfile, name='stdout'):
    pass