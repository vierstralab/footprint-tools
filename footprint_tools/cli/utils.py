def chunkify(l, n):
	"""Splits an iterable list into n chunks"""
	return [l[i::nchunks] for i in range(nchunks)]

def add_log_filehandler(logger, outfile, name='stdout'):
	pass

from argh.exceptions import CommandError

def tuple_ints(arg):
	"""
	Function to parser a tuple of integers from command line
	"""
	try:
		fw, rev = list(map(int, arg.split(',')))
		return (fw, rev)
	except:
		raise CommandError("Argument must be a comma-delimited tuple -- i.e., 0,-1")

def list_ints(arg):
	"""
	Function to parser a list of integers from command line
	"""
	try:
		return list(map(int, arg.split(',')))
	except:
		raise CommandError("Argument must be a comma-delimited list -- i.e., 0,1,2,3")
