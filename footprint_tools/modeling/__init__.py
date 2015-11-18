__all__ = ["bias", "dispersion", "smoothing"]

from .. import cutcounts

import numpy as np
import pyfaidx

def reverse_complement(seq):
	"""
	Computes the reverse complement of a genomic sequence
	
	Parameters
	----------
	seq : string 
		FASTA DNA sequence
	
	Returns
	-------
	value: string
		Reverse complement of input sequence
	"""

	compl = { 'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'}
	return ''.join([ compl[base] for base in seq ])[::-1]

def predict_interval(reads, seq, interval, bm, half_window_width = 10, smoothing_class = None):
	"""
	Creates an expected distribution of cleavage counts within a genomic interval
	using observed data, a sequence preference model (bias), and a windowed smoothing
	function.
	
	Parameters
	----------
	reads :	bamfile (or equivalent class)
		Raw sequence alignment file
	seq : faidx
		FASTA index instance	
	interval : genomic_interval
		Genomic interval 
	bm : bias_model

	half_window_width : int

	smoothing_class : smoothing class or None

	Returns
	-------
	A dictionary
	"""
	
	obs_counts = {'+': None, '-': None}
	exp_counts = {'+': None, '-': None}
	win_counts = {'+': None, '-': None}

	# Pad the interval sequence by the half-windown width and the bm model offset
	padding = half_window_width + bm.offset
	interval_seq = seq[interval.chrom][interval.start-padding:interval.end+padding].seq.upper()

	# Pad the cut counts by the half window width, and the smoothing class half window width
	padding = half_window_width
	if smoothing_class:
		padding = padding + smoothing_class.half_window_width

	counts = reads[interval.widen(padding)]

	# Generalized function to predict a cleavage count
	predict_func = lambda x: bm.predict(x[1:], x[0])[half_window_width]

	for strand in ['+', '-']:

		# Pre-calculate the sequence bias propensity table from bias model
		if strand == '+':
			probs = bm.probs(interval_seq)
		else:
			probs = bm.probs(reverse_complement(interval_seq))[::-1]

		# Get indicies for overlaping 1-bp stepped windows
		win_idx = np.vstack( np.arange(i-half_window_width, i+half_window_width+1) for i in np.arange(half_window_width, len(counts[strand])-half_window_width) )

		# Count tags in window
		w_counts = np.apply_along_axis(np.sum, 1, counts[strand][win_idx])

		# Smooth windows if necesary
		if smoothing_class:
			win_counts[strand] = smoothing_class.smooth(w_counts, half_window_width)
		else:
			win_counts[strand] = w_counts
		
		z = np.concatenate((win_counts[strand][:, np.newaxis], probs[win_idx[0:len(interval), :]]), axis = 1)

		obs_counts[strand] = counts[strand][padding:-padding]
		exp_counts[strand] = np.apply_along_axis(predict_func, 1, z)

	return { "obs": obs_counts, "exp": exp_counts, "win": win_counts }


def build_histogram(reads, seq, intervals, bm, half_window_width = 5, size = (200, 1000)):
	"""Creates a histogram of expected vs. observed cleavages

	Returns
	-------
	(np.matrix): matrix of size (m, n) that holds counts used to fit the dispersion model
	"""

	h = np.zeros(size)
	for interval in intervals:
		res = modeling.predict_interval(reads, seq, interval, bm, half_window_width)

		obs = res["obs"]['+'][1:] + res["obs"]['-'][:-1]
		exp = res["exp"]['+'][1:] + res["exp"]['-'][:-1]

		for o, e in zip(obs, exp):
			try:
				h[e, o] += 1.0
			except IndexError:
				pass

	return h


def chunks_list(x, chunksize):
	"""Chunk a list into a list of lists
	"""
	
	n = max(1, chunksize)
	return [ x[i:i+n] for i in range(0, len(x), n) ]

def build_histogram_parallel_wrapper(intervals, read_filepath, fa_filepath, bm, half_window_width, size):
	"""	Helper function to create a new bamfile instance
	in the multicore/gridmap case
	"""
	
	#open a bamfile instance
	reads = cutcounts.bamfile(read_filepath)
	seq = pyfaidx.Fasta(fa_filepath)
	#canonical function
	return build_histogram(reads, seq, intervals, bm, half_window_width, size)

def build_histogram_multicore(read_filepath, fa_filepath, intervals, bm, half_window_width = 5, size = (200, 1000), processes = 8, chunksize = 250):
	"""Multiprocessing parallel implementation
	"""
	
	# make a partial function to wrap arguments
	from functools import partial
	wrapper_func = partial(build_histogram_parallel_wrapper, read_filepath = read_filepath, fa_filepath = fa_filepath, bm = bm, half_window_width = half_window_width, size = size)

	# chunk the intervals
	chunks = chunks_list(intervals, chunksize)

	# run a pool of threads
	import multiprocessing

	pool = multiprocessing.Pool(processes = processes)
	res_chunks = pool.map(wrapper_func, chunks)

	# combine the chunks
	res = res_chunks[0]
	for i in np.arange(1, len(res_chunks)):
		res += res_chunks[i]

	return res

def build_histogram_gridmap(read_filepath, intervals, bm, half_window_width = 5, size = (200, 1000), processes = 8, chunksize = 250):
	"""Gridmap parallel implementation
	"""

	# make a partial function to wrap arguments
	from functools import partial
	wrapper_func = partial(build_histogram_parallel_wrapper, read_filepath = read_filepath, bm = bm, half_window_width = half_window_width, size = size)

	# chunk the intervals
	chunks = chunks_list(intervals, chunksize)

	# run a pool of threads
	import gridmap

	joblist = [ gridmap.Job(wrapper_func, [chunk]) for chunk in chunks ]
	res_chunks = gridmap.process_jobs(joblist)

	# combine the chunks
	res = res_chunks[0]
	for i in np.arange(1, len(res_chunks)):
		res += res_chunks[i]

	return res
