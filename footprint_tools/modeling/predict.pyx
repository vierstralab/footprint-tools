# cython: boundscheck=False
# cython: wraparound=False
# cython: nonecheck=False

from libc.stdlib cimport free, malloc
cimport numpy as np

import numpy as np

cdef extern from "predict.h":
	struct result:
		double* exp
		double* win	
	void free_result_t(result*)
	result* fast_predict(double*, double*, int, int, int, double)

ctypedef np.float64_t data_type_t

cdef predict(data_type_t [:] obs, data_type_t [:] probs, int half_window_width, int smoothing_half_window_width, double smoothing_clip):

	cdef int i
	cdef int l = obs.shape[0]

	# Predict
	cdef result* res = fast_predict(&obs[0], &probs[0], l, half_window_width, smoothing_half_window_width, smoothing_clip)

	# Copy back into a python object
	cdef np.ndarray[data_type_t, ndim = 1, mode = 'c'] exp = np.zeros(l, dtype = np.float64, order = 'c')
	cdef data_type_t [:] exp_view = exp

	cdef np.ndarray[data_type_t, ndim = 1, mode = 'c'] win = np.zeros(l, dtype = np.float64, order = 'c')
	cdef data_type_t [:] win_view = win

	for i in range(l): 
		exp_view[i] = res.exp[i]
		win_view[i] = res.win[i]

	# Free up the memory from C
	free_result_t(res)

	return exp, win

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

def predict_interval(reads, fasta, interval, bm, half_window_width = 10, smoothing_half_window_width = 50, smoothing_clip = 0.05):
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
		Genomic interval to be computed
	bm : bias_model
		A bias model to be used for resampling data in windows

	half_window_width : int

	smoothing_class : smoothing class or None

	Returns
	-------
	A dictionary
	"""
	
	obs_counts = {'+': None, '-': None}
	exp_counts = {'+': None, '-': None}
	win_counts = {'+': None, '-': None}
	win_probs = {'+': None, '-': None}

	padding = half_window_width + smoothing_half_window_width
	# Pad the interval sequence by the half-windown width and the bm model offset
	seq = fasta[interval.chrom][interval.start-padding-bm.offset():interval.end+padding+bm.offset()].seq.upper()
	
	# Get counts
	counts = reads[interval.widen(padding)]

	for strand in ['+', '-']:

		# Pre-calculate the sequence bias propensity table from bias model
		if strand == '+':
			probs = bm.probs(seq)
		else:
			probs = bm.probs(reverse_complement(seq))[::-1]

		exp, win = predict(np.ascontiguousarray(counts[strand]), np.ascontiguousarray(probs), half_window_width, smoothing_half_window_width, smoothing_clip)

		w = counts[strand].shape[0] - padding

		obs_counts[strand] = counts[strand][padding:w]
		exp_counts[strand] = exp[padding:w]
		win_counts[strand] = win[padding:w]

	return { "obs": obs_counts, "exp": exp_counts, "win": win_counts }



