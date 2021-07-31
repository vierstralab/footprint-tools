# Copyright 2015 Jeff Vierstra

# cython: boundscheck=False
# cython: wraparound=False
# cython: nonecheck=False
# cython: embedsignature=True

cimport numpy as np

import numpy as np

from genome_tools import genomic_interval, genomic_interval_set

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
	"""Computes reverse complement of a DNA sequence
	
	Parameters
	----------
	seq : str
	    DNA sequence string
 	
	Returns
	-------
	str
	    Reverse complement of ``seq``
	"""
	compl = { 'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'}
	return ''.join([ compl.get(base, 'N') for base in seq ])[::-1]

class prediction(object):

	"""Class that holds a wrapper function to 
		compute the expected cleavage counts
	
	Attributes
	----------
	bm : bias.bias_model
	    Instance of bias model class
	fasta_func : pysam.FastaFile
		Instance of FastaFile
	half_window_width : int
	    Window width to apply bias model (final windows size = 2W+1)
	padding : int
	    Padding applied to ``interval`` when retrieving per-nucleotide data
	read_func : cutcounts.bamfile
		Instance of counts reader
	smoothing_clip : float
	    Fraction of nucleotides to trim when computing smoothed mean
	smoothing_half_window_width : int
	    Desc
	"""
	
	def __init__(self, read_func, fasta_func, bm, half_window_width = 5, smoothing_half_window_width = 0, smoothing_clip = 0.01):
		"""
		
		Parameters
		----------
		read_func : TYPE
		    Description
		fasta_func : TYPE
		    Description
		bm : TYPE
		    Description
		half_window_width : int, optional
		    Description
		smoothing_half_window_width : int, optional
		    Description
		smoothing_clip : float, optional
		    Description
		"""
		self.read_func = read_func
		self.fasta_func = fasta_func
		self.bm = bm

		# Smoothing/windowing parameters
		self.half_window_width = half_window_width
		self.smoothing_half_window_width = smoothing_half_window_width
		self.smoothing_clip = smoothing_clip

		# Base padding
		self.padding = self.half_window_width + smoothing_half_window_width

		# self.interval = interval

		# pad_interval = interval.widen(self.padding)
		
		# Note: We clip the first base when recombining the positive 
		# and negative strand, so add an extra base upfront
		# pad_interval.start -= 1

		# self.counts = reads[pad_interval]
		# self.seq = fasta[pad_interval.chrom][pad_interval.start-bm.offset():pad_interval.end+bm.offset()].seq.upper()

	def compute(self, x):
		"""Summary
		
		Parameters
		----------
		x : genome_toools.genomic_interval,  genome_toools.genomic_interval_set, or list
		    Genomic region(s) to generate predicted cleavages
		
		Returns
		-------
		tuple (or list of tuples)
		    Observed, expected and windowed cleavage counts
		"""

		# if isinstance(x, (list, genomic_interval_set)):
		# 	print("LIST")
		# 	for i in x:
		# 		yield self.compute(i)
	
		# Note: We clip the first base when recombining the positive 
		# and negative strand, so add an extra base upfront
		pad_interval = x.widen(self.padding)
		pad_interval.start -= 1

		# Get the raw cleavage counts and FASTA sequence
		raw_counts = self.read_func[pad_interval]
		raw_seq = self.fasta_func.fetch(pad_interval.chrom, pad_interval.start-self.bm.offset(), pad_interval.end+self.bm.offset()).upper()

		obs_counts = {'+': None, '-': None}
		exp_counts = {'+': None, '-': None}
		win_counts = {'+': None, '-': None}

		for strand in ['+', '-']:

			# Pre-calculate the sequence bias propensity table from bias model
			if strand == '+':
				probs = self.bm.probs(raw_seq)
			else:
				probs = self.bm.probs(reverse_complement(raw_seq))[::-1]

			exp, win = predict(np.ascontiguousarray(raw_counts[strand]), np.ascontiguousarray(probs), self.half_window_width, self.smoothing_half_window_width, self.smoothing_clip)

			w = raw_counts[strand].shape[0] - self.padding

			obs_counts[strand] = raw_counts[strand][self.padding:w]
			exp_counts[strand] = exp[self.padding:w]
			win_counts[strand] = win[self.padding:w]

		return obs_counts, exp_counts, win_counts


