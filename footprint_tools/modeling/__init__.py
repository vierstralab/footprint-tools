__all__ = ["bias", "dispersion", "smoothing"]

'''
JDV Aug.2015

footprint toolkit/cleavage_model

Class used to perform various shufflings of cleavage data
in larger windows (e.g., hotspots). These functions aid in
generating shuffled data combining a local bias model and
a smoothing function. At some point I should put the data
smoothing functions into a separate class.
'''

import numpy as np

def reverse_complement(seq):
	"""Function that computes the reverse complement of a genomic sequence
	Args:
		seq (str): A DNA sequence
	Returns:
		str: A reverse complemented version of seq
	"""

	compl = { 'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'}
	return ''.join([ compl[base] for base in seq ])[::-1]

def predict_interval(reads, seq, interval, bm, half_window_width = 10, smoothing_class = None):
	"""Creates an expected distribution of cleavage counts within a genomic interval
		using observed data, a sequence preference model (bias), and a windowed smoothing
		function.
	Args:
		reads:
		seq:
		interval:
		bm:
		half_window_width (optional:
		smoothing_class (optional):
	Returns:
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
			#sm_w_counts = smoothing_class.smooth(w_counts, half_window_width)
			#win_counts[strand] = sm_w_counts
			#z = np.concatenate((sm_w_counts[:, np.newaxis], probs[idx[0:len(interval),:]]), axis = 1)
		else:
			win_counts[strand] = w_counts
			#z = np.concatenate((w_counts[:, np.newaxis], probs[idx[0:len(interval),:]]), axis = 1)

		z = np.concatenate((win_counts[strand][:, np.newaxis], probs[win_idx[0:len(interval), :]]), axis = 1)

		obs_counts[strand] = counts[strand][padding:-padding]
		exp_counts[strand] = np.apply_along_axis(predict_func, 1, z)

	return { "obs": obs_counts, "exp": exp_counts, "win": win_counts }

