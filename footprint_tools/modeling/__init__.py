__all__ = ["bias", "dispersion"]

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
import scipy.stats

import pyfaidx

def reverse_complement(seq):

		compl = { 'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'}

		res = ''.join([ compl[base] for base in seq ])

		return res[::-1]

'''
    0  1  2  3  4  5
5'-[]-[]-[]-[]-[]-[]-3'
3'-[]-[]-[]-[]-[]-[]-5'
    5  4  3  2  1  0
'''

def predict_interval(reads, faidx, interval, bm, half_window_width = 10, smoothing_class = None):

	obs_counts = {'+': None, '-': None}
	exp_counts = {'+': None, '-': None}
	win_counts = {'+': None, '-': None}

	# get sequence (adding extra tails depending on the bias model)
	seq = faidx[interval.chrom][interval.start-half_window_width-bm.offset:interval.end+half_window_width+bm.offset].seq.upper()

	# 
	padding = half_window_width
	if smoothing_class:
		padding = padding + smoothing_class.half_window_width

	counts = reads[interval.widen(padding)]

	predict_func = lambda x: bm.predict(x[1:], x[0])[half_window_width]

	for strand in ['+', '-']:

		# Pre-calculate the sequence bias probability table
		if strand == '+':
			probs = bm.probs(seq)#[bm.offset:-bm.offset]
		else:
			probs = bm.probs(reverse_complement(seq))[::-1]

		#print len(interval), len(seq), len(probs)
		#print probs

		# get indicies for windowns
		idx = np.vstack( np.arange(i-half_window_width, i+half_window_width+1) for i in np.arange(half_window_width, len(counts[strand])-half_window_width) )

		# count tags in window
		w_counts = np.apply_along_axis(np.sum, 1, counts[strand][idx])

		# smooth windows if necesary
		if smoothing_class:
			
			sm_w_counts = smoothing_class.smooth(w_counts, half_window_width)
			win_counts[strand] = sm_w_counts

			z = np.concatenate((sm_w_counts[:, np.newaxis], probs[idx[0:len(interval),:]]), axis = 1)

		else:

			win_counts[strand] = w_counts
			z = np.concatenate((w_counts[:, np.newaxis], probs[idx[0:len(interval),:]]), axis = 1)

		obs_counts[strand] = counts[strand][padding:-padding]
		exp_counts[strand] = np.apply_along_axis(predict_func, 1, z)

	return { "obs": obs_counts, "exp": exp_counts, "smoothed": win_counts }



