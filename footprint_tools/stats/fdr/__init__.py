from bisect import bisect

import numpy as np

def sample_null(arr, model, windowed_pval_func, times):
	"""Sample data from the expected distribution for emperical FDR

	Parameters
	----------
		arr (array): 
			Expected counts
	Returns
	-------
	"""

	# Sample counts
	sampled_counts = np.vstack([ model.resample(arr) for i in range(times)])
	
	# Compute P-values
	pvals = np.apply_along_axis(lambda x: model.p_value(arr, x), 1, sampled_counts)
	
	# Compute windowed pvals
	w_pvals = np.apply_along_axis(windowed_pval_func, 1, pvals)

	return w_pvals

def emperical_fpr(arr, pvals):
	"""
	"""
	sorted_arr = np.sort(arr, axis = 1)
	sorted_pvals_idx = np.argsort(pvals)

	n = np.arange(1, pvals.shape[0]+1, dtype = np.float64)

	counts = np.apply_along_axis(lambda x: bisect(x, pvals[sorted_pvals_idx]), 1, sorted_arr)
	false_positive_rates = np.median(counts, axis = 0) / n

	return false_positive_rates[np.argsort(sorted_pvals_idx)]
