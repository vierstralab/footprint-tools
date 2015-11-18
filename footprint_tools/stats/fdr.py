import numpy as np

from .windowing import *
from .bisect import *

def sample(arr, model, times):
	"""
	"""
	sampled_counts = np.vstack([ model.resample(arr) for i in range(times)])
	
	# compute pvals
	pval_func = lambda x: model.p_value(arr, x)
	pvals = np.apply_along_axis(pval_func, 1, sampled_counts)

	# compute windowed pvals
	windowed_pval_func = lambda x: windowed_p_value(x, 3, stouffers_z)[1]
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

def segment(arr, threshold):

	pass