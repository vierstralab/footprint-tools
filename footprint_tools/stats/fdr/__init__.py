# Copyright 2015 Jeff Vierstra

from __future__ import division

from .bisect import *

import numpy as np

def emperical_fpr(pvals_null, pvals):
	"""
	"""
	sorted_pvals_null = np.sort(pvals_null, axis = 0)
	sorted_pvals_idx = np.argsort(pvals)

	n = np.arange(1, pvals.shape[0]+1, dtype = np.float64)

	counts = np.apply_along_axis(lambda x: bisect(x, np.asarray(pvals)[sorted_pvals_idx]), 0, sorted_pvals_null)
	false_positive_rates = np.median(counts, axis = 1) / n
	
	false_positive_rates[false_positive_rates > 1] = 1

	return false_positive_rates[np.argsort(sorted_pvals_idx)]
