# Copyright 2015 Jeff Vierstra

from __future__ import division

from .bisect import *

import numpy as np
import scipy.stats


def emperical_fdr(pvals_null, pvals):
	"""
	"""
	sorted_pvals_null = np.sort(pvals_null, axis = 0)
	sorted_pvals_idx = np.argsort(pvals)

	n = np.arange(1, pvals.shape[0]+1, dtype = np.float64)

	counts = np.apply_along_axis(lambda x: bisect(x, np.asarray(pvals)[sorted_pvals_idx]), 0, sorted_pvals_null)
	false_positive_rates = np.mean(counts, axis = 1) / n
	
	false_positive_rates[false_positive_rates > 1] = 1

	return false_positive_rates[np.argsort(sorted_pvals_idx)]

# Not sure if Storey FDR is best (appears to be too conservative)

def pi0est(pvals, lamb = None):
	""" """

	lamb = np.arange(0.05, 1, 0.05) if lamb is None else lamb
	
	n = len(pvals)

	pi0 = np.array(map(lambda l: np.mean(pvals>=l)/(1-l), lamb))

	#bootstrap method
	minpi0 = np.percentile(pi0, q = 10)
	W = np.array(map(lambda l: np.sum(pvals >= l), lamb))
	mse = (W / (n**2 * (1-lamb)**2)) * (1 - W/n) + (pi0 - minpi0)**2

	return min(pi0[mse== np.min(mse)], 1)
	#return lamb

def qvalue(pvals):
	""" """

	pi0 = pi0est(pvals)

	n = len(pvals)
	u = np.argsort(pvals)
	v = scipy.stats.rankdata(pvals, method = 'max')
	#qvals = (pi0*n*pvals) / v
	qvals = (pi0*n*pvals) / (v * (1 - (1-pvals)**n))
	qvals[u[n-1]] = min(qvals[u[n-1]], 1)
	for i in range(0, n-1)[::-1]:
		qvals[u[i]] = min(qvals[u[i]], qvals[u[i+1]])

	return qvals

