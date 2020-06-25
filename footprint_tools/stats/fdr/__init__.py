# Copyright 2015 Jeff Vierstra

from .bisect import *

import numpy as np
import scipy.stats


def emperical_fdr(pvals_null, pvals):
	"""
	"""

	''' DEPRECATED
	sorted_pvals_null = np.sort(pvals_null, axis = 0)
	sorted_pvals_idx = np.argsort(pvals)

	n = np.arange(1, pvals.shape[0]+1, dtype = np.float64)

	counts = np.apply_along_axis(lambda x: bisect(x, np.asarray(pvals)[sorted_pvals_idx]), 0, sorted_pvals_null)
	false_positive_rates = np.mean(counts, axis = 1) / n
	
	false_positive_rates[false_positive_rates > 1] = 1

	return false_positive_rates[np.argsort(sorted_pvals_idx)]
	'''

	sorted_pvals_null = np.sort(np.ravel(pvals_null))
	sorted_pvals_idx = np.argsort(pvals)

	counts = bisect(sorted_pvals_null, pvals[sorted_pvals_idx])
	false_positive_rates = counts / len(sorted_pvals_null)
	false_positive_rates[false_positive_rates > 1] = 1
	return false_positive_rates[np.argsort(sorted_pvals_idx)]

# Not sure if Storey FDR is best (appears to be too conservative)

def pi0est(pvals, lamb = None):
	""" """

	lamb = np.arange(0.05, 1, 0.05) if lamb is None else lamb
	
	n = len(pvals)

	pi0 = np.array([np.mean(pvals>=l)/(1-l) for l in lamb])

	#bootstrap method
	minpi0 = np.percentile(pi0, q = 10)
	W = np.array([np.sum(pvals >= l) for l in lamb])
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

def bh_qvalue(pvals):
    """
    Return Benjamini-Hochberg FDR q-values corresponding to p-values C{pv}.

    This function implements an algorithm equivalent to L{bh_rejected} but
    yields a list of 'adjusted p-values', allowing for rejection decisions
    based on any given threshold.

    @type pv: list
    @param pv: p-values from a multiple statistical test

    @rtype: list
    @return: adjusted p-values to be compared directly with the desired FDR
      level
    """
    m = len(pvals)
    rank, sorted_pvals = list(zip(*sorted(enumerate(pvals), None, operator.itemgetter(1))))
    if pvals[0] < 0 or pvals[-1] > 1:
        raise ValueError("p-values must be between 0 and 1")
    qvalues = np.zeros(m)
    mincoeff = sorted_pvals[-1]
    qvalues[rank[-1]] = mincoeff
    for j in range(m-2, -1, -1):
        coeff = m*sorted_pvals[j]/float(j+1)
        if coeff < mincoeff:
            mincoeff = coeff
        qvalues[rank[j]] = mincoeff
    return qvalues

