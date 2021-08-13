"""
This module contains functions to perform multiple test correction.
"""
import numpy as np
import scipy.stats

import operator

from ..utils import bisect

def emperical_fdr(pvals_null, pvals):
    """Computes emperical FDR from a distribution of null (sampled) p-values

    Parameters
    ---------
    pvals_null : array_like
        Sampled p-values from null expectation
    pvals : array_like
        Observed p-values
    
    Returns
    -------
    out : ndarray
        Empirically adjusted p-values to be compared directly with the desired FDR
    """
    sorted_pvals_null = np.sort(np.ravel(pvals_null))
    sorted_pvals_idx = np.argsort(pvals)

    counts = bisect(sorted_pvals_null, pvals[sorted_pvals_idx])
    false_positive_rates = counts / len(sorted_pvals_null)
    false_positive_rates[false_positive_rates > 1] = 1
    return false_positive_rates[np.argsort(sorted_pvals_idx)]

# Not sure if Storey FDR is best (appears to be too conservative)

def pi0est(pvals, lamb = None):
    """
    
    Parameters
    ---------
    pvals : :class:`np.array'
        Observed p-values
    lambda : :class:`np.array`, optional
        Lambda 

    Returns
    -------

    """

    n = len(pvals)
    lamb = np.arange(0.05, 1, 0.05) if lamb is None else lamb
    pi0 = np.array([np.mean(pvals>=l)/(1-l) for l in lamb])

    #bootstrap method
    minpi0 = np.percentile(pi0, q = 10)
    W = np.array([np.sum(pvals >= l) for l in lamb])
    mse = (W / (n**2 * (1-lamb)**2)) * (1 - W/n) + (pi0 - minpi0)**2

    return min(pi0[mse== np.min(mse)], 1)

def qvalue(pvals):
    """Storey q-values
    
    Parameters
    ----------

    Returns
    -------

    """

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
    """Return Benjamini-Hochberg FDR q-values corresponding to p-values C{pv}.

    This function implements an algorithm equivalent to L{bh_rejected} but
    yields a list of 'adjusted p-values', allowing for rejection decisions
    based on any given threshold.

    Parameters
    ----------
    pvals : array_like
        P-values from a multiple statistical test

    Raises
    ------
    ValueError

    Returns
    -------
    out : array_like
        Adjusted p-values to be compared directly with the desired FDR
    """
    m = len(pvals)
    rank, sorted_pvals = list(zip(*sorted(enumerate(pvals), None, operator.itemgetter(1))))
    if pvals[0] < 0 or pvals[-1] > 1:
        raise ValueError("P-values must be between 0 and 1")
    qvalues = np.zeros(m)
    mincoeff = sorted_pvals[-1]
    qvalues[rank[-1]] = mincoeff
    for j in range(m-2, -1, -1):
        coeff = m*sorted_pvals[j]/float(j+1)
        if coeff < mincoeff:
            mincoeff = coeff
        qvalues[rank[j]] = mincoeff
    return qvalues

