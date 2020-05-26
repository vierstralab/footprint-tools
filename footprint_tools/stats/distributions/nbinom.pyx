# Copyright 2015 Jeff Vierstra

# cython: embedsignature=True


import scipy.special
import scipy.optimize
import numpy as np

cimport cython
cimport numpy as np

cdef extern from "cephes.h":
	double c_exp(double) nogil
	double c_log(double) nogil
	double c_log1p(double) nogil
	double c_lgamma(double) nogil
	double c_incbet(double, double, double) nogil


def mle(par, data, sm):

	"""Objective function for MLE estimate according to
	https://en.wikipedia.org/wiki/Negative_binomial_distribution#Maximum_likelihood_estimation
	
	Parameters
	----------
	par : TYPE
	    Description
	data : TYPE
	    Description
	sm : TYPE
	    Description
	data: the points to be fit
	sm: \sum data / len(data)
	
	Returns
	-------
	TYPE
	    Description
	"""

	p = par[0]
	r = par[1]
	n = len(data)
	f0 = sm / (r+sm) - p
	f1 = np.sum( scipy.special.psi(data+r) ) - n * scipy.special.psi(r) + n * np.log( r / (r+sm) )
	
	return np.array([f0, f1])


def fit(data, p = None, r = None):
	"""Summary
	
	Parameters
	----------
	data : TYPE
	    Description
	p : None, optional
	    Description
	r : None, optional
	    Description
	
	Returns
	-------
	TYPE
	    Description
	"""
	if p is None or r is None:
		av = np.average(data)
		va = np.var(data)
		r = (av*av) / (va-av)
		p = (va-av) / (va)

	sm = np.sum(data)/len(data)
	x = scipy.optimize.fsolve(mle, np.array([p, r]), args = (data, sm))

	return (x[0], x[1])

cpdef data_type_t logpmf(int k, data_type_t p, data_type_t r) nogil:
	cdef double coeff = c_lgamma(k+r) - c_lgamma(k+1) - c_lgamma(r)
	return coeff + r * c_log(p) + k * c_log1p(-p)


cpdef data_type_t pmf(int k, data_type_t p, data_type_t r) nogil:
	return c_exp( logpmf(k, p, r) )

cpdef data_type_t cdf(int k, data_type_t p, data_type_t r) nogil:
	return c_incbet(r, k+1, p)

cpdef data_type_t mean(data_type_t p, data_type_t r) nogil:
	return p*r/(1-p)

"""Not yet implemented
def rvs(data_type_t p, data_type_t r, int n = 1):
	return np.random.negative_binomial(r, p, n)
"""