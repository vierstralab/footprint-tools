import scipy.special
import scipy.optimize
import numpy as np
 
cimport cython
cimport numpy as np

ctypedef np.float64_t data_type_t

cdef extern from "cephes.h":
	double c_incbet(double, double, double)

cdef class _nbinom(object):

	def MLE(self, par, data, sm):
		
		'''
		Objective function for MLE estimate according to
		https://en.wikipedia.org/wiki/Negative_binomial_distribution#Maximum_likelihood_estimation
 
		Keywords:
		data -- the points to be fit
		sm -- \sum data / len(data)
		'''

		p = par[0]
		r = par[1]
		n = len(data)
		f0 = sm / (r+sm) - p
		f1 = np.sum( scipy.special.psi(data+r) ) - n * scipy.special.psi(r) + n * np.log( r / (r+sm) )
		
		return np.array([f0, f1])

	def fit(self, data, p = None, r = None):
		
		if p is None or r is None:
			av = np.average(data)
			va = np.var(data)
			r = (av*av)/(va-av)
			p = (va-av)/(va)

		sm = np.sum(data)/len(data)
		x = scipy.optimize.fsolve(self.MLE, np.array([p, r]), args=(data, sm))

		return (x[0], x[1])

	cpdef pmf(self, int k, data_type_t p, data_type_t r):
		#return negbinom_pmf(k, p, r)
		return 0

	cpdef cdf(self, int k, data_type_t p, data_type_t r):
		return c_incbet(r, k+1, p)

	cpdef mean(self, data_type_t p, data_type_t r):
		return p*r/(1-p)

	cpdef rvs(self, data_type_t p, data_type_t r, int n = 1):
		return np.random.negative_binomial(r, p, n)

nbinom = _nbinom()
