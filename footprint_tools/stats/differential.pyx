
# Copyright 2017 Jeff Vierstra

# cython: boundscheck=False
# cython: wraparound=False
# cython: nonecheck=False
# cython: embedsignature=True


import numpy as np
cimport numpy as np

import scipy.stats
import scipy.optimize

ctypedef np.float64_t data_type_t

from ..modeling cimport dispersion

ctypedef dispersion.dispersion_model dm_t

cdef extern from "cephes.h":
	double c_lgamma(double) nogil
	double c_log(double) nogil
	double c_sqrt(double) nogil
	double PI

cpdef np.ndarray[data_type_t, ndim = 3, mode = 'c'] compute_logpmf_values(dm, np.ndarray[data_type_t, ndim = 2, mode = 'c'] obs, np.ndarray[data_type_t, ndim = 2, mode = 'c'] exp,  data_type_t lo, data_type_t hi, int nslices):
	
	cdef int i 		
	cdef int n = obs.shape[0]
	cdef int m = obs.shape[1]

	bins = np.power(2, np.linspace(lo, hi, nslices))

	cdef np.ndarray[data_type_t, ndim = 3, mode = 'c'] e = np.ascontiguousarray(exp[:,:,np.newaxis] * bins)
	cdef np.ndarray[data_type_t, ndim = 3, mode = 'c'] o = np.ascontiguousarray(np.repeat(obs[:,:,np.newaxis], nslices, axis = 2))
	cdef np.ndarray[data_type_t, ndim = 3, mode = 'c'] res = np.zeros((n, m, nslices), order = 'c')

	for i in range(n):
		<dm_t>(dm[i]).log_pmf_values_0(np.ravel(e[i,:,:]), np.ravel(o[i,:,:]), np.ravel(res[i,:,:]))

	return res.T


#various stats functions (should be move to ..stats/distributions folder)

def invchi2_logpdf(x, nu, tau2):
	return (nu/2.0)*np.log(tau2*nu/2.0) + (-nu*tau2/2.0/x) - scipy.special.gammaln(nu/2.0) - (1.0+nu/2.0)*np.log(x)

def invchi2_likelihood(par, x):
	return -np.sum(invchi2_logpdf(x, par[0], par[1]))

cdef double t_logpmf(double x, double nu, double mu, double sig2) nogil:
	return c_lgamma((nu+1.0)/2.0) - c_lgamma(nu/2.0) - c_log(c_sqrt(PI*nu*sig2)) - ((nu+1.0)/2.0 * c_log(1.0+((x-mu)**2)/sig2/nu ))


cpdef np.ndarray[data_type_t, ndim = 2, mode = 'c'] compute_log_prior_t(np.ndarray[data_type_t, ndim = 2, mode = 'c'] obs_over_exp , data_type_t nu_0, data_type_t sig2_0,  data_type_t lo, data_type_t hi, int nslices):
	
	""" 
	Computes the log posterior distribution of the hyperparameters (normal * inv chi2 = students t)
	"""

	cdef int i, j
	cdef int n = obs_over_exp.shape[0]
	cdef int m = obs_over_exp.shape[1]

	cdef data_type_t mu, ssqdev, nu_1, sig2_1

	cdef np.ndarray[data_type_t, ndim = 2, mode = 'c'] res = np.zeros((m, nslices), dtype = np.float64, order = 'c')
	cdef data_type_t [:,:] res_view = res

	cdef double h = (hi-lo)/<double>nslices

	for i in range(m):
		
		mu = np.nanmean(obs_over_exp[:,i])
		ssqdev = np.nansum((obs_over_exp[:,i]-mu)**2)

		nu_1 = nu_0 + n
		sig2_1 = (nu_0*sig2_0 + ssqdev) / (nu_0 + n)

		with nogil:

			for j in range(0, nslices):
				#res_view[i, j] = scipy.stats.t.logpdf(bins[j], df = nu_1, loc = mu, scale = np.sqrt(sig2_1))
				res_view[i, j] = t_logpmf(lo+(j*h), nu_1, mu, sig2_1)

		res[i, :] += np.log(h)

	return res.T # - res.sum(axis = 1)[:, np.newaxis]).T




