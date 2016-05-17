# Copyright 2016 Jeff Vierstra

# cython: boundscheck=False
# cython: wraparound=False
# cython: nonecheck=False

cimport cython
cimport numpy as np

from .distributions cimport nbinom

ctypedef np.float64_t data_type_t

import numpy as np
import statsmodels.api as sm

#from .distributions import nbinom

cdef extern from "cephes.h":
	double c_exp(double) nogil

def estimate_size_factors(counts):

	logcounts = np.log(counts)

	geo_mu = np.mean(logcounts, axis = 0)

	sf = np.apply_along_axis(lambda x: np.exp(np.median((x - geo_mu)[np.isfinite(geo_mu) & (x > 0)])), 1, logcounts)

	return sf

def base_mean(counts, size_factors):
	
	return np.mean(counts / size_factors[:, np.newaxis], axis = 0)

def base_variance(counts, size_factors):
	
	return np.var(counts / size_factors[:, np.newaxis], axis = 0)

def estimate_dispersion(mu, var, size_factors):

	xim = np.mean(1/size_factors)
	disp = (var - xim * mu) / (mu*mu)

	disp[np.isnan(disp)] = 1e-8
	disp[disp <= 0] = 1e-8

	return disp

def fit_dispersion(mu, disp):
	"""Fits the observed per-element dispersion to a global model using GLM regression"""

	u = np.where(disp > 1e-8)

	x = 1/mu[u]
	y = disp[u]

	x = sm.add_constant(x)
	fit = sm.GLM(y, x, family=sm.families.Gamma(sm.genmod.families.links.identity)).fit()

	return lambda x: fit.params[0] + fit.params[1]/x

cpdef data_type_t nbinom_exact_test(int kx, data_type_t px, data_type_t rx, int ky, data_type_t py, data_type_t ry, int lower_tail = True) nogil:

	cdef unsigned int i
	cdef int ks = kx + ky

	if ks == 0:
		return 1.

	#cdef data_type_t [:] ps = np.ones(ks, dtype = np.float64, order = 'c')

	cdef double d = 0.0
	cdef double n = 0.0
	cdef double pmf

	for i in range(ks):
		pmf = c_exp(nbinom.logpmf(i, px, rx) + nbinom.logpmf(ks-i, py, ry))
		if lower_tail and i < kx:
			n += pmf
		elif not lower_tail and i >= kx:
			n += pmf
		d += pmf

	#cdef data_type_t d = np.sum(ps)
	#cdef data_type_t n = np.sum(ps[:kx]) if lower_tail else np.sum(ps[kx:])
	cdef data_type_t p = n/d if d > 0 else 1.

	return min(2*p, 1.)

cpdef compute_pvalues(np.ndarray[data_type_t, ndim = 2] a, np.ndarray[data_type_t, ndim = 2] b, np.ndarray[data_type_t, ndim = 1] sfa, np.ndarray[data_type_t, ndim = 1] sfb, lower_tail = True):

	assert a.shape[1] == b.shape[1]

	cdef unsigned int l = a.shape[1]

	cdef np.ndarray[data_type_t, ndim = 1] ma = base_mean(a, sfa)
	cdef np.ndarray[data_type_t, ndim = 1] mb = base_mean(b, sfb)

	cdef np.ndarray[data_type_t, ndim = 1] va = base_variance(a, sfa)
	cdef np.ndarray[data_type_t, ndim = 1] vb = base_variance(b, sfb)

	# Estimate the dispersion
	cdef np.ndarray[data_type_t, ndim = 1] est_dispa = estimate_dispersion(ma, va, sfa)
	cdef np.ndarray[data_type_t, ndim = 1] est_dispb = estimate_dispersion(mb, vb, sfb)

	# Fit dispersion using a GLM (see above)
	cdef np.ndarray[data_type_t, ndim = 1] fit_dispa = fit_dispersion(ma, est_dispa)(ma)
	cdef np.ndarray[data_type_t, ndim = 1] fit_dispb = fit_dispersion(mb, est_dispb)(ma)

	# Take the maximum of the estimated or fit dispersion value
	cdef np.ndarray[data_type_t, ndim = 1] dispa = np.maximum(est_dispa, fit_dispa)
	cdef np.ndarray[data_type_t, ndim = 1] dispb = np.maximum(est_dispb, fit_dispb)

	#
	cdef np.ndarray[data_type_t, ndim = 1] mu = np.mean(np.vstack([a, b]) / np.concatenate([sfa, sfb])[:, np.newaxis], axis = 0)

	# Recompute scaled variance and dispersions
	cdef np.ndarray[data_type_t, ndim = 1] fva = (mu  * sum(sfa)) + (dispa * mu*mu * sum(sfa*sfa))
	cdef np.ndarray[data_type_t, ndim = 1] fvb = (mu  * sum(sfb)) + (dispa * mu*mu * sum(sfb*sfb))

	cdef np.ndarray[data_type_t, ndim = 1] fda = (fva - (mu*sum(sfa))) / (mu * sum(sfa))**2
	cdef np.ndarray[data_type_t, ndim = 1] fdb = (fvb - (mu*sum(sfb))) / (mu * sum(sfb))**2

	#
	cdef np.ndarray[data_type_t, ndim = 1] ka = np.sum(a, axis = 0)
	cdef np.ndarray[data_type_t, ndim = 1] kb = np.sum(b, axis = 0)

	cdef np.ndarray[data_type_t, ndim = 1] ra = 1/fda
	cdef np.ndarray[data_type_t, ndim = 1] rb = 1/fdb

	cdef np.ndarray[data_type_t, ndim = 1] pa = ra / (ra + mu * sum(sfa))
	cdef np.ndarray[data_type_t, ndim = 1] pb = rb / (rb + mu * sum(sfb))

	cdef np.ndarray[data_type_t, ndim = 1] p = np.ones(l, dtype = np.float64)
	cdef data_type_t[:] mvp = p

	cdef unsigned int i
	for i in range(l):
		mvp[i] = nbinom_exact_test(<int>ka[i], pa[i], ra[i], <int>kb[i], pb[i], rb[i], lower_tail)

	return p
		