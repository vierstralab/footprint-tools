# Copyright 2015 Jeff Vierstra

# cython: boundscheck=False
# cython: wraparound=False
# cython: nonecheck=False

from ..stats.distributions import nbinom
import numpy as np

cimport numpy as np
cimport cython

from ..stats.distributions cimport nbinom

ctypedef np.float64_t data_type_t

cdef class dispersion_model:
	"""
	Dispersion model base class
	"""

	cdef np.ndarray _h
	cdef data_type_t _mu_params_a, _mu_params_b
	cdef data_type_t _r_params_a, _r_params_b

	def __init__(self):
		self._h = None		
		self._mu_params_a = self._mu_params_b = 0.0
		self._r_params_a = self._r_params_b = 0.0
		pass

	# Pickling function
	def __reduce__(self):
		x = {}
		#x['h'] = self.h
		x['mu_params'] = self.mu_params
		x['r_params'] = self.r_params
		return (dispersion_model, (), x)

	# Pickling function
	def __setstate__(self, x):
		#self.h = x['h']
		self.mu_params = x['mu_params']
		self.r_params = x['r_params']

	property h:
		def __get__(self):
			return self._h
		def __set__(self, x):
			self._h = x

	property mu_params:
		def __get__(self):
			return self._mu_params_a, self._mu_params_b
		def __set__(self, x):
			self._mu_params_a = x[0]
			self._mu_params_b = x[1]

	property r_params:
		def __get__(self):
			return self._r_params_a, self._r_params_b
		def __set__(self, x):
			self._r_params_a = x[0]
			self._r_params_b = x[1]

	cpdef mu(self, data_type_t x):
		"""Computes the mu term for the negative binomial.
		
		:param x: predicted values to be converted to observed means
		
		:returns mu: (float)
		"""

		cdef double res = self._mu_params_a + self._mu_params_b * x
		return res if res > 0.0 else 0.0

	cpdef r(self, data_type_t x):
		"""Computes the dispersion term for the negative binomial.
		Note that the model parameters estimate the inverse.

		:param x: predicted values to be converted to dispersion parameter r

		:returns r: (float)
		"""

		cdef double res = 1.0 / (self._r_params_a + self._r_params_b * x)
		return res if res > 0.0 else 1e-6

	def __str__(self):
		"""Print model to string"""

		res = "mu = %0.4f + %0.4fx\n" % (self.mu_params)
		res += "r = %0.4f + %0.4fx" % (self.r_params)
		return res

	cpdef log_pmf_values(self, data_type_t [:] exp, data_type_t [:] obs):
		"""Computing the probability mass function"""

		cdef int i, n = exp.shape[0]
		cdef double r, mu
		cdef data_type_t [:] res = np.ones(n, dtype = np.float64, order = 'c')

		for i in range(n):
			r = self.r(exp[i])
			mu = self.mu(exp[i])
			res[i] = nbinom.logpmf(<int>obs[i], (r/(r+mu)), r)

		return res

	cpdef p_values(self, data_type_t [:] exp, data_type_t [:] obs):
		"""Computes CDF p-value from negative binomial

		:param exp: 
		:param obs:

		:return p:
		"""
		
		cdef int i, n = exp.shape[0]
		cdef double r, mu
		cdef data_type_t [:] res = np.ones(n, dtype = np.float64, order = 'c')

		for i in range(n):
			r = self.r(exp[i])
			mu = self.mu(exp[i])
			res[i] = nbinom.cdf(<int>obs[i], r/(r+mu), r)

		return res

	cpdef resample_p_values(self, data_type_t [:] x, int times):

		cdef int i, j, n = x.shape[0]
		cdef data_type_t r, mu, k

		cdef long [:] vals

		cdef data_type_t [:,:] res_pvals = np.ones((n, times), dtype = np.float64, order = 'c')
		cdef long [:,:] res_vals = np.zeros((n, times), dtype = np.int_, order = 'c')

		for i in range(n):
			r = self.r(x[i])
			mu = self.mu(x[i])
		
			vals = np.random.negative_binomial(r, r/(r+mu), times)
			
			for j in range(times):
				res_pvals[i, j] = nbinom.cdf(<int>vals[j], r/(r+mu), r)
				res_vals[i, j] = vals[j]

		return res_vals, res_pvals

def learn_dispersion_model(h, cutoff = 100, trim = [2.5, 97.5]):
	"""Learn a dispersion model from the
	expected vs. observed histogram
	"""

	size = h.shape[0]
	n = np.zeros(size)
	p = np.zeros(size)
	r = np.zeros(size)

	thresholded = []

	# Make an initial negative binomial fit
	for i in np.arange(size):
	
		pos = 0
		x = np.zeros(np.sum(h[i,:]))
		for j in np.arange(len(h[i,:])):
			num = int(h[i, j])
			x[pos:pos+num] = j
			pos += num 

		lower = int(np.floor(pos*(trim[0]/100.0)))
		upper = int(np.ceil(pos*(trim[1]/100.0)))

		n[i] = len(x)

		if n[i] >= cutoff:
			
			thresholded.append(i)
			mu = np.mean(x[lower:upper])
			var = np.var(x[lower:upper])

			est_r = (mu * mu) / (var - mu)
			if est_r <= 0.0: est_r = 1.0

			est_p = est_r / (est_r + mu)
			
			(p[i], r[i]) = nbinom.fit(x[lower:upper], p = est_p, r = est_r)

	# Back-compute the mean values from the negative binomial parameters
	mus = p*r/(1-p)

	# Extrapolate the fit
	import statsmodels.api as sm

	x = sm.add_constant(np.arange(size))

	# 
	model_mu = sm.WLS(mus[thresholded], x[thresholded], n[thresholded]).fit()

	# Kind of a 'hacked' solution to fit the dispersion values; remove the first 5 thresholded datapoints
	model_r = sm.OLS(1/r[thresholded[5:]], x[thresholded[5:]]).fit()

	# Create a dispersion model class
	res = dispersion_model()
	res.h = h
	res.mu_params = model_mu.params
	res.r_params = model_r.params

	return res

# Read and write functions for making a portable dispersion model

import json
import base64

def base64encode(x):
	
	return [str(x.dtype), base64.b64encode(x), x.shape]

def base64decode(x):

	dtype = np.dtype(x[0])
	arr = np.frombuffer(base64.decodestring(x[1]), dtype)
	if len(x) > 2:
		return arr.reshape(x[2])
	return arr

def read_dispersion_model(filename):
	
	import json

	model = dispersion_model()
	with open(filename) as f:
		params = json.load(f)
		model.mu_params = params["mu_params"]
		model.r_params = params["r_params"]
		model.h = base64decode(params["h"])
		
	return model

def write_dispersion_model(model):

	import json

	out = { "mu_params": [model.mu_params[0], model.mu_params[1]], 
			"r_params": [model.r_params[0], model.r_params[1]],
			"h": base64encode(np.asarray(model.h)) }
	return json.dumps(out, indent = 4)

