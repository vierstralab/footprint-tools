# cython: boundscheck=False
# cython: wraparound=False
# cython: nonecheck=False

from .. import cutcounts, stats, modeling

import numpy as np

cimport cython
cimport numpy as np

ctypedef np.float64_t data_type_t

cdef class dispersion_model:
	"""
	Dispersion model base class
	"""

	cdef data_type_t _mu_params_a, _mu_params_b
	cdef data_type_t _r_params_a, _r_params_b

	def __init__(self):
		#self.h = None		
		self._mu_params_a = self._mu_params_b = 0.0
		self._r_params_a = self._r_params_b = 0.0
		pass

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

	cdef _mu(self, data_type_t [:] x):
		cdef int i, n = x.shape[0]
		cdef data_type_t [:] res = np.empty(n, dtype = np.float64)
		for i in range(n):
			res[i] = self._mu_params_a + self._mu_params_b * x[i]
		return res

	def mu(self, x):
		"""Computes the mu term for the negative binomial.
		
		Parameters
		----------
		
		Returns
		-------
		(float): mu
		"""
		return self._mu(x)

	cdef _r(self, data_type_t [:] x):
		cdef int i, n = x.shape[0]
		cdef data_type_t [:] res = np.empty(n, dtype = np.float64)
		for i in range(n):
			res[i] = 1.0 / (self._r_params_a + self._r_params_b * x[i])
		return res

	def r(self, x):
		
		"""Computes the dispersion term for the negative binomial.
		Note that the model parameters estimate the inverse.

		Parameters
		----------

		Returns
		-------
		r : float
			Dispersion term r
		"""
		return self._r(x)

	def __str__(self):
		
		res = "mu = %0.4f + %0.4fx\n" % (self.mu_params)
		res += "r = %0.4f + %0.4fx" % (self.r_params)
		return res

	cpdef p_value(self, data_type_t [:] exp, data_type_t [:] obs):
		
		"""Computes log p-value from negative binomial
		
		Parameters
		----------

		Returns
		-------

		"""
		cdef data_type_t [:] mu = self._mu(exp)
		cdef data_type_t [:] r = self._r(exp)
		
		cdef int i, n = exp.shape[0]
		cdef data_type_t [:] p = np.empty(n, dtype = np.float64)

		for i in range(n):
			p[i] = r[i]/(r[i]+mu[i])

		cdef np.ndarray[data_type_t, ndim = 1] pvals = np.ones(n, dtype = np.float64)
		cdef data_type_t [:] pvals_view = pvals

		for i in range(n):
			pvals_view[i] = stats.nbinom.cdf(obs[i], p[i], r[i])

		return pvals

	cpdef resample(self, data_type_t [:] x):
		
		"""Resamples cleavage counts from negative binomial
		
		Parameters
		----------

		Returns
		-------

		Todo
		----
		Replace numpy RVS with a home-made solution for speed

		"""
		cdef data_type_t [:] mu = self._mu(x)
		cdef data_type_t [:] r = self._r(x)
		
		cdef int i, n = x.shape[0]
		cdef data_type_t [:] p = np.empty(n, dtype = np.float64)
		for i in range(n):
			p[i] = r[i]/(r[i]+mu[i])

		cdef np.ndarray[data_type_t, ndim = 1] samp = np.zeros(n, dtype = np.float64)
		cdef data_type_t [:] samp_view = samp

		for i in range(n):
			samp_view[i] = stats.nbinom.rvs(p[i], r[i], 1)[0]
		
		return samp


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
			num = float(h[i, j])
			x[pos:pos+num] = j
			pos += num 

		lower = np.floor(pos*(trim[0]/100.0))
		upper = np.ceil(pos*(trim[1]/100.0))

		n[i] = len(x)

		if n[i] >= cutoff:
			thresholded.append(i)
			mu = np.mean(x[lower:upper])
			var = np.var(x[lower:upper])
			est_r = (mu * mu) / (var - mu)
			est_p = est_r / (est_r + mu)
			(p[i], r[i]) = stats.nbinom.fit(x[lower:upper], p = est_p, r = est_r)

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

def load_dispersion_model(filename):
	
	import json

	model = dispersion_model()
	with open(filename) as f:
		params = json.load(f)
		model.mu_params = params["mu_params"]
		model.r_params = params["r_params"]
		#model.h = base64decode(params["h"])
	return model

def write_dispersion_model(model):

	import json

	out = { "mu_params": [model.mu_params[0], model.mu_params[1]], 
			"r_params": [model.r_params[0], model.r_params[1]],
			"h": base64encode(model.h) }
	return json.dumps(out, indent = 4)

