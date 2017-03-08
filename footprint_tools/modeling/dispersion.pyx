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

from scipy import optimize

from functools import partial

def piecewise_three(x, x0, y0, x1, k1, k2, k3):

	# x0,y0 : first breakpoint
	# x1 : second breakpoint
	# k1,k2,k3 : 3 slopes.

	y1=y0+ k2*(x1-x0) # for continuity

	return (
	(x<x0)              *   (y0 + k1*(x-x0))      +
	((x>=x0) & (x<x1))  *   (y0 + k2*(x-x0))      +
	(x>=x1)             *   (y1 + k3*(x-x1)))


cdef class dispersion_model:
	"""
	Dispersion model base class
	"""

	cdef np.ndarray _h
	cdef np.ndarray _p, _r

	#cdef data_type_t _mu_params_a, _mu_params_b
	#cdef data_type_t _r_params_a, _r_params_b

	cdef np.ndarray _mu_params, _r_params
	cdef public object _mu_func, _r_func

	def __init__(self):
		
		# Histograms
		self._h = None
		
		# Emperical fits
		self._p = None
		self._r = None

		# Regression parameters
		#self._mu_params_a = self._mu_params_b = 0.0
		#self._r_params_a = self._r_params_b = 0.0
		
		self._mu_params = self._r_params = None
		self._mu_func = self._r_func = None

		pass

	# Pickling function
	def __reduce__(self):
		x = {}
		#x['h'] = self.h
		#x['p'] = self.p
		#x['r'] = self.r

		x['mu_params'] = self.mu_params
		x['r_params'] = self.r_params
		return (dispersion_model, (), x)

	# Pickling function
	def __setstate__(self, x):
		#self.h = x['h']
		#self.p = x['p']
		#self.r = x['r']
		
		self.mu_params = x['mu_params']
		self.r_params = x['r_params']

	property h:
		def __get__(self):
			return self._h
		def __set__(self, x):
			self._h = x

	property p:
		def __get__(self):
			return self._p
		def __set__(self, x):
			self._p = x

	property r:
		def __get__(self):
			return self._r
		def __set__(self, x):
			self._r = x

	property mu_params:
		def __get__(self):
			#return self._mu_params_a, self._mu_params_b
			return self._mu_params
		def __set__(self, x):
			#self._mu_params_a = x[0]
			#self._mu_params_b = x[1]
			self._mu_params = x
			#self._mu_func = np.poly1d(self._mu_params)
			self._mu_func = partial(piecewise_three, x0 = x[0], y0 = x[1], x1 = x[2], k1 = x[3], k2 = x[4], k3 = x[5])

	property r_params:
		def __get__(self):
			#return self._r_params_a, self._r_params_b
			return self._r_params
		def __set__(self, x):
			#self._r_params_a = x[0]
			#self._r_params_b = x[1]
			self._r_params = x
			#self._r_func = np.poly1d(self._r_params)
			self._r_func = partial(piecewise_three, x0 = x[0], y0 = x[1], x1 = x[2], k1 = x[3], k2 = x[4], k3 = x[5])


	cpdef fit_mu(self, data_type_t x):
		"""Computes the mu term for the negative binomial.
		
		:param x: predicted values to be converted to observed means
		
		:returns mu: (float)
		"""

		#cdef double res = self._mu_params_a + self._mu_params_b * x
		cdef double res = self._mu_func(x)

		return res if res > 0.0 else 0.0

	cpdef fit_r(self, data_type_t x):
		"""Computes the dispersion term for the negative binomial.
		Note that the model parameters estimate the inverse.

		:param x: predicted values to be converted to dispersion parameter r

		:returns r: (float)
		"""

		#cdef double res = 1.0 / (self._r_params_a + self._r_params_b * x)
		cdef double res = 1.0/self._r_func(x)
		return res if res > 0.0 else 1e-6

	def __str__(self):
		"""Print model to string"""

		#res = "mu = %0.4f + %0.4fx\n" % (self.mu_params)
		#res += "r = %0.4f + %0.4fx" % (self.r_params)
		res = ""
		return res

	cpdef log_pmf_values(self, data_type_t [:] exp, data_type_t [:] obs):
		"""Computing the probability mass function"""

		cdef int i, n = exp.shape[0]
		cdef double r, mu
		cdef data_type_t [:] res = np.ones(n, dtype = np.float64, order = 'c')

		for i in range(n):
			r = self.fit_r(exp[i])
			mu = self.fit_mu(exp[i])
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
			r = self.fit_r(exp[i])
			mu = self.fit_mu(exp[i])
			res[i] = nbinom.cdf(<int>obs[i], r/(r+mu), r)

		return res

	cpdef resample_p_values(self, data_type_t [:] x, int times):

		cdef int i, j, n = x.shape[0]
		cdef data_type_t r, mu, k

		cdef long [:] vals

		cdef data_type_t [:,:] res_pvals = np.ones((n, times), dtype = np.float64, order = 'c')
		cdef long [:,:] res_vals = np.zeros((n, times), dtype = np.int_, order = 'c')

		for i in range(n):
			r = self.fit_r(x[i])
			mu = self.fit_mu(x[i])
		
			vals = np.random.negative_binomial(r, r/(r+mu), times)
			
			for j in range(times):
				res_pvals[i, j] = nbinom.cdf(<int>vals[j], r/(r+mu), r)
				res_vals[i, j] = vals[j]

		return res_vals, res_pvals


@cython.boundscheck(True)
@cython.wraparound(True)

def learn_dispersion_model(h, cutoff = 250, trim = [2.5, 97.5]):
	"""Learn a dispersion model from the
	expected vs. observed histogram
	"""

	size = h.shape[0]
	n = np.zeros(size)
	p = np.zeros(size)
	r = np.zeros(size)

	#thresholded = []

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
			
			#thresholded.append(i)
			mu = np.mean(x[lower:upper])
			var = np.var(x[lower:upper])

			est_r = (mu * mu) / (var - mu)
			if est_r <= 0.0: est_r = 1.0

			est_p = est_r / (est_r + mu)
			
			(p[i], r[i]) = nbinom.fit(x[lower:upper], p = est_p, r = est_r)

		else:

			p[i] = r[i] = np.nan

	# Back-compute the mean values from the negative binomial parameters
	mus = p*r/(1-p)

	# Extrapolate using polynomial fit
	x = np.arange(size)
	sele = np.isfinite(mus)

	#w = np.log10(n)

	# Extrapolate the fit
	#import statsmodels.api as sm
	
	#xx = sm.add_constant(x)
	
	# Kind of a 'hacked' solution to fit the dispersion values; remove the first 5 thresholded datapoints

	#model_mu = sm.WLS(mus[sele], xx[sele], w[sele]).fit()
	#model_mu = np.polyfit(x[valid], mus[valid], deg = 1, w = w[valid])p0=(2,2,35, 1, 1, 1)
	
	
	x0 = 5
	y0 = mus[x0]
	k1 = (y0 - mus[0]) / (x0)

	x1 = 10
	y1 = mus[x1]
	k2 = (y1 - y0) / (x1 - x0)

	x2 = x[sele][-1]
	y2 = mus[sele][-1]
	k3 = (y2 - y1) / (x2 - x1)

	model_mu, e = optimize.curve_fit(piecewise_three, x[sele], mus[sele], p0 = (x0, y0, x1, k1, k2, k3))

	# remove non-linear part of r parameters
	# sele[:10] = False

	#model_r = sm.WLS((1/r)[sele], xx[sele], w[sele]).fit()
	#model_r = np.polyfit(x[valid], 1.0/r[valid], deg = 1, w = w[valid])
	
	
	x0 = 2
	y0 = r[x0]
	k1 = (1/y0 - 1/r[0]) / (x0)

	x1 = 5
	y1 = r[x1]
	k2 = (1/y1 - 1/y0) / (x1 - x0)

	x2 = x[sele][-1]
	y2 = r[sele][-1]
	k3 = (1/y2 - 1/y1) / (x2 - x1)

	model_r, e = optimize.curve_fit(piecewise_three, x[sele], 1.0/r[sele], p0 = (x0, y0, x1, k1, k2, k3))

	# Create a dispersion model class
	res = dispersion_model()
	res.h = h
	res.p = p
	res.r = r

	#res.mu_params = model_mu.params[::-1]
	res.mu_params = model_mu
	#res.r_params = model_r.params[::-1]
	res.r_params = model_r
	
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
		
		#model.mu_params = params["mu_params"]
		#model.r_params = params["r_params"]

		model.mu_params = base64decode(params["mu_params"])
		model.r_params = base64decode(params["r_params"])

		if "h" in params:
			model.h = base64decode(params["h"])
		if "p" in params:
			model.p = base64decode(params["p"])
		if "r" in params:
			model.r = base64decode(params["r"])
		
	return model

def write_dispersion_model(model):

	import json

	#out = { "mu_params": [model.mu_params[0], model.mu_params[1]], 
	#		"r_params": [model.r_params[0], model.r_params[1]],
	#		"h": base64encode(np.asarray(model.h)),
	#		"p": base64encode(np.asarray(model.p)),
	#		"r": base64encode(np.asarray(model.r))
	#	}

	out = { "mu_params": base64encode(np.asarray(model.mu_params, order = 'C')), 
			"r_params": base64encode(np.asarray(model.r_params, order = 'C')),
			"h": base64encode(np.asarray(model.h, order = 'C')),
			"p": base64encode(np.asarray(model.p, order = 'C')),
			"r": base64encode(np.asarray(model.r, order = 'C'))
		}

	return json.dumps(out, indent = 4)

