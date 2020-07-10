# Copyright 2015 Jeff Vierstra

# cython: boundscheck=False
# cython: wraparound=False
# cython: nonecheck=False
# cython: embedsignature=True

from footprint_tools.stats.distributions import nbinom
import numpy as np
from scipy import optimize

import pwlf

cimport cython
cimport numpy as np

from footprint_tools.stats.distributions cimport nbinom


cpdef data_type_t piecewise_three(x, data_type_t x0, data_type_t x1, data_type_t x2,
	 								data_type_t y0, data_type_t y1, data_type_t y2,
	 								data_type_t k0, data_type_t k1, data_type_t k2):
	
	return (
        (x<x0) * (y0 + k0*(x)) +
        ((x>=x0) & (x<x1)) * (y1 + k1*(x)) +
        (x>=x1) * (y2 + k2*(x))
    )

cpdef data_type_t piecewise_four(x, data_type_t x0, data_type_t x1, data_type_t x2, data_type_t x3, 
	 								data_type_t y0, data_type_t y1, data_type_t y2, data_type_t y3,
	 								data_type_t k0, data_type_t k1, data_type_t k2, data_type_t k3):
	
	return (
        (x<x0) * (y0 + k0*(x)) +
        ((x>=x0) & (x<x1)) * (y1 + k1*(x)) +
        ((x>=x1) & (x<x2)) * (y2 + k2*(x)) +
        (x>=x2) * (y3 + k3*(x))
    )

cpdef data_type_t piecewise_five(x, data_type_t x0, data_type_t x1, data_type_t x2, data_type_t x3, data_type_t x4, 
	 								data_type_t y0, data_type_t y1, data_type_t y2, data_type_t y3, data_type_t y4,
	 								data_type_t k0, data_type_t k1, data_type_t k2, data_type_t k3, data_type_t k4):
	
	return (
        (x<x0) * (y0 + k0*(x)) +
        ((x>=x0) & (x<x1)) * (y1 + k1*(x)) +
        ((x>=x1) & (x<x2)) * (y2 + k2*(x)) +
        ((x>=x2) & (x<x3)) * (y3 + k3*(x)) +
        (x>=x3) * (y4 + k4*(x))
    )

cdef class dispersion_model:
	"""
	Dispersion model base class
	"""

	def __init__(self):
		
		# Histograms
		self._h = None
		
		# Emperical fits
		self._p = None
		self._r = None

		self._mu_params = self._r_params = None

		pass

	# Pickling function
	def __reduce__(self):
		x = {}
		x['mu_params'] = self.mu_params
		x['r_params'] = self.r_params
		return (dispersion_model, (), x)

	# Pickling function
	def __setstate__(self, x):
		self.mu_params = x['mu_params']
		self.r_params = x['r_params']

	property h:
		"""Histrogram of observed cleavages at each predicted cleavage rate"""
		def __get__(self):
			return self._h
		def __set__(self, x):
			self._h = x

	property p:
		"""Array of the negative binomial ML fit parameter `p`"""
		def __get__(self):
			return self._p
		def __set__(self, x):
			self._p = x

	property r:
		"""Array of the negative binomial ML fit parameter `r`"""
		def __get__(self):
			return self._r
		def __set__(self, x):
			self._r = x

	property mu_params:
		
		def __get__(self):
			return self._mu_params
		def __set__(self, x):
			self._mu_params = np.array(x, order = 'c')

	property r_params:
		def __get__(self):
			return self._r_params
		def __set__(self, x):
			self._r_params = np.array(x, order = 'c')


	cpdef data_type_t fit_mu(self, data_type_t x):
		"""Computes the mu term for the negative binomial.
		
		:param x: predicted values to be converted to observed means
		
		:returns mu: (float)
		"""

		cdef data_type_t [:] par = self._mu_params
		cdef data_type_t res = piecewise_three(x, *par)

		return res if res > 0.0 else 0.1

	cpdef data_type_t fit_r(self, data_type_t x):
		"""Computes the dispersion term for the negative binomial.
		Note that the model parameters estimate the inverse.

		:param x: predicted values to be converted to dispersion parameter r

		:returns: A fitted paramter for the negative bionomial distribution
		:rtype: float
		"""

		cdef data_type_t [:] par = self._r_params
		cdef data_type_t res = 1.0/piecewise_five(x, *par)

		return res if res > 0.0 else 1e-6

	def __str__(self):
		"""Print model to string"""
		res = "Oops...the functionality has not yet been implemented!"
		return res

	@cython.cdivision(True)
	cpdef data_type_t [:] log_pmf_values(self, data_type_t [:] exp, data_type_t [:] obs):
		"""Computing the probability mass function"""

		cdef int i, n = exp.shape[0]
		cdef double r, mu
		
		cdef data_type_t [:] res = np.zeros(n, dtype = np.float64, order = 'c')
	
		for i in range(n):
			r = self.fit_r(exp[i])
			mu = self.fit_mu(exp[i])
			res[i] = nbinom.logpmf(<int>obs[i], (r/(r+mu)), r)

		return res

	@cython.cdivision(True)
	cpdef data_type_t [:] pmf_values(self, data_type_t [:] exp, data_type_t [:] obs):
		"""Computing the probability mass function"""

		cdef int i, n = exp.shape[0]
		cdef double r, mu
		
		cdef data_type_t [:] res = np.zeros(n, dtype = np.float64, order = 'c')

		for i in range(n):
			r = self.fit_r(exp[i])
			mu = self.fit_mu(exp[i])
			res[i] = nbinom.pmf(<int>obs[i], (r/(r+mu)), r)

		return res


	@cython.cdivision(True)
	cpdef data_type_t [:] pmf_values_0(self, data_type_t [:] exp, data_type_t [:] obs, data_type_t [:] res):
		"""Computing the probability mass function (but writes to a matrix pointer)"""

		cdef int i, n = exp.shape[0]
		cdef double r, mu
		
		#cdef data_type_t [:] res = np.zeros(n, dtype = np.float64, order = 'c')

		for i in range(n):
			r = self.fit_r(exp[i])
			mu = self.fit_mu(exp[i])
			res[i] = nbinom.pmf(<int>obs[i], (r/(r+mu)), r)

		return res

	@cython.cdivision(True)
	cpdef data_type_t [:] log_pmf_values_0(self, data_type_t [:] exp, data_type_t [:] obs, data_type_t [:] res):
		"""Computing the probability mass function (but writes to a matrix pointer)"""

		cdef int i, n = exp.shape[0]
		cdef double r, mu
		
		#cdef data_type_t [:] res = np.zeros(n, dtype = np.float64, order = 'c')

		for i in range(n):
			r = self.fit_r(exp[i])
			mu = self.fit_mu(exp[i])
			res[i] = nbinom.logpmf(<int>obs[i], (r/(r+mu)), r)

		return res

	cpdef data_type_t [:] p_values(self, data_type_t [:] exp, data_type_t [:] obs):
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


def learn_dispersion_model(h, cutoff = 250, trim = (2.5, 97.5)):
	"""Learn a dispersion model from the expected vs. observed histogram
	
	:param h: 2-D :class:`numpy.array` containing the distribution of observerd cleavages at each expected cleavage rate
	:type h: numpy.array
	:param cutoff: Mininum number of observed cleavages tp perform ML negative binomial fit
	:type cutoff: int
	:param trim: Percent of data to trim  from the observed cleavage count (to mitigate outlier effects)
	:type trim: tuple of int

	:return: A dispersion model class :class:`dispersion_model`
	:rtype: dispersion_model

	:todo: Add exceptions for failure to fit, etc.
	"""

	size = int(h.shape[0])
	p = np.zeros(size)
	r = np.zeros(size)

	# Make an initial negative binomial fit
	for i in range(0, size):
	
		pos = 0
		x = np.zeros(int(np.sum(h[i,:])))
		for j in range(len(h[i,:])):
			num = int(h[i, j])
			x[pos:pos+num] = j
			pos += num

		# If more than 500k points downsample to
		# make curve-fitting tractable
		if len(x)>1e5:
			x=np.random.choice(x, size=int(1e5))
			x=np.sort(x)

		if len(x) >= cutoff:		

			lower = int(np.floor(x.shape[0]*(trim[0]/100.0)))
			upper = int(np.ceil(x.shape[0]*(trim[1]/100.0)))

			mu = np.mean(x[lower:upper])
			var = np.var(x[lower:upper])

			est_r = (mu * mu) / (var - mu)
			if est_r <= 0.0: 
				est_r = 10.0
			est_p = est_r / (est_r + mu)
			
			(p[i], r[i]) = nbinom.fit(x[lower:upper], p = est_p, r = est_r)

		else:
			p[i] = r[i] = np.nan

	# Back-compute the mean values from the negative binomial parameters
	mus = p*r/(1-p)

	# hack
	r[r>200]=200.

	# Extrapolate using polynomial fit
	x = np.arange(size)
	sele = np.isfinite(mus)

	first_x = np.min(x[sele])
	last_x = np.max(x[sele])*0.75

	# fit mu with a 3 segments
	fit_mu=pwlf.PiecewiseLinFit(x[sele], mus[sele])
	res=fit_mu.fit_with_breaks_force_points(np.linspace(first_x, last_x, 4), [0], [mus[0]])

	# fit r with 5 s segments
	fit_r=pwlf.PiecewiseLinFit(x[sele], 1.0/r[sele])
	res = optimize.minimize(fit_r.fit_with_breaks_opt, [3.0, 7.0, 15.0, 25.0])
	
	x0=np.zeros(6)
	x0[0]=first_x
	x0[-1]=last_x
	x0[1:-1]=res.x

	#res=fit_r.fit_with_breaks_force_points([first_x, 2, 7, 15, 25, last_x], [1], [1.0/r[1]])
	res=fit_r.fit_with_breaks_force_points(x0, [1], [1.0/r[1]])


	# Create a dispersion model class
	res = dispersion_model()
	res.h = h
	res.p = p
	res.r = r

	res.mu_params = list(fit_mu.fit_breaks[1:]) + list(fit_mu.intercepts) + list(fit_mu.slopes)
	res.r_params = list(fit_r.fit_breaks[1:]) + list(fit_r.intercepts) + list(fit_r.slopes)
	
	return res

# Read and write functions for making a portable dispersion model

import json
import base64

def base64encode(x):
	
	return [str(x.dtype), base64.b64encode(x), x.shape]

def base64decode(x):

	dtype = np.dtype(x[0])
	arr = np.frombuffer(base64.b64decode(x[1]), dtype)
	if len(x) > 2:
		return arr.reshape(x[2])
	return arr

def read_dispersion_model(filename):
	
	import json
	import urllib.request as request

	if filename.startswith('http'):
		file = request.urlopen(filename)
	else:
		file = open(filename, 'r')

	params = json.load(file)

	file.close()

	model = dispersion_model()
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
	"""Write a JSON format dispersion model

	:param model: An instance of :class:`dispersion_model`
	:return: JSON-formatted dump of dispersion model
	:rtype: str
	"""

	import json

	out = { "mu_params": base64encode(np.asarray(model.mu_params, order = 'C')), 
			"r_params": base64encode(np.asarray(model.r_params, order = 'C')),
			"h": base64encode(np.asarray(model.h, order = 'C')),
			"p": base64encode(np.asarray(model.p, order = 'C')),
			"r": base64encode(np.asarray(model.r, order = 'C'))
		}

	return json.dumps(out, indent = 4)

