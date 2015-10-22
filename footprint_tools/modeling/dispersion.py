import numpy as np
import scipy.stats

from footprint_tools import cutcounts, stats, modeling

import pyfaidx

'''
Dispersion model base class

'''

class dispersion_model:

	def __init__(self):

		self.h = None
		
		self.mu_params = None
		self.r_params = None

	def mu(self, x):

		res = self.mu_params[0] + self.mu_params[1] * x
		return res if res > 0.1 else 0.1

	def r(self, x):

		return (self.r_params[0] + (self.r_params[1] * x))

	def __str__(self):

		res = "mu: y = %0.4f + %0.4fx\n" % (self.mu_params[0], self.mu_params[1])
		res += "r: y = %0.4f + %0.4fx" % (self.r_params[0], self.r_params[1])
		return res

	def p_value(self, exp, obs):

		mu = self.mu(exp)
		r = self.r(exp)
		p = r/(r+mu)

		#print exp, obs, r, p, mu
		return stats.nbinom.cdf(obs, 1-p, r)

'''
 Functions to create histogram of expected vs. observed cleavages	
'''

'''
Creates a histogram of the observed data 
at each of the the epxected counts within
a set of intervals

Returns a (m,n) matrix that holds counts
Used to fit the dispersion model
'''

def build_histogram(reads, faidx, intervals, bm, half_window_width = 5, size = (200, 1000)):

	h = np.zeros(size)

	for interval in intervals:

		res = modeling.predict_interval(reads, faidx, interval, bm, half_window_width)

		obs = res["obs"]['+'][1:] + res["obs"]['-'][:-1]
		exp = res["exp"]['+'][1:] + res["exp"]['-'][:-1]

		for o, e in zip(obs, exp):
			try:
				h[e, o] += 1.0
			except IndexError:
				pass

	return h

'''
Helper function to create a new bamfile instance in the multicore/gridmap case
'''

def build_histogram_parallel_wrapper(intervals, read_filepath, fa_filepath, bm, half_window_width, size):
	
	#open a bamfile instance
	reads = cutcounts.bamfile(read_filepath)
	faidx = pyfaidx.Fasta(fa_filepath)
	#canonical function
	return build_histogram(reads, faidx, intervals, bm, half_window_width, size)

'''
Helper function to chunk the intervals
'''

def chunks_list(x, chunksize):
	n = max(1, chunksize)
	return [ x[i:i+n] for i in range(0, len(x), n) ]


'''
Multiprocessing parallel implementation
'''

def build_histogram_multicore(read_filepath, fa_filepath, intervals, bm, half_window_width = 5, size = (200, 1000), processes = 8, chunksize = 250):

	# make a partial function to wrap arguments
	from functools import partial
	wrapper_func = partial(build_histogram_parallel_wrapper, read_filepath = read_filepath, fa_filepath = fa_filepath, bm = bm, half_window_width = half_window_width, size = size)

	# chunk the intervals
	chunks = chunks_list(intervals, chunksize)

	# run a pool of threads
	import multiprocessing

	pool = multiprocessing.Pool(processes = processes)
	res_chunks = pool.map(wrapper_func, chunks)

	# combine the chunks
	res = res_chunks[0]
	for i in np.arange(1, len(res_chunks)):
		res += res_chunks[i]

	return res

'''
Gridmap parallel implementation
'''

def build_histogram_gridmap(read_filepath, intervals, bm, half_window_width = 5, size = (200, 1000), processes = 8, chunksize = 250):

	# make a partial function to wrap arguments
	from functools import partial
	wrapper_func = partial(build_histogram_parallel_wrapper, read_filepath = read_filepath, bm = bm, half_window_width = half_window_width, size = size)

	# chunk the intervals
	chunks = chunks_list(intervals, chunksize)

	# run a pool of threads
	import gridmap

	joblist = [ gridmap.Job(wrapper_func, [chunk]) for chunk in chunks ]
	res_chunks = gridmap.process_jobs(joblist)

	# combine the chunks
	res = res_chunks[0]
	for i in np.arange(1, len(res_chunks)):
		res += res_chunks[i]

	return res

'''
Function to learn the dispersion model from the expected vs. observed histogram
'''

def learn_dispersion_model(h, cutoff = 100, trim = [2.5, 97.5]):

	size = h.shape[0]

	n = np.zeros(size)
	mu = np.zeros(size)
	var = np.zeros(size)
		
	p = np.zeros(size)
	r = np.zeros(size)

	thresholded = []

	# Make an initial negative binomial fit

	for i in np.arange(size):
		
		# Un-furl H array

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

			mu[i] = np.mean(x[lower:upper])
			var[i] = np.var(x[lower:upper])

			est_r = (mu[i] * mu[i]) / (var[i] - mu[i])
			est_p = est_r / (est_r + mu[i])

			(p[i], r[i]) = stats.nbinom.fit(x[lower:upper], p = est_p, r = est_r)

	# Extrapolate the fit

	import statsmodels.api as sm

	x = sm.add_constant(np.arange(size))
	model_mu = sm.WLS(mu[thresholded], x[thresholded], n[thresholded]).fit()
	model_r = sm.OLS(r[thresholded], x[thresholded]).fit()

	# Create a dispersion model class

	res = dispersion_model()
	res.h = h

	res.mu_params = model_mu.params
	res.r_params = model_r.params

	# Return

	return res


# Read and write functions for making a portable dispersion model

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

	res = dispersion_model()

	with open(filename) as f:
	
		params = json.load(f)

		res.mu_params = params["mu_params"]
		res.r_params = params["r_params"]
		res.h = base64decode(params["h"])

	return res

def write_dispersion_model(model):

	import json

	out = { "mu_params": [model.mu_params[0], model.mu_params[1]], 
			"r_params": [model.r_params[0], model.r_params[1]],
			"h": base64encode(model.h) }

	return json.dumps(out, indent = 4)

