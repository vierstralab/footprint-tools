from .. import cutcounts, stats, modeling

import numpy as np
import pyfaidx

class dispersion_model:

	'''
	Dispersion model base class

	'''

	def __init__(self):
		self.h = None		
		self.mu_params = None
		self.r_params = None

	def mu(self, x):
		return self.mu_params[0] + self.mu_params[1] * x

	def r(self, x):
		return (self.r_params[0] + (self.r_params[1] * x))

	def __str__(self):
		res = "mu = %0.4f + %0.4fx\n" % (self.mu_params[0], self.mu_params[1])
		res += "r = %0.4f + %0.4fx" % (self.r_params[0], self.r_params[1])
		return res

	def p_value(self, exp, obs, greater = False):
		mu = self.mu(exp)
		r = self.r(exp)
		p = r/(r+mu)
		if greater:
			return 1.0 - stats.nbinom.cdf(obs-1, 1-p, r)
		else:
			return stats.nbinom.cdf(obs, 1-p, r)

def build_histogram(reads, seq, intervals, bm, half_window_width = 5, size = (200, 1000)):
	
	'''
 	Functions to create histogram of expected vs. observed cleavages	
	Creates a histogram of the observed data at each of the the expected
	counts within a set of intervals

	Returns:
	A (m,n) matrix that holds counts used to fit the dispersion model
	'''

	h = np.zeros(size)
	for interval in intervals:
		res = modeling.predict_interval(reads, seq, interval, bm, half_window_width)

		obs = res["obs"]['+'][1:] + res["obs"]['-'][:-1]
		exp = res["exp"]['+'][1:] + res["exp"]['-'][:-1]

		for o, e in zip(obs, exp):
			try:
				h[e, o] += 1.0
			except IndexError:
				pass

	return h

def build_histogram_parallel_wrapper(intervals, read_filepath, fa_filepath, bm, half_window_width, size):
	
	'''
	Helper function to create a new bamfile instance in the multicore/gridmap case
	'''
	
	#open a bamfile instance
	reads = cutcounts.bamfile(read_filepath)
	seq = pyfaidx.Fasta(fa_filepath)
	#canonical function
	return build_histogram(reads, seq, intervals, bm, half_window_width, size)

def chunks_list(x, chunksize):
	
	'''
	Helper function to chunk the intervals
	'''
	
	n = max(1, chunksize)
	return [ x[i:i+n] for i in range(0, len(x), n) ]


def build_histogram_multicore(read_filepath, fa_filepath, intervals, bm, half_window_width = 5, size = (200, 1000), processes = 8, chunksize = 250):
	'''
	Multiprocessing parallel implementation
	'''
	
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

def build_histogram_gridmap(read_filepath, intervals, bm, half_window_width = 5, size = (200, 1000), processes = 8, chunksize = 250):
	
	'''
	Gridmap parallel implementation
	'''

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

def learn_dispersion_model(h, cutoff = 100, trim = [2.5, 97.5]):
	
	'''
	Function to learn the dispersion model from the expected vs. observed histogram
	'''

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
	model_mu = sm.WLS(mus[thresholded], x[thresholded], n[thresholded]).fit()
	model_r = sm.OLS(r[thresholded], x[thresholded]).fit()

	# Create a dispersion model class
	res = dispersion_model()
	res.h = h
	res.mu_params = model_mu.params
	res.r_params = model_r.params

	return res

# Read and write functions for making a portable dispersion model
def base64encode(x):
	
	import base64
	
	return [str(x.dtype), base64.b64encode(x), x.shape]

def base64decode(x):
	
	import base64

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

