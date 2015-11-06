import scipy.special
import scipy.optimize
import numpy as np
import mpmath as mp
 
class _nbinom(object):

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

	def pmf(self, k, p, r):
		try:
			p = np.float64(mp.gamma(k+r) / (mp.gamma(k+1) * mp.gamma(r)) * np.power(1-p, r) * np.power(p, k))
		except:
			p = 0.0
		return p

	def cdf(self, k, p, r):
		try:
			p = 1.0 - np.float64(mp.betainc(k+1, r, x1 = 0, x2 = p, regularized = True))
		except:
			p = 0.0
		return p

	def mean(self, p, r):
		return p*r/(1-p)

nbinom = _nbinom()
