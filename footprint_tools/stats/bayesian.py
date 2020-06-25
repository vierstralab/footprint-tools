# Copyright 2015 Jeff Vierstra

import numpy as np
import scipy.stats

from footprint_tools.stats import windowing

def compute_prior_weighted(fdr, w, cutoff = 0.05, pseudo = 0.5):
	"""Returns prior of whether a TF is bound or not 
	
	Parameters
	----------
	fdr : numpy.array
	    2D :class:`np.array` with FPR values for each datset and position
	w : numpy.array
	    Binary array whether samples is DHS
	cutoff : float, optional
	    FPR cutoff value to label as occupied when building prior
	pseudo : float, optional
	    Psuedocount to add as prior to Beta distribution
	
	Returns
	-------
	TYPE
	    Description
	"""

	k = np.sum(fdr <= cutoff, axis = 0) # num of datasets w/ fp
	n = np.sum(w, axis = 0) #np.sum(w) # num of datasets w/ HS
	A = n-k + pseudo
	B = k + pseudo
	pr = A/(A+B)
	
	res = np.ones(fdr.shape)
	res *= pr[np.newaxis,:]
	res[w==0] = 1

	return res

def compute_delta_prior(obs, exp, fdr, beta_prior, cutoff = 0.05):
	"""Returns a point estimate of exepected cleavage depletion with a footprint per nucleotide 
	
	Parameters
	----------
	obs : TYPE
	    Description
	exp : TYPE
	    Description
	fdr : TYPE
	    Description
	beta_prior : TYPE
	    Description
	cutoff : float, optional
	    Description
	
	Returns
	-------
	TYPE
	    Description
	"""

	(n, w) = obs.shape

	mus = np.ones((n, w))
	ws = np.ones((n, w))
	
	for i in range(n):
		
		k = obs[i,:]
		n = np.max(np.vstack([exp[i,:], obs[i,:]]), axis = 0)

		mu, v = scipy.stats.beta.stats(k+beta_prior[i][0], n-k+beta_prior[i][1], 
			loc = 0, scale = 1, moments = 'mv')

		mus[i,:] = mu
		ws[i,:] = 1/np.sqrt(v)

	ws[fdr > cutoff] = 0

	delta = np.sum(ws*mus, axis = 0) / np.sum(ws, axis = 0)
	delta[np.isnan(delta)] = 1

	return delta

def log_likelihood(obs, exp, dm, delta = 1, w = 3):
	"""
	
	Parameters
	----------
	obs : TYPE
	    Description
	exp : TYPE
	    Description
	dm : TYPE
	    Description
	delta : int, optional
	    Description
	w : int, optional
	    Description
	
	Returns
	-------
	TYPE
	    Description
	"""
	res = np.ones((obs.shape[0], obs.shape[1]), order = 'c')

	n = obs.shape[0]
	for i in range(n):
		res[i,:] = windowing.sum(dm[i].log_pmf_values(exp[i,:] * delta, obs[i,:]), w)
	
	return res

def posterior(prior, ll_on, ll_off):
	"""
	Parameters
	----------
	prior : TYPE
	    Description
	ll_on : TYPE
	    Description
	ll_off : TYPE
	    Description
	
	Returns
	-------
	TYPE
	    Description
	"""
	prior_on = np.log(1-prior)
	prior_off = np.log(prior)

	p_off = prior_off + ll_off
	p_on = prior_on + ll_on
	denom = np.logaddexp(p_on, p_off)

	return (p_off - denom)




