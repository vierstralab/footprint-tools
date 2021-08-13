"""
This sub-module contains functions compute the posterior probability that
a nucleotide is footprinted.
"""

import numpy as np
import scipy.stats

from footprint_tools.stats import windowing

def compute_prior_weighted(fdr, w, cutoff = 0.05, pseudocount = 0.5):
    """Returns prior of whether a nucleotide is footprinted
    
    Parameters
    ----------
    fdr : ndarray
        FDR values for each datset and position
    w : ndarray
        Binary array whether nucleotide is within a DHS
    cutoff : float, optional
        FPR cutoff value to label as occupied when building prior
    pseudocount : float, optional
        Psuedocount to add as prior to Beta distribution
    
    Returns
    -------
    :class:`np.array`
        Per-nucleotide occupancy prior
    """

    k = np.sum(fdr <= cutoff, axis = 0) # num of datasets w/ fp
    n = np.sum(w, axis = 0) #np.sum(w) # num of datasets w/ HS
    a = n-k + pseudocount
    b = k + pseudocount
    pr = a/(a+b)
    
    res = np.ones(fdr.shape)
    res *= pr[np.newaxis,:]
    res[w==0] = 1

    return res

def compute_delta_prior(obs, exp, fdr, beta_prior, cutoff = 0.05):
    """Returns a point estimate of exepected cleavage 
    depletion with a footprint at each nucleotide 
    
    Parameters
    ----------
    obs : ndarray
        Array with observed cleavages counts (2D - nucleotides x samples)
    exp : ndarray
        Array with expectd cleavages counts (2D - nucleotides x samples)
    fdr : ndarray
        Array with adjusted p-values (2D - nucleotides x samples)
    beta_prior : ndarray
        Description
    cutoff : float, optional
        FDR cutoff value to label a nucleotide as occupied when building prior
    
    Returns
    -------
    delta: ndarray
        Per-nucleotide priors of expected nucleotide protection at footprinted sites 
        (values explicitly in [0, 1]). Arrays is 1-D corresponding to DHS
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
    """Likelihood function of observed counts given the bias corrected data
    and an expected protection
    
    Parameters
    ----------
    obs : ndarray
        Observed cleavage counts
    exp : ndarray
        Expected cleavage counts
    dm : :class:`dispersion_model`
        Dispersional model to compute cleavage statistics
    delta : array_like, optional
        Values to scale exp counts
    w : int, optional
        Half window width to compute log-likelihood. Default 3bp (combined 7bp window)
    
    Returns
    -------
    log_likelihood: ndarray
        Array of log-likelihoods for each nucloetides and samples
    """
    res = np.ones((obs.shape[0], obs.shape[1]), order = 'c')

    n = obs.shape[0]
    for i in range(n):
        res[i,:] = windowing.sum(dm[i].log_pmf_values(exp[i,:] * delta, obs[i,:]), w)
    
    return res

def posterior(prior, ll_on, ll_off):
    """Compute the posterior probability of a nucleotide

    Parameters
    ----------
    prior : ndarray
        Per-nucleotide prior that a site is footprinted
    ll_on : ndarray
        Log-likelihood that a nucleotide is occupied
    ll_off : ndarray
        Log-likelihood that a nucleotide is unoccupied
    
    Returns
    -------
    posteriors : ndarray
        log posterior that a nucleotde is footprinted in sample. Array is
        2-D corresponding to nucleotides and samples. 
    """
    prior_on = np.log(1-prior)
    prior_off = np.log(prior)

    p_off = prior_off + ll_off
    p_on = prior_on + ll_on
    denom = np.logaddexp(p_on, p_off)

    return (p_off - denom)
