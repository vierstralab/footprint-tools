"""
Implementation of negative binomial distribution
"""
# cython: embedsignature=True

import scipy.special
import scipy.optimize
import numpy as np

cimport cython
cimport numpy as np

cdef extern from "math.h":
    double log(double) nogil
    double sqrt(double) nogil
    double exp(double) nogil

cdef extern from "hcephes.h":
    double hcephes_log1p(double) nogil
    double hcephes_lgam(double) nogil
    double hcephes_incbet(double, double, double) nogil

import warnings

def mle(par, data, sm):
    """Objective function for MLE estimate according to
    https://en.wikipedia.org/wiki/Negative_binomial_distribution#Maximum_likelihood_estimation
    
    Parameters
    ----------
    par : list
        Distribution parameters that are iteratively updated
    data : iterable list or :class:`np.array`
        Data points to be fit
    sm: \sum data / len(data)
    
    Returns
    -------
    out : ndarray
        Likelihood of parameters given data
    """

    p = par[0]
    r = par[1]
    n = len(data)
    f0 = sm / (r+sm) - p
    f1 = np.sum( scipy.special.psi(data+r) ) - n * scipy.special.psi(r) + n * np.log( r / (r+sm) )
    
    return np.array([f0, f1])

def fit(data, p = None, r = None):
    """Function to fit the parameters of a negative binomial to data
    
    Parameters
    ----------
    data : :class:`np.array`
        Data to fit (1-D array)
    p : float, optional
        Initial guess for `p`
    r : float, optional
        Initial guess for `r`
    
    Returns
    -------
    out: tuple
       Maximum likelihood fit parameters (p, r)
    """
    if p is None or r is None:
        av = np.average(data)
        va = np.var(data)
        r = (av*av) / (va-av)
        p = (va-av) / (va)

    sm = np.sum(data)/len(data)
    
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        x = scipy.optimize.fsolve(mle, np.array([p, r]), args = (data, sm))

    return (x[0], x[1])

cpdef data_type_t logpmf(int k, data_type_t p, data_type_t r) nogil:
    """log of probability mass function for negative binomial

    Parameters
    ----------
    k : int
        Number of successes
    p : float
        Probability of success
    r: float
        Number of failures

    Returns
    -------
    p : float
        log of probability mass at value `k`
    """
    cdef double coeff = hcephes_lgam(k+r) - hcephes_lgam(k+1) - hcephes_lgam(r)
    return coeff + r * log(p) + k * hcephes_log1p(-p)

cpdef data_type_t pmf(int k, data_type_t p, data_type_t r) nogil:
    """Probability mass function for negative binomial

    Parameters
    ----------
    k : int
        Number of successes
    p : float
        Probability of success
    r: float
        Number of failures

    Returns
    -------
    p : float
        Probability mass at value `k`
    """
    return exp(logpmf(k, p, r))

cpdef data_type_t cdf(int k, data_type_t p, data_type_t r) nogil:
    """Cumulative distribution function for negative binomial

    Parameters
    ----------
    k : int
        Number of successes
    p : float
        Probability of success
    r: float
        Number of failures

    Returns
    -------
    p : float
        Cumulative density at value `k`
    """
    return hcephes_incbet(r, k+1, p)

cpdef data_type_t mean(data_type_t p, data_type_t r) nogil:
    """Mean of negative binomial distribution

    Parameters
    ----------
    p : float
        Probability of success
    r: float
        Number of failures

    Returns
    -------
    mu : float
        Mean of distribution
    """
    return p*r/(1-p)

cpdef data_type_t var(data_type_t p, data_type_t r) nogil:
    """Variance of negative binomial distribution

    Parameters
    ----------
    p : float
        Probability of success
    r: float
        Number of failures

    Returns
    -------
    var : float
        Variance of distribution
    """
    return (p*r)/((1-p)*(1-p))

cpdef data_type_t rvs(data_type_t p, data_type_t r) nogil:
    """Sample a random variate from negative binomial distribution
    
    Parameters
    ----------
    p : float
        Probability of success
    r: float
        Number of failures

    Returns
    -------
    x : float
        Random variate
    """
    raise NotImplementedError