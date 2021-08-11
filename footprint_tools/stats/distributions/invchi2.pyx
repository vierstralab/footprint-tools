"""
Implementation of inverse Chi-squared distribution
"""

cimport cython
cimport numpy as np

cdef extern from "cephes.h":
    double c_log(double) nogil
    double c_lgamma(double) nogil
     
cpdef data_type_t logpmf(data_type_t x, data_type_t nu, data_type_t tau2) nogil:
    """Log probability mass function for inverse Chi-squared distribution

    Parameters
    ----------
    x : float

    nu : float

    tau2 : float

    Returns
    -------
    p : float

    """
    #return (nu/2.0)*np.log(tau2*nu/2.0) + (-nu*tau2/2.0/x) - scipy.special.gammaln(nu/2.0) - (1.0+nu/2.0)*np.log(x)
    return nu/2.0 * c_log(tau2*nu/2.0) + (-nu*tau2/2.0/x) - c_lgamma(nu/2.0) - (1.0+nu/2.0)*c_log(x)

cpdef data_type_t log_likelihood(data_type_t [:] x, data_type_t nu, data_type_t tau2) nogil:
    """Log-likelihood for inverse Chi-squared distribution

    Parameters
    ----------
    data : ndarray
        Data values
    nu : float

    tau2 : float

    Returns
    -------
    log_likelihood : float
    """
    cdef int n = x.shape[0]
    cdef data_type_t res = 0

    for i in range(n):
        ll += logpmf(x[i], nu, tau2)
    return ll
