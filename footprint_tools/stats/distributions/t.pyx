"""
Implementation of Student's t distribution
"""
# cython: embedsignature=True

cimport cython
cimport numpy as np

cdef extern from "cephes.h":
    double c_lgamma(double) nogil
    double c_log(double) nogil
    double c_sqrt(double) nogil
    double PI

ctypedef np.float64_t data_type_t

cpdef data_type_t logpmf(data_type_t x, data_type_t nu, data_type_t mu, data_type_t sig2) nogil:
    """Log probability mass function of Student's t distribution
    
    Parameters
    ----------
    x : float

    nu : float

    mu : float

    sig2 : float

    Returns
    -------
    p : float
    """
    return c_lgamma((nu+1.0)/2.0) - c_lgamma(nu/2.0) - c_log(c_sqrt(PI*nu*sig2)) - ((nu+1.0)/2.0 * c_log(1.0+((x-mu)**2)/sig2/nu ))
