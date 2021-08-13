"""
Implementation of Student's t distribution
"""
# cython: embedsignature=True

cimport cython
cimport numpy as np

cdef extern from "math.h":
    double log(double) nogil
    double sqrt(double) nogil

cdef extern from "hcephes.h":
    double hcephes_lgam(double) nogil
    double HCEPHES_PI

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
    return hcephes_lgam((nu+1.0)/2.0) - hcephes_lgam(nu/2.0) - log(sqrt(HCEPHES_PI*nu*sig2)) - ((nu+1.0)/2.0 * log(1.0+((x-mu)**2)/sig2/nu ))
