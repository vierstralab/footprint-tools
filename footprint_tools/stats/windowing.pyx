# cython: boundscheck=False
# cython: wraparound=False
# cython: nonecheck=False

import numpy as np
import scipy.stats

cimport cython
cimport numpy as np

ctypedef np.float64_t data_type_t

cdef extern from "cephes.h":
    double chdtr(double, double)
    double ndtri(double)
    double ndtr(double)

cpdef fishers_combined(data_type_t [:] x):
    cdef data_type_t chi, p
  
    chi = -2.0 * np.sum(np.log(x))
    p = 1.0-chdtr(2*x.shape[0], chi)
    
    return chi, p

cpdef stouffers_z(data_type_t [:] x):
    
    cdef int i, n = x.shape[0]
    cdef data_type_t z, p, s = 0.0

    for i in range(0, n):
        s += ndtri(1.0-x[i])
    
    z = s / sqrt(n)
    p = 1.0-ndtr(z)

    return z, p

cpdef windowed_p_value(data_type_t [:] x, int w, func_ptr):

    """Compute z-scores from groups of p-values
    Args:
        x (array): p-values
        w (int): window size to smooth
        func: test statistic function
    Returns:
        (chi, p): tuple of test statistics and p-values
    """
    
    cdef int i, n = x.shape[0]
    cdef np.ndarray[data_type_t, ndim = 1] z = np.zeros(n, dtype = np.float64)
    cdef np.ndarray[data_type_t, ndim = 1] p = np.ones(n, dtype = np.float64)

    for i in range(w, n-w+1):
        z[i], p[i] = func_ptr(x[i-w:i+w+1])
    return z, p
