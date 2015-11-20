# cython: boundscheck=False
# cython: wraparound=False
# cython: nonecheck=False

import numpy as np
import scipy.stats

cimport cython
cimport numpy as np

ctypedef np.float64_t data_type_t

cdef extern from "cephes.h":
    double c_chdtr(double, double)
    double c_ndtri(double)
    double c_ndtr(double)
    double c_sqrt(double)

cpdef fishers_combined(data_type_t [:] x):
    
    cdef data_type_t chi, p
  
    chi = -2.0 * np.sum(np.log(x))
    p = 1.0-c_chdtr(2*x.shape[0], chi)

    return p

cpdef stouffers_z(data_type_t [:] x):
    
    cdef int i, n = x.shape[0]
    cdef data_type_t z, p, s = 0.0

    for i in range(n):
        s += c_ndtri(1.0-x[i])
    
    z = s / c_sqrt(n)
    p = 1.0-c_ndtr(z)

    return p

cpdef windowed_p_value(data_type_t [:] x, int w, func_ptr):

    """Compute z-scores from groups of p-values
    Args:
        x (array): p-values
        w (int): window size to smooth
        func: test statistic function
    Returns:
        p (float): p-values
    """
    
    cdef int i, n = x.shape[0]
    cdef np.ndarray[data_type_t, ndim = 1] p = np.ones(n, dtype = np.float64)
    cdef data_type_t [:] p_view = p
    
    for i in range(w, n-w+1):
        p_view[i] = func_ptr(x[i-w:i+w+1])
    
    return p

