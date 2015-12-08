# Copyright 2015 Jeff Vierstra

# cython: boundscheck=False
# cython: wraparound=False
# cython: nonecheck=False

import numpy as np

cimport numpy as np

ctypedef np.float64_t data_type_t

ctypedef double (*func_t)(const double *, int);

cdef extern from "windowing.h":
    double fast_fishers_combined(const double* const, int)
    double fast_stouffers_z(const double* const, int)
    double* fast_windowing_func(const double* const, int, int, func_t)

cdef windowing_func(data_type_t [:] x, int w, func_t func_ptr):
    
    cdef int i
    cdef int n = x.shape[0]

    cdef double* res = fast_windowing_func(&x[0], n, w, func_ptr)
    cdef data_type_t [:] win = np.ones(n, dtype = np.float64, order = 'c')
    
    for i in range(w, n-w):
        win[i] = res[i]
    
    free(res)

    return win

cpdef stouffers_z(data_type_t [:] x, int w):
    
    """Compute p-value for a window using
        Stouffer's method
    Args:
        x (array): p-values
        w (int): half window size to smooth
    Returns:
        p (float): p-values
    """ 

    return windowing_func(x, w, fast_stouffers_z)

cpdef fishers_combined(data_type_t [:] x, int w):
    
    """Compute p-value for a window using
        Fisher combined method
    Args:
        x (array): p-values
        w (int): half window size to smooth
    Returns:
        p (float): p-values
    """

    return windowing_func(x, w, fast_fishers_combined)