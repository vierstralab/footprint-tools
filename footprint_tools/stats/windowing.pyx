# Copyright 2015 Jeff Vierstra

# cython: boundscheck=False
# cython: wraparound=False
# cython: nonecheck=False
# cython: embedsignature=True


from libc.stdlib cimport free


import numpy as np

cimport numpy as np

ctypedef np.float64_t data_type_t

ctypedef double (*func_t)(const double *, int)
ctypedef double (*weighted_func_t)(const double* const, const double* const, int)

cdef extern from "windowing.h":
    #
    double fast_sum(const double* const, int)
    double fast_product(const double* const, int)
    double fast_fishers_combined(const double* const, int)
    double fast_stouffers_z(const double* const, int)
    double* fast_windowing_func(const double* const, int, int, func_t)
    #
    double fast_weighted_stouffers_z(const double* const, const double* const, int)
    double* fast_weighted_windowing_func(const double* const, const double* const, int, int, weighted_func_t)

cdef windowing_func(data_type_t [:] x, int hw, func_t func_ptr):
    
    cdef int i
    cdef int n = x.shape[0]

    cdef double* res = fast_windowing_func(&x[0], n, hw, func_ptr)
    cdef data_type_t [:] win = np.ones(n, dtype = np.float64, order = 'c')
    
    for i in range(hw, n-hw):
        win[i] = res[i]
    
    free(res)

    return np.asarray(win)

cpdef sum(data_type_t [:] x, int hw):
    """Sum of values in sliding local window

    :param x: (array) values
    :param hw: (int) half window size to smooth
    
    :return: (float) summed values
    """

    return windowing_func(x, hw, fast_sum)

cpdef product(data_type_t [:] x, int hw):
    """Product of values in sliding local window

    :param x: (array) values
    :param hw: (int) half window size to smooth
    
    :return: (float) product of values
    """
    
    return windowing_func(x, hw, fast_product)

cpdef fishers_combined(data_type_t [:] x, int hw):
    """Compute p-value for a window using Fisher's combined method

    :param x: (array) p-values
    :param hw: (int) half window size to smooth
    
    :return p: (float) combined p-values
    """


    return windowing_func(x, hw, fast_fishers_combined)

cpdef stouffers_z(data_type_t [:] x, int hw):
    """Compute p-value for a window using Stouffer's method

    :param x: (array) p-values
    :param hw: (int) half window size to smooth
    
    :return p: (float) combined p-values
    """

    return windowing_func(x, hw, fast_stouffers_z)

cdef weighted_windowing_func(data_type_t [:] x, data_type_t [:] w, int hw, weighted_func_t func_ptr):
    
    cdef int i
    cdef int n = x.shape[0]

    cdef double* res = fast_weighted_windowing_func(&x[0], &w[0], n, hw, func_ptr)
    cdef data_type_t [:] win = np.ones(n, dtype = np.float64, order = 'c')
    
    for i in range(hw, n-hw):
        win[i] = res[i]
    
    free(res)

    return np.asarray(win)

cpdef weighted_stouffers_z(data_type_t [:] x, data_type_t [:] w, int hw):
    """Compute p-value for a window using weighted Stouffer's method

    :param x: (array) p-values
    :param w: (array) weights
    :param hw: (int) half window size to smooth
    
    :return p: (float) combined p-values
    """

    return weighted_windowing_func(x, w, hw, fast_weighted_stouffers_z)
