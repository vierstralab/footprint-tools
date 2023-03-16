"""
This sub-module contains functions for performing sliding window functions
written in native C. 
"""

# Copyright 2015-2021 Jeff Vierstra

# cython: boundscheck=False
# cython: wraparound=False
# cython: nonecheck=False
# cython: embedsignature=True

import numpy as np

from libc.stdlib cimport free
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
    """Windowing function wrapper that passes ndarray to
    fast functions written in native C

    Parameters
    ---------
    x : :class:`numpy.ndarray`
        Array of values to perform sliding window function
    hw : int
        Half window width
    func_ptr: callable
        Function to apply to each window 
    """
    cdef int i
    cdef int n = x.shape[0]

    cdef double* res = fast_windowing_func(&x[0], n, hw, func_ptr)
    cdef data_type_t [:] win = np.ones(n, dtype = np.float64, order = 'c')
    
    for i in range(hw, n-hw):
        win[i] = res[i]
    
    free(res)

    return np.asarray(win)

cpdef sum(data_type_t [:] x, int hw):
    """Sum of values in sliding local window across array

    Paramters
    ---------
    x : :class:`numpy.ndarray`
        Values to perform windowed sum. Must be memory
        contiguous (i.e., np.ascontiguousarray(x))
    hw: int
        Half window with to apply sum

    Returns
    -------
    out : ndarray
       Array of windowed summed values
    """
    return windowing_func(x, hw, fast_sum)

cpdef product(data_type_t [:] x, int hw):
    """Sum of values in sliding local window across array

    Paramters
    ---------
    x : :class:`numpy.ndarray`
        Values to perform windowed product. Must be memory
        contiguous (i.e., np.ascontiguousarray(x))
    hw: int
        Half window with to apply sum

    Returns
    -------
    out : :class:`numpy.ndarray`
       Array of windowed product values
    """
    return windowing_func(x, hw, fast_product)

cpdef fishers_combined(data_type_t [:] x, int hw):
    """Compute p-values for a window using Fisher's combined method

    Paramters
    ---------
    x : :class:`numpy.ndarray`
        P-values to perform Fisher's combined method. Array must 
        be memory contiguous (i.e., np.ascontiguousarray(x))
    hw: int
        Half window with to combined p-values

    Returns
    -------
    out : :class:`numpy.ndarray`
       Array of combined p-values
    """
    return windowing_func(x, hw, fast_fishers_combined)

cpdef stouffers_z(data_type_t [:] x, int hw):
    """Compute p-values for a window using Stouffer's Z-score method

    Paramters
    ---------
    x : :class:`numpy.ndarray`
        P-values to perform Stouffer's Z-score method. Array must 
        be memory contiguous (i.e., np.ascontiguousarray(x))
    hw: int
        Half window with to combined p-values

    Returns
    -------
    out : :class:`numpy.ndarray`
       Array of combined p-values
    """
    return windowing_func(x, hw, fast_stouffers_z)

cdef weighted_windowing_func(data_type_t [:] x, data_type_t [:] w, int hw, weighted_func_t func_ptr):
    """Weighted windowing function wrapper that passes ndarray to
    fast functions written in native C

    Parameters
    ---------
    x : :class:`numpy.ndarray`
        Array of values to perform sliding window function
    w: ndarray
        Weights for each array element
    hw : int
        Half window width
    func_ptr: callable
        Function to apply to each window 
    """
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

    Parameters
    ----------
    x : :class:`numpy.ndarray`
        P-values to perform Stouffer's Z-score method. Array must 
        be memory contiguous (i.e., np.ascontiguousarray(x))
    w: ndarray
        Weights for each elements in `x`
    hw: int
        Half window with to combined p-values

    Returns
    -------
    out : :class:`numpy.ndarray`
       Array of combined p-values
    """
    return weighted_windowing_func(x, w, hw, fast_weighted_stouffers_z)
