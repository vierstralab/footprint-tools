"""
This sub-module contains utility functions for stats module 
"""

# cython: boundscheck=False
# cython: wraparound=False
# cython: nonecheck=False
# cython: embedsignature=True

cimport numpy as np
import numpy as np

ctypedef np.float64_t data_type_t

cpdef segment(data_type_t [:] x, data_type_t threshold, int w = 1, bint decreasing = 0):
    """Segment an array into continuous elements passing a threshold

    Parameters
    ----------
    x : array_like
        Array of values to segment
    threshold : float
        Threshold for grouping elements
    w : int
        Window size
    decreasing : bool

    Returns
    -------
    out : list
        List of tuples which specific the start and end index 
        of contiguous intervals that pass threshold

    """
    cdef double dir = -1 if decreasing else 1

    cdef list ret = []
    cdef int i, curr_start = -1
    for i in range(x.shape[0]):
        if curr_start < 0:
            if dir*x[i] >= dir*threshold:
                curr_start = i-w+1
        else:
            if dir*x[i] < dir*threshold:
                if len(ret) > 0 and curr_start <= ret[-1][1]:
                    ret[-1][1] = i-1+w
                else:
                    ret.append( [curr_start, i-1+w] )
                curr_start = -1
    return ret

cpdef bisect(data_type_t [:] a, data_type_t [:] b):
    """Bisect arrays

    Parameters
    ----------
    a : array_like
        
    b : array_like

    Returns
    -------
    out : ndarray
        Bisected list
    """
    cdef int i
    cdef int n = b.shape[0]

    cdef int lo = 0
    cdef int hi = a.shape[0]

    cdef np.ndarray[data_type_t, ndim = 1] ind = np.zeros(n, dtype = np.float64)

    for i in range(n):
        while lo < hi:
            if b[i] < a[lo]: break
            else: lo = lo + 1
        ind[i] = lo

    return ind