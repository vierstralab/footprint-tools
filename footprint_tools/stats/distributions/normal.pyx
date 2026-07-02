"""
Implementation of normal distribution
"""
# cython: embedsignature=True

import numpy as np

cimport cython
cimport numpy as np

cdef extern from "math.h":
    double pow(double) nogil
    double log(double) nogil
    double sqrt(double) nogil

cdef extern from "hcephes.h":
    double HCEPHES_PI

cpdef data_type_t logpdf(
    data_type_t x,
    data_type_t mu,
    data_type_t sig2,
) noexcept nogil:
    cdef data_type_t d = x - mu
    return -0.5 * (d * d / sig2 + log(2.0 * HCEPHES_PI * sig2))
