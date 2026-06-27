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

cpdef data_type_t logpmf(data_type_t x, data_type_t mu, data_type_t sig2) nogil:
	cdef data_type_t sig = sqrt(sig2)
	return -0.5 * pow((x-mu)/sig), 2) - log(sig*sqrt(2*HCEPHES_PI))
