cimport cython
cimport numpy as np

ctypedef np.float64_t data_type_t

cpdef data_type_t logpmf(int k, data_type_t p, data_type_t r) nogil
cpdef data_type_t pmf(int k, data_type_t p, data_type_t r) nogil
cpdef data_type_t cdf(int k, data_type_t p, data_type_t r) nogil
cpdef data_type_t mean(data_type_t p, data_type_t r) nogil
