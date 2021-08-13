cimport cython
cimport numpy as np

ctypedef np.float64_t data_type_t

cpdef data_type_t logpmf(data_type_t x, data_type_t nu, data_type_t mu, data_type_t sig2) nogil
