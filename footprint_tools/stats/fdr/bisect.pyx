# Copyright 2015 Jeff Vierstra

# cython: boundscheck=False
# cython: wraparound=False
# cython: nonecheck=False
# cython: embedsignature=True


import numpy as np

cimport cython
cimport numpy as np

ctypedef np.float64_t data_type_t

cpdef bisect(data_type_t [:] a, data_type_t [:] b):
	"""
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