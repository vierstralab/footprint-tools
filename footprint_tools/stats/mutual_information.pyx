# Copyright 2015 Jeff Vierstra

# cython: boundscheck=False
# cython: wraparound=False
# cython: nonecheck=False

from libc.stdlib cimport calloc, free

import numpy as np
import ctypes

cimport numpy as np

ctypedef np.float64_t data_type_t

cdef extern from "mutual_information.h":
	int pairwise_mutual_information(int*, int, int, double *)


cpdef mutual_information(np.ndarray[int, ndim = 2, mode = 'c'] mat):

	cdef int nx = mat.shape[0]
	cdef int ny = mat.shape[1]

	cdef int [:,:] mat_view = mat

	cdef np.ndarray[data_type_t, ndim = 2, mode = 'c'] result = np.zeros((ny, ny), dtype = np.float64, order = 'c')
	cdef data_type_t [:,:] result_view = result

	cdef int success = pairwise_mutual_information(&mat_view[0,0], nx, ny, &result_view[0,0])

	return result
