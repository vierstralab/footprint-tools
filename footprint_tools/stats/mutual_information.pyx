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

	cdef np.ndarray[data_type_t, ndim = 2, mode = 'c'] result = np.zeros((ny, ny), dtype = np.float64, order = 'c')

	cdef double* c_result = <double*>calloc(ny * ny, sizeof(double))

	cdef int success = pairwise_mutual_information(&mat[0,0], nx, ny, &c_result[0])

	cdef int i, j
	for i in range(ny):
		for j in range(ny):
			result[i, j] = c_result[i * ny + j]

	return result
