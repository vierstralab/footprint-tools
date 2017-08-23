# Copyright 2015 Jeff Vierstra

# cython: boundscheck=False
# cython: wraparound=False
# cython: nonecheck=False

from libc.stdlib cimport calloc, free

import numpy as np
cimport numpy as np

from scipy.spatial import cKDTree
from scipy.special import digamma

from bisect import bisect_left

ctypedef np.float64_t data_type_t

cdef extern from "mutual_information.h":
	int pairwise_mutual_information(double*, int, int, int, double*, int, int, int)
	void shuffle_columns_fast(double*, int, int)
	int pairwise_mutual_information_kraskov(double*, int, int, double*, int, int, int)

cpdef mi_kraskov(np.ndarray[data_type_t, ndim = 2, mode = 'c'] mat, int k = 3):
	
	cdef int m = mat.shape[0]
	cdef int n = mat.shape[1]

	cdef data_type_t [:,:] mat_view = mat

	cdef np.ndarray[data_type_t, ndim = 2, mode = 'c'] result = np.zeros((n, n), dtype = np.float64, order = 'c')
	cdef data_type_t [:,:] result_view = result

	cdef int success = pairwise_mutual_information_kraskov(&mat_view[0,0], m, n, &result_view[0,0], n, 1, 0)

	return result

cpdef mutual_information(np.ndarray[data_type_t, ndim = 2, mode = 'c'] mat, int bins = 3):

	cdef int nx = mat.shape[0]
	cdef int ny = mat.shape[1]

	cdef data_type_t [:,:] mat_view = mat

	cdef np.ndarray[data_type_t, ndim = 2, mode = 'c'] result = np.zeros((ny, ny), dtype = np.float64, order = 'c')
	cdef data_type_t [:,:] result_view = result

	cdef int success = pairwise_mutual_information(&mat_view[0,0], nx, ny, bins, &result_view[0,0], ny, 1, 0)

	return result


cpdef shuffle_columns(np.ndarray[data_type_t, ndim = 2, mode = 'c'] mat):
	
	cdef int nx = mat.shape[0]
	cdef int ny = mat.shape[1]

	cdef data_type_t [:,:] mat_view = mat

	shuffle_columns_fast(&mat_view[0, 0], nx, ny)

	return mat

cpdef mutual_information_p(np.ndarray[data_type_t, ndim = 2, mode = 'c'] mat, int bins = 3, int ntimes = 10000):
	
	cdef int i, j
	cdef int n = mat.shape[1] 
	
	cdef np.ndarray[data_type_t, ndim = 2, mode = 'c'] m = mutual_information(mat, bins = bins)

	cdef np.ndarray[data_type_t, ndim = 2, mode = 'c'] mat_rand = np.copy(mat, order = 'c')
	cdef data_type_t [:,:] mat_rand_view = mat_rand

	cdef np.ndarray[data_type_t, ndim = 3, mode = 'c'] m_null = np.zeros((n, n, ntimes), dtype = np.float64, order = 'c')
	cdef data_type_t [:,:,:] m_null_view = m_null

	cdef np.ndarray[data_type_t, ndim = 2, mode = 'c'] p = np.ones((n, n), dtype = np.float64, order = 'c')

	for i in range(ntimes):
		shuffle_columns_fast(&mat_rand_view[0, 0], mat_rand.shape[0], mat_rand.shape[1])
#		m_null[:,:,i] = mutual_information(mat_rand, bins = 2)
		pairwise_mutual_information(&mat_rand_view[0,0], mat_rand.shape[0], mat_rand.shape[1], bins, &m_null_view[0,0,0], ntimes*n, ntimes, i)
	
	m_null.sort(axis = 2)

	cdef int pos

	for i in range(n):
		for j in range(n):
			pos = bisect_left(m_null[i, j,:], m[i, j])
			p[i,j] = 1.0 - (<float>pos / <float>ntimes)

	p[p<=0] = 1.0/<float>ntimes
	return (m, m_null, p)





