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
	#int pairwise_mutual_information_kraskov(double*, int, int, double*, int)

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

	cdef np.ndarray[data_type_t, ndim = 2, mode = 'c'] mat_rand = np.copy(mat, order = 'c')
	cdef data_type_t [:,:] mat_rand_view = mat_rand

	shuffle_columns_fast(&mat_rand_view[0, 0], mat_rand.shape[0], mat_rand.shape[1])

	return mat_rand

cpdef mutual_information_p(np.ndarray[data_type_t, ndim = 2, mode = 'c'] mat, int bins = 3, int ntimes = 10000):
	
	cdef int i, j
	cdef int n = mat.shape[1] 
	
	cdef np.ndarray[data_type_t, ndim = 2, mode = 'c'] mi = mutual_information(mat, bins = bins)

	cdef np.ndarray[data_type_t, ndim = 2, mode = 'c'] mat_rand = np.copy(mat, order = 'c')
	cdef data_type_t [:,:] mat_rand_view = mat_rand

	cdef np.ndarray[data_type_t, ndim = 3, mode = 'c'] mi_null = np.zeros((n, n, ntimes), dtype = np.float64, order = 'c')
	cdef data_type_t [:,:,:] mi_null_view = mi_null

	cdef np.ndarray[data_type_t, ndim = 2, mode = 'c'] p = np.ones((n, n), dtype = np.float64, order = 'c')

	for i in range(ntimes):
		shuffle_columns_fast(&mat_rand_view[0, 0], mat_rand.shape[0], mat_rand.shape[1])
#		m_null[:,:,i] = mutual_information(mat_rand, bins = 2)
		pairwise_mutual_information(&mat_rand_view[0,0], mat_rand.shape[0], mat_rand.shape[1], bins, &mi_null_view[0,0,0], ntimes*n, ntimes, i)
	
	mi_null.sort(axis = 2)

	cdef int pos

	for i in range(n):
		for j in range(n):
			pos = bisect_left(mi_null[i, j,:], mi[i, j])
			p[i,j] = 1.0 - (<float>pos / <float>ntimes)

	p[p<=0] = 1.0/<float>ntimes
	return (mi, mi_null, p)

from scipy.spatial import cKDTree
from scipy.special import digamma

cpdef compute_kraskov(data_type_t [:] x, data_type_t [:] y, treex, treey, int k = 3):
	
	cdef int n = x.shape[0]
	cdef double ptx, pty

	cdef double avgi = 0.0
	cdef double avgj = 0.0

	treexy = cKDTree(np.column_stack([x, y]))

	for i in range(n):

		ptx = x[i]
		pty = y[i]

		d = treexy.query([ptx, pty], k=k+1, p = np.inf)[0][k]

		avgi += digamma(len(treex.query_ball_point([ptx], d-1e-15, p = np.inf))) / n
		avgj += digamma(len(treey.query_ball_point([pty], d-1e-15, p = np.inf))) / n

	return (-avgi-avgj+digamma(3.0)+digamma(n)) / np.log(2)

cpdef mutual_information_kraskov(np.ndarray[data_type_t, ndim = 2, mode = 'c'] mat, int k = 3):

	cdef int i, j, l
	cdef int m = mat.shape[0]
	cdef int n = mat.shape[1]

	trees = [ cKDTree(mat[:,i][:,np.newaxis]) for i in range(n) ]

	cdef np.ndarray[data_type_t, ndim = 2, mode = 'c'] res = np.ones((n, n), dtype = np.float64)
	cdef np.ndarray[data_type_t, ndim = 3, mode = 'c'] res_rand = np.ones((n, n, 10000), dtype = np.float64) * -1.0e6

	cdef double best_rand_mi

	for i in range(n):
		for j in range(i+1, n):			
			
			res[i, j] = compute_kraskov(mat[:,i], mat[:,j], trees[i], trees[j], k = k)

			best_rand_mi = -1.0e6

			for s in range(1000):

				res_rand[i, j, s] = compute_kraskov(mat[:,i], np.random.permutation(mat[:,j]), trees[i], trees[j], k = k)
				best_rand_mi = res_rand[i, j, s] if res_rand[i, j, s] > best_rand_mi else best_rand_mi

				if s == 10 and best_rand_mi > res[i, j]:
					break

				if s == 100 and best_rand_mi > res[i, j]:
					break

	return res, res_rand

	

