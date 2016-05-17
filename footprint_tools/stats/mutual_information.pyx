# Copyright 2016 Jeff Vierstra

# cython: boundscheck=False
# cython: wraparound=False
# cython: nonecheck=False

from libc.stdlib cimport malloc, free

import numpy as np
##
cimport numpy as np

ctypedef np.float64_t data_type_t

cdef extern from "mutual_information.h":
    #
	double entropy(double)
	double* joint_probability(int*, int*, int)

cpdef mutual_info(data_type_t [:,:] x, thresh):
	
	cdef int i
	cdef int j
	
	cdef int m = x.shape[0]
	cdef int l = x.shape[1]

	cdef int[:,:] t = np.zeros((l, m), dtype = np.intc, order = 'c')
	for i in range(m):
		for j in range(l):
			t[j, i] = 1 if x[i, j] <= thresh else 0

	cdef double[:] p = np.sum(t, axis = 1)/<double>m

	cdef double[:,:] pxy = np.zeros((l, l), dtype = np.float64, order = 'c')
	cdef double *z

	for i in range(l):
		for j in range(l):
			print(np.array(t[i,:]))
			print(np.array(t[j,:]))
			z = joint_probability(&t[i,0], &t[j,0], m)
			print(z[0], z[1],z[2],z[3])
	return pxy