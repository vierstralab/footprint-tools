# Copyright 2015 Jeff Vierstra

# cython: embedsignature=True



cimport numpy as np

ctypedef np.float64_t data_type_t

cpdef segment(data_type_t [:] x, data_type_t threshold, int w = 1, bint decreasing = 0):
	"""Segment an array into continuous elements passing a threshhold

	"""

	cdef double dir = -1 if decreasing else 1

	cdef list ret = []
	cdef int i, curr_start = -1
	for i in range(x.shape[0]):
		if curr_start < 0:
			if dir*x[i] >= dir*threshold:
				curr_start = i-w+1
		else:
			if dir*x[i] < dir*threshold:
				if len(ret) > 0 and curr_start <= ret[-1][1]:
					ret[-1][1] = i-1+w
				else:
					ret.append( [curr_start, i-1+w] )
				curr_start = -1
	return ret
