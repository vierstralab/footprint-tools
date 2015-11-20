
cimport numpy as np

ctypedef np.float64_t data_type_t

cpdef segment(data_type_t [:] x, data_type_t threshold):
	"""Segment an array into continuous elements passing a threshhold

	Returns
	-------
	segments (list of tuples): start and end points to regions that pass a threshold
	"""
	cdef list ret = []
	cdef int i, curr_start = -1
	for i in range(x.shape[0]):
		if curr_start < 0:
			if x[i] >= threshold:
				curr_start = i-3
		else:
			if x[i] < threshold:
				ret.append( (curr_start, i-1+3) )
				curr_start = -1
	return ret
