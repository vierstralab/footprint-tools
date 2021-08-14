# Copyright 2015-2021 Jeff Vierstra

cimport numpy as np

ctypedef np.float64_t data_type_t

cpdef data_type_t piecewise_four(x, data_type_t x0, data_type_t x1, data_type_t x2, data_type_t x3, 
                                     data_type_t y0, data_type_t y1, data_type_t y2, data_type_t y3,
                                     data_type_t k0, data_type_t k1, data_type_t k2, data_type_t k3)


cpdef data_type_t piecewise_three(x, data_type_t x0, data_type_t x1, data_type_t x2,
                                     data_type_t y0, data_type_t y1, data_type_t y2,
                                     data_type_t k0, data_type_t k1, data_type_t k2)

cdef class dispersion_model:
    
    cdef np.ndarray _h
    cdef np.ndarray _p
    cdef np.ndarray _r
    cdef np.ndarray _r_params
    cdef np.ndarray _mu_params
    
    cpdef str _metadata

    cpdef data_type_t fit_mu(self, data_type_t x)
    cpdef data_type_t fit_r(self, data_type_t x)

    cpdef data_type_t [:] log_pmf_values(self, data_type_t [:] exp, data_type_t [:] obs)
    cpdef data_type_t [:] pmf_values(self, data_type_t [:] exp, data_type_t [:] obs)

    cpdef data_type_t [:] log_pmf_values_0(self, data_type_t [:] exp, data_type_t [:] obs, data_type_t [:] res)
    cpdef data_type_t [:] pmf_values_0(self, data_type_t [:] exp, data_type_t [:] obs, data_type_t [:] res)

    cpdef data_type_t [:] p_values(self, data_type_t [:] exp, data_type_t [:] obs)
    cpdef sample(self, data_type_t [:] x, int times)
    