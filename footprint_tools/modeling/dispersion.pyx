"""
This module contains classess and functions that implement a dispersion model.
"""

# cython: boundscheck=False
# cython: wraparound=False
# cython: nonecheck=False
# cython: embedsignature=True

cimport cython
cimport numpy as np
import numpy as np

import footprint_tools

from footprint_tools.stats.distributions cimport nbinom
from footprint_tools.stats.distributions import nbinom

import warnings
from scipy import optimize
import pwlf

import logging
logger = logging.getLogger(__name__)

cpdef data_type_t piecewise_three(x, data_type_t x0, data_type_t x1, data_type_t x2,
                                     data_type_t y0, data_type_t y1, data_type_t y2,
                                     data_type_t k0, data_type_t k1, data_type_t k2):
    
    return (
        (x<x0) * (y0 + k0*(x)) +
        ((x>=x0) & (x<x1)) * (y1 + k1*(x)) +
        (x>=x1) * (y2 + k2*(x))
    )

cpdef data_type_t piecewise_four(x, data_type_t x0, data_type_t x1, data_type_t x2, data_type_t x3, 
                                     data_type_t y0, data_type_t y1, data_type_t y2, data_type_t y3,
                                     data_type_t k0, data_type_t k1, data_type_t k2, data_type_t k3):
    
    return (
        (x<x0) * (y0 + k0*(x)) +
        ((x>=x0) & (x<x1)) * (y1 + k1*(x)) +
        ((x>=x1) & (x<x2)) * (y2 + k2*(x)) +
        (x>=x2) * (y3 + k3*(x))
    )

cpdef data_type_t piecewise_five(x, data_type_t x0, data_type_t x1, data_type_t x2, data_type_t x3, data_type_t x4, 
                                     data_type_t y0, data_type_t y1, data_type_t y2, data_type_t y3, data_type_t y4,
                                     data_type_t k0, data_type_t k1, data_type_t k2, data_type_t k3, data_type_t k4):
    
    return (
        (x<x0) * (y0 + k0*(x)) +
        ((x>=x0) & (x<x1)) * (y1 + k1*(x)) +
        ((x>=x1) & (x<x2)) * (y2 + k2*(x)) +
        ((x>=x2) & (x<x3)) * (y3 + k3*(x)) +
        (x>=x3) * (y4 + k4*(x))
    )

cdef class dispersion_model:
    """
    Dispersion model class
    """
    def __init__(self):
        # Histograms
        self._h = None
        
        # Emperical fits
        self._p = None
        self._r = None

        self._mu_params = self._r_params = None

        self._metadata = ''

    # Pickling function
    def __reduce__(self):
        x = {}
        x['mu_params'] = self.mu_params
        x['r_params'] = self.r_params
        return (dispersion_model, (), x)

    # Pickling function
    def __setstate__(self, x):
        self.mu_params = x['mu_params']
        self.r_params = x['r_params']

    property h:
        """Histrogram of observed cleavages at each predicted cleavage rate"""
        def __get__(self):
            return self._h
        def __set__(self, x):
            self._h = x

    property p:
        """Array of the negative binomial MLE fit parameters `p`"""
        def __get__(self):
            return self._p
        def __set__(self, x):
            self._p = x

    property r:
        """Array of the negative binomial MLE fit parameters `r`"""
        def __get__(self):
            return self._r
        def __set__(self, x):
            self._r = x
    
    property metadata:
        def __get__(self):
            return self._metadata
        def __set__(self, x):
            self._metadata = x
    
    property mu_params:
        def __get__(self):
            return self._mu_params
        def __set__(self, x):
            self._mu_params = np.array(x, order = 'c')

    property r_params:
        def __get__(self):
            return self._r_params
        def __set__(self, x):
            self._r_params = np.array(x, order = 'c')


    cpdef data_type_t fit_mu(self, data_type_t x):
        """Computes the fitted mu term for the negative binomial
        from a piece-wise linear fit.
        
        Parameters
        ----------
        x : float

        Returns
        -------
        mu : float
            mu computed from the regression fit
        """

        cdef data_type_t [:] par = self._mu_params
        cdef data_type_t res = piecewise_three(x, *par)

        return res if res > 0.0 else 0.1

    cpdef data_type_t fit_r(self, data_type_t x):
        """Computes the dispersion term for the negative binomial
        from a piece-wise linear fit. Note that the model parameters 
        estimate the inverse.
        
        Parameters
        ----------
        x : float
        
        Returns
        -------
        r : float
            r computed from the regression fit
        """
        cdef data_type_t [:] par = self._r_params
        cdef data_type_t res = 1.0/piecewise_five(x, *par)

        return res if res > 0.0 else 1e-6

    def __str__(self):
        """Print model to string"""
        raise NotImplementedError

    @cython.cdivision(True)
    cpdef data_type_t [:] log_pmf_values(self, data_type_t [:] exp, data_type_t [:] obs):
        """Compute the log probability mass function

        Parameters
        ----------
        exp: :class:`numpy.ndarray`
            Expected cleavage counts
        obs: :class:`numpy.ndarray`
            Observed cleavage counts

        Returns
        -------
        logp : :class:`numpy.ndarray`
            Array of log probability mass function values computed from 
            the expected cleavage distributions
        """
        cdef int i, n = exp.shape[0]
        cdef double r, mu
        
        cdef data_type_t [:] res = np.zeros(n, dtype = np.float64, order = 'c')
    
        for i in range(n):
            r = self.fit_r(exp[i])
            mu = self.fit_mu(exp[i])
            res[i] = nbinom.logpmf(<int>obs[i], (r/(r+mu)), r)
        
        return res

    @cython.cdivision(True)
    cpdef data_type_t [:] pmf_values(self, data_type_t [:] exp, data_type_t [:] obs):
        """Compute the probability mass function

        Parameters
        ----------
        exp: :class:`numpy.ndarray`
            Expected cleavage counts
        obs: :class:`numpy.ndarray`
            Observed cleavage counts

        Returns
        -------
        p : :class:`numpy.ndarray`
            Array of probability mass function values computed from 
            the expected cleavage distributions
        """
        cdef int i, n = exp.shape[0]
        cdef double r, mu
        
        cdef data_type_t [:] res = np.zeros(n, dtype = np.float64, order = 'c')

        for i in range(n):
            r = self.fit_r(exp[i])
            mu = self.fit_mu(exp[i])
            res[i] = nbinom.pmf(<int>obs[i], (r/(r+mu)), r)

        return res

    @cython.cdivision(True)
    cpdef data_type_t [:] log_pmf_values_0(self, data_type_t [:] exp, data_type_t [:] obs, data_type_t [:] res):
        """Computing the log probability mass function to pointer

        Parameters
        ----------
        exp: :class:`numpy.ndarray`
            Expected cleavage counts
        obs: :class:`numpy.ndarray`
            Observed cleavage counts

        Returns
        -------
        logp : :class:`numpy.ndarray` (memoryview)
            Array pointer to log probability mass function values computed from 
            the expected cleavage distributions

        Notes
        -----
        This function is equivalent to `log_pmf_values`, except it stores values to
        a matrix pointer
        """
        cdef int i, n = exp.shape[0]
        cdef double r, mu
        
        for i in range(n):
            r = self.fit_r(exp[i])
            mu = self.fit_mu(exp[i])
            res[i] = nbinom.logpmf(<int>obs[i], (r/(r+mu)), r)

        return res
    
    @cython.cdivision(True)
    cpdef data_type_t [:] pmf_values_0(self, data_type_t [:] exp, data_type_t [:] obs, data_type_t [:] res):
        """Compute the probability mass function to pointer

        Parameters
        ----------
        exp: :class:`numpy.ndarray`
            Expected cleavage counts
        obs: :class:`numpy.ndarray`
            Observed cleavage counts

        Returns
        -------
        p : :class:`numpy.ndarray` (memoryview)
            Array pointer to probability mass function values computed from 
            the expected cleavage distributions

        Notes
        -----
        This function is equivalent to `pmf_values`, except it stores values to
        a matrix pointer
        """
        cdef int i, n = exp.shape[0]
        cdef double r, mu
        
        for i in range(n):
            r = self.fit_r(exp[i])
            mu = self.fit_mu(exp[i])
            res[i] = nbinom.pmf(<int>obs[i], (r/(r+mu)), r)

        return res

    cpdef data_type_t [:] p_values(self, data_type_t [:] exp, data_type_t [:] obs):
        """Compute cumulative distribution (lower-tail p-value)
        from negative binomial

        Parameters
        ----------
        exp: :class:`numpy.ndarray`
            Expected cleavage counts
        obs: :class:`numpy.ndarray`
            Observed cleavage counts

        Returns
        -------
        pvals : :class:`numpy.ndarray`
            Array of p-values
        """     
        cdef int i, n = exp.shape[0]
        cdef double r, mu
        cdef data_type_t [:] res = np.ones(n, dtype = np.float64, order = 'c')

        for i in range(n):
            r = self.fit_r(exp[i])
            mu = self.fit_mu(exp[i])
            res[i] = nbinom.cdf(<int>obs[i], r/(r+mu), r)

        return res

    cpdef sample(self, data_type_t [:] x, int times):
        """Sample counts from negative binomial distribution and
        compute p-values

        Parameters
        ----------
        x : :class:`numpy.ndarray`
            Count values to specifying from which distribution
            to resample. This typically expected count values.
        times : int
            Number of times to sample (per element)

        Returns
        -------
        sampled_counts : :class:`numpy.ndarray`
            Array of sample counts (2-D array - positions by number of samples)
        sampled_pvals : :class:`numpy.ndarray`
            Array of sample counts (2-D array - positions by number of samples)
        """
        cdef int i, j, n = x.shape[0]
        cdef data_type_t r, mu, k

        cdef long [:,:] sampled_vals = np.zeros((n, times), dtype = np.int_, order = 'c')
        cdef data_type_t [:,:] sampled_pvals = np.ones((n, times), dtype = np.float64, order = 'c')

        cdef long [:] vals

        for i in range(n):
            r = self.fit_r(x[i])
            mu = self.fit_mu(x[i])
        
            vals = np.random.negative_binomial(r, r/(r+mu), times)
            
            for j in range(times):
                sampled_vals[i, j] = vals[j]
                sampled_pvals[i, j] = nbinom.cdf(<int>vals[j], r/(r+mu), r)

        return sampled_vals, sampled_pvals

def learn_dispersion_model(h, cutoff = 250, trim = (2.5, 97.5)):
    """Learn a dispersion model from the expected 
    vs. observed histogram

    Parameters
    ----------
    h : :class:`numpy.ndarray`
        A 2-dimemsional array containing the distribution of
        observerd cleavages at each expected cleavage rate
    cutoff : int
        Mininum number of observed cleavages to perform ML 
        negative binomial fit at each value of expected cleavages
    trim : tuple (float)
        Percent of data to trim  from the observed cleavage count 
        (to mitigate outlier effects)

    Returns
    -------
    model : :class:`dispersion_model`
        A dispersion model learned from observed and expected counts

    Todo
    ----
    Add exceptions for failure to fit, etc.
    """
    size = int(h.shape[0])
    p = np.zeros(size)
    r = np.zeros(size)

    logger.info("Analyzing {:,} cleavages".format(np.sum(h)))

    # Make an initial negative binomial fit
    for i in range(0, size):
    
        # Unpack histogram
        pos = 0
        x = np.zeros(int(np.sum(h[i,:])))
        for j in range(len(h[i,:])):
            num = int(h[i, j])
            x[pos:pos+num] = j
            pos += num

        # If more than 500k points downsample to
        # make curve-fitting tractable
        if len(x)>1e5:
            x=np.random.choice(x, size=int(1e5))
            x=np.sort(x)

        if len(x) >= cutoff:		
            # Find data points to trim
            lower = int(np.floor(x.shape[0]*(trim[0]/100.0)))
            upper = int(np.ceil(x.shape[0]*(trim[1]/100.0)))

            # Compute intial estimates of mean and variance
            mu = np.mean(x[lower:upper])
            var = np.var(x[lower:upper])

            est_r = (mu * mu) / (var - mu)
            if est_r <= 0.0: 
                est_r = 10.0
            est_p = est_r / (est_r + mu)
            
            # Fit negative binomial
            (p[i], r[i]) = nbinom.fit(x[lower:upper], p = est_p, r = est_r)

        else:
            # Not enough data points to fit - will exclude 
            # from further analysis
            p[i] = r[i] = np.nan

    # Back-compute the mean values from the negative binomial parameters
    mus = p*r/(1-p)

    # Hack to control for problematic fits
    r[r>200] = 200.

    # Extrapolate using polynomial fit
    x = np.arange(size)
    sele = np.isfinite(mus)

    first_x = np.min(x[sele])
    last_x = np.max(x[sele])*0.75

    # fit mu with a 3 segments
    fit_mu = pwlf.PiecewiseLinFit(x[sele], mus[sele])
    res = fit_mu.fit_with_breaks_force_points(np.linspace(first_x, last_x, 4), [0], [mus[0]])

    # fit r with 5 s segments
    fit_r = pwlf.PiecewiseLinFit(x[sele], 1.0/r[sele])
    res = optimize.minimize(fit_r.fit_with_breaks_opt, [3.0, 7.0, 15.0, 25.0])
    
    x0 = np.zeros(6)
    x0[0] = first_x
    x0[-1] = last_x
    x0[1:-1] = res.x

    res = fit_r.fit_with_breaks_force_points(x0, [1], [1.0/r[1]])
    
    # Create a dispersion model class
    res = dispersion_model()
    res.h = h
    res.p = p
    res.r = r

    res.mu_params = list(fit_mu.fit_breaks[1:]) + list(fit_mu.intercepts) + list(fit_mu.slopes)
    res.r_params = list(fit_r.fit_breaks[1:]) + list(fit_r.intercepts) + list(fit_r.slopes)
    
    return res

import base64

def base64encode(x):
    return [str(x.dtype), base64.b64encode(x), x.shape]

def base64decode(x):
    dtype = np.dtype(x[0])
    arr = np.frombuffer(base64.b64decode(x[1]), dtype)
    if len(x) > 2:
        return arr.reshape(x[2])
    return arr

def load_dispersion_model(filename):
    """Load a dispersion model encoded in JSON format

    Parameters
    ----------
    filename : str
        Path to JSON-format dispersion model

    Returns
    -------
    model : :class:`dispersion_model`
        A dispersion model loaded from file
    """
    import simplejson as json
    import urllib.request as request

    if filename.startswith('http'):
        file = request.urlopen(filename)
    else:
        file = open(filename, 'r')

    params = json.load(file)

    file.close()

    model = dispersion_model()
    model.mu_params = base64decode(params['mu_params'])
    model.r_params = base64decode(params['r_params'])

    if 'h' in params:
        model.h = base64decode(params['h'])
    if 'p' in params:
        model.p = base64decode(params['p'])
    if 'r' in params:
        model.r = base64decode(params['r'])
    if 'metadata' in params:
        model.metadata = params['metadata']

    return model

def write_dispersion_model(model):
    """Write a JSON format dispersion model

    Parameters
    ----------
    model : :class:`dispersion_model`
        An instance of dispersion_model

    Returns
    -------
    out : str
        JSON-formatted dump of dispersion model
    """
    import simplejson as json
    from datetime import datetime

    out = { 'mu_params': base64encode(np.asarray(model.mu_params, order = 'C')), 
            'r_params': base64encode(np.asarray(model.r_params, order = 'C')),
            'h': base64encode(np.asarray(model.h, order = 'C')),
            'p': base64encode(np.asarray(model.p, order = 'C')),
            'r': base64encode(np.asarray(model.r, order = 'C')),
            'metadata': f"Created with {footprint_tools.__name__} {footprint_tools.__version__} "
                         + f"on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
    }

    return json.dumps(out, indent = 4)

