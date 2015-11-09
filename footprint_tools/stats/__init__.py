from .distributions import *

import scipy.stats

def windowed_chi_squared(x, w = 3):
    ''' Compute chi-squared values from groups of ln p-values
    Args:
    	x (array): natural log p-values
    	w (int): window size to smooth
    Returns:
    	(chi, p): tuple of chi2 values and p-values
    '''
    
    chi = np.zeros(len(x))
    p = np.zeros(len(x))
    for i in np.arange(w, len(x)-w+1):
        chi[i] = -2 * np.sum(x[i-w:i+w+1])
        p[i] = scipy.stats.chi2.logsf(chi[i], w*2)

    return (chi, p)


