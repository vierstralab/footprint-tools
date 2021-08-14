import numpy as np

import matplotlib.pyplot as mpl
import matplotlib.ticker as mtick

def set_spines(ax, remove=['top', 'right']):
    """Clear spines"""
    [ax.spines[loc].set_color('none') for loc in remove]

def format_axis_default(ax):
    """Format axis with defaults"""
    ax.xaxis.set_ticks_position('bottom')
    ax.xaxis.set_tick_params(direction='out')
    ax.xaxis.set(major_locator = mtick.MaxNLocator(4))

    ax.yaxis.set_ticks_position('left')
    ax.yaxis.set_tick_params(direction='out')
    ax.yaxis.set(major_locator = mtick.MaxNLocator(4))

def plot_model_fit(dm, axs=None, xlim=(0, 100)):
    """Plot model mu parameters

    Parameters
    ----------
    dm : :class:`dispersion_model`
        A dispersion model
    ax : :class:`matplotlib.Axes` or NoneType
        Plot axes
    xlim : tuple
    """
    if not axs:
        fig, axs = mpl.subplots(nrows=1, ncols=2, sharex=True, sharey=False)
    else:
        fig = mpl.gcf()

    x = np.arange(max(1, xlim[0]), xlim[1])

    # Raw parameters
    r = np.array([dm.r[i] for i in x])
    p = np.array([dm.p[i] for i in x])
    mu = p*r/(1.0-p)

    # Smoothed parameters
    fit_mu = np.array([dm.fit_mu(i) for i in x])
    fit_r = np.array([dm.fit_r(i) for i in x])

    # Plot functions & appearance
    set_spines(axs[0])
    format_axis_default(axs[0])

    axs[0].plot(x, mu, label='MLE neg. binomial fit')
    axs[0].plot(x, fit_mu, label='Smoothed fit', ls='dashed')
    axs[0].plot(xlim, xlim, label='y=x', color='grey', ls='dashed', zorder=-10)

    axs[0].set_xlabel('Expected cleavage count')
    axs[0].set_ylabel('Observed cleavages')

    axs[0].legend()

    # Plot functions & appearance
    set_spines(axs[1])
    format_axis_default(axs[1])

    axs[1].plot(x, 1/r, label='MLE neg. binomial fit')
    axs[1].plot(x, 1/fit_r, label='Smoothed fit', ls='dashed')

    axs[1].set_xlabel('Expected cleavage count')
    axs[1].set_ylabel('1/r')

    axs[1].legend()

    return fig, axs

def plot_model_histogram(dm, n=25, ax=None, xlim=(0, 125), show_poisson=True):
    """Plots a density histogram of the observed cleavage counts
    at an expected cleavage rate <n>.

    Parameters
    ---------
    dm : :class:`dispersion_model`
        A dispersion model
    ax : :class:'matplotlib.Axes`
        Plotting axes
    """

    from scipy.stats import nbinom, poisson

    if not ax:
        fig, ax = mpl.subplots()
    else:
        fig = mpl.gcf()

    x = np.arange(xlim[0], xlim[1])

    mu = dm.fit_mu(n)
    r = dm.fit_r(n)

    dens = dm.h[n,x[0]:(x[-1]+1)].copy()
    dens /= dm.h[n,:].sum()

    set_spines(ax)
    format_axis_default(ax)

    # Raw observed counts
    ax.bar(x, dens, width=1, color='lightgrey', label="Observed")

    # NB fit
    y = nbinom.pmf(x, r, r/(r+mu))
    ax.plot(x, y, color="red", label="Negative binomial")

    # Poisson
    if show_poisson:
        y = poisson.pmf(x, mu=n)
        ax.plot(x, y, color="blue", label="Poisson")

    ax.set_xlabel("Observed cleavage count")
    ax.set_ylabel("Density")

    ax.set_title(f"{n} expected cleavages")

    ax.legend()

    return fig, ax