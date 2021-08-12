import math

import click

import numpy as np
import scipy.stats

from matplotlib.pylab import rcParams
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MaxNLocator

from footprint_tools.modeling import dispersion

from footprint_tools.cli.utils import list_args

import logging
logger = logging.getLogger(__name__)

def plot_model_mu(dm, ax=None, xlim=(0, 100)):
    """
    Plot model mu parameters
    """
    x = np.arange(xlim[0], xlim[1])

    # Raw parameters
    r = np.array([dm.r[i] for i in x])
    p = np.array([dm.p[i] for i in x])
    mu = p*r/(1.0-p)

    # Smoothed parameters
    fit_mu = np.array([dm.fit_mu(i) for i in x])

    # Plot functions & appearance
    ax.plot(x, mu, label='MLE neg. binomial fit')
    ax.plot(x, fit_mu, label='Smoothed fit', ls='dashed')
    ax.plot(xlim, xlim, label='y=x', color='grey', ls='dashed', zorder=-10)

    ax.set_xlabel('Expected cleavage count')
    ax.set_ylabel('Observed cleavages')

    [ax.spines[loc].set_color('none') for loc in ['top', 'right']]

    ax.xaxis.set_ticks_position('bottom')
    ax.xaxis.set_tick_params(direction='out')
    ax.xaxis.set(major_locator = MaxNLocator(4))

    ax.yaxis.set_ticks_position('left')
    ax.yaxis.set_tick_params(direction = 'out')
    ax.yaxis.set(major_locator = MaxNLocator(4))

    ax.legend()

def plot_model_r(dm, ax=None, xlim=(1, 100)):
    """
    Plot model dispersion parameters
    """
    x = np.arange(xlim[0], xlim[1])

    # Raw parameters
    r = np.array([dm.r[i] for i in x])

    # Smoothed parameters
    fit_r = np.array([dm.fit_r(i) for i in x])

    ax.plot(x, 1/r, label='MLE neg. binomial fit')
    ax.plot(x, 1/fit_r, label='Smoothed fit', ls='dashed')

    ax.set_xlabel("Expected cleavage count")
    ax.set_ylabel("1/r")

    [ax.spines[loc].set_color('none') for loc in ['top', 'right']]

    ax.xaxis.set_ticks_position('bottom')
    ax.xaxis.set_tick_params(direction='out')
    ax.xaxis.set(major_locator = MaxNLocator(4))

    ax.yaxis.set_ticks_position('left')
    ax.yaxis.set_tick_params(direction = 'out')
    ax.yaxis.set(major_locator = MaxNLocator(4))

    ax.legend()

def plot_histogram(dm, n=25, show_poisson=True, ax=None, xlim=(0, 125)):
    """
    Plots a density histogram of the observed cleavage counts
    at an expected cleavage rate (n).
    """
    x = np.arange(xlim[0], xlim[1])

    mu = dm.fit_mu(n)
    r = dm.fit_r(n)

    # Raw observed counts
    ax.bar(x, dm.h[n,x[0]:(x[-1]+1)]/np.sum(dm.h[n,:]), width=1, color='lightgrey', label="Observed")

    # NB fit
    y_nbinom=scipy.stats.nbinom.pmf(x, r, r/(r+mu))
    ax.plot(x, y_nbinom, color="red", label="Negative binomial")

    # Poisson
    if show_poisson:
        y_pois=scipy.stats.poisson.pmf(x, mu=n)
        ax.plot(x, y_pois, color="blue", label="Poisson")

    ax.set_xlim(x[0], x[-1])

    [ax.spines[loc].set_visible(False) for loc in ["top", "right"]]
    ax.set_xlabel("Observed cleavage count")
    ax.set_ylabel("Density")

    ax.set_title(f"{n} expected cleavages")

    ax.legend()

@click.command(name='plot_dm')
@click.argument('dispersion_model_file')
@click.option('--histograms',
    type=click.STRING, default="5,25,50,75", callback=list_args(int),
    help='Plot histograms of observed counts at site with expected counts (comma-seperated list)')
@click.option('--outfile',
    type=click.STRING, default='dm.pdf',
    help='Output file path for plot (suffix determines image format)')
def run(dispersion_model_file, histograms=[15,25,50,75], outfile='dm.pdf'):
    """Diagnostic plotting of a dispersion model
    
    Output:
        dm.pdf - a PDF file with plots
    """
    
    dm = dispersion.load_dispersion_model(dispersion_model_file)

    plt_params = {
        'legend.fontsize': 'xx-small',
        'axes.labelsize': 'small',
        'axes.titlesize': 'small',
        'xtick.labelsize': 'x-small',
        'ytick.labelsize': 'x-small'}
    rcParams.update(plt_params)

    npanels = len(histograms)+2
    ncols = 2
    nrows = math.ceil(npanels/ncols)

    fig = plt.figure()
    gs = gridspec.GridSpec(nrows, ncols,  wspace=0.75, hspace=0.75)

    logger.info("Plotting model parameters")

    ax = fig.add_subplot(gs[0,0])
    plot_model_mu(dm, ax)

    ax = fig.add_subplot(gs[0,1])
    plot_model_r(dm, ax)

    logger.info(f"Plotting histograms - {histograms}")

    for i, n in enumerate(histograms):
        row_index = (i // ncols) + 1
        col_index = i % ncols
        ax = fig.add_subplot(gs[row_index, col_index])
        plot_histogram(dm, n=n, ax=ax)

    fig.set_size_inches(2.5*ncols, 2*nrows)
    
    logger.info(f"Saving plots to {outfile}")	

    plt.savefig(outfile, transparent=True)

    return 0
