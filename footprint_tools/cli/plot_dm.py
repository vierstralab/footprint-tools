import click
from click_option_group import optgroup

import math

from matplotlib.pylab import rcParams
import matplotlib.pyplot as plt
import matplotlib.gridspec as mgridspec

from footprint_tools.modeling import dispersion
from footprint_tools.plotting import (plot_model_fit, plot_model_histogram)

from footprint_tools.cli.utils import list_args

import logging
logger = logging.getLogger(__name__)

@click.command(name='plot_dm')
@click.argument('dispersion_model_file')
@optgroup.group('Output options')
@optgroup.option('--histograms',
    type=click.STRING, default="5,25,50,75", callback=list_args(int),
    help='Plot histograms of observed counts at site with expected counts (comma-seperated list)')
@optgroup.option('--outfile',
    type=click.STRING, default='dm.pdf',
    help='Output file path for plot (suffix determines image format)')
def run(dispersion_model_file, histograms=[15,25,50,75], outfile='dm.pdf'):
    """Diagnostic plotting of a dispersion model
    
    Using the ``--histogram`` argument specifies the histograms to plot. An arbitry 
    number of histograms can be plotted.

    Outputs a PDF with plots
    """
    try:
        dm = dispersion.load_dispersion_model(dispersion_model_file)
    except IOError as e:
        logger.critical(e)
        return 1

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
    gs = mgridspec.GridSpec(nrows, ncols,  wspace=0.75, hspace=0.75)

    logger.info("Plotting model parameters")

    ax0 = fig.add_subplot(gs[0,0])
    ax1 = fig.add_subplot(gs[0,1])
    
    plot_model_fit(dm, axs=(ax0,ax1))

    logger.info(f"Plotting histograms - {histograms}")

    for i, n in enumerate(histograms):
        row_index = (i // ncols) + 1
        col_index = i % ncols
        ax = fig.add_subplot(gs[row_index, col_index])
        plot_model_histogram(dm, n=n, ax=ax)

    fig.set_size_inches(2.5*ncols, 2*nrows)
    
    logger.info(f"Saving plots to {outfile}")	

    plt.savefig(outfile, transparent=True)

    return 0
