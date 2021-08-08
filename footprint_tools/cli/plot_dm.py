import sys
import os
import math

import argh
from argh.decorators import named, arg

import numpy as np
import scipy.stats

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MaxNLocator

from footprint_tools.modeling import dispersion

from footprint_tools.cli.utils import list_ints

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
	fit_r = np.array([dm.fit_r(i) for i in x])

	# Plot functions & appearance
	ax.plot(x, mu, label='MLE neg. binomial fit')
	ax.plot(x, fit_mu, label='Smoothed fit', ls='dashed')
	ax.plot(xlim, xlim, label='y=x', color='grey', ls='dashed', zorder=-10)

	ax.set_xlabel('Expected cleavage count')
	ax.set_ylabel('Observed cleavages (mean)')

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
	ax.plot(x, 1/fit_r, label='Smooth fit', ls='dashed')

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
	ax.bar(x, dm.h[n,x[0]:x[-1]]/np.sum(dm.h[n,:]), width=1, color='lightgrey', label="Observed")

	# NB fit
	y_nbinom=scipy.stats.nbinom.pmf(x, r, r/(r+mu))
	ax.plot(x, y_nbinom, color="red", label="Negative binomial")

	# Poisson
	if show_poisson:
		y_pois=scipy.stats.poisson.pmf(x, mu=n)
		ax.plot(x, y_pois, color="blue", label="Poisson" % n)

	ax.set_xlim(x[0], x[-1])

	[ax.spines[loc].set_visible(False) for loc in ["top", "right"]]
	ax.set_xlabel("Observed DNase I cleavage counts")
	ax.set_ylabel("Density")

	ax.set_title("Observed cleavage counts at positions with %d expected cleavages" % n)

	ax.legend()

@named('plot_dm')
@arg('dispersion_model_file',
	type=str,
	help='Dispersion model file (can be a remote URL -- http protocol)')
@arg('--histograms',
	type=list_ints,
	default=[15,25,50,75],
	help='')
def run(dispersion_model_file, histograms=[15,25,50,75]):
	"""
	Diagnostic plotting of a dispersion model
	"""
	
	dm = dispersion.read_dispersion_model(dispersion_model_file)

	npanels = len(histograms)+2
	
	ncols = 2
	nrows = math.ceil(npanels/ncols)

	fig = plt.figure()
	gs = gridspec.GridSpec(nrows, ncols)

	ax = fig.add_subplot(gs[0,0])
	plot_model_mu(dm, ax)

	ax = fig.add_subplot(gs[0,1])
	plot_model_r(dm, ax)

	outfile = os.path.abspath(os.path.join(os.getcwd(), 'dm.pdf'))
	plt.savefig(outfile, transparent=True)

	return 0