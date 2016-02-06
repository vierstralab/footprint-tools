# Copyright 2015 Jeff Vierstra

import numpy as np

def step(arr, xaxis = False, interval = 0):
	"""
	"""

	if xaxis and interval == 0:
		interval = abs(arr[1] - arr[0]) / 2.0
	newarr = np.array(zip(arr - interval, arr + interval)).ravel()
	return newarr


def step_plot(arr, ax, **kwargs):
	ax.fill_between(step(np.arange(len(arr)), xaxis = True), step(np.zeros(len(arr))), step(arr), **kwargs)


def gradient_step_plot(arr, lim, ax, **kwargs):

	from matplotlib.path import Path
	from matplotlib.patches import PathPatch

	x = step(np.arange(len(arr)), xaxis = True)
	y = step(arr)

	p = PathPatch(Path(np.array([x, y]).transpose()), facecolor = 'none', lw = 0)
	ax.add_patch(p)

	im = ax.imshow(x.reshape(y.size, 1), interpolation = "bicubic", origin = "lower", cmap = "Purples",
					aspect = "auto", clip_path = p, clip_on = True,  extent = [0, len(arr), lim[0], lim[-1]], **kwargs)

def segments_plot(arr, thresholds, w, labels, ax, **kwargs):

	from ..stats import segment
	from matplotlib.patches import Rectangle

	levels = [segment(arr, thresh, w) for thresh in thresholds]

	for i in np.arange(len(levels)):
		for s in levels[i]:
			p = Rectangle((s[0], i-0.35), s[1]-s[0], 0.7, **kwargs)
			ax.add_patch(p)

	ax.set_yticks(np.arange(len(levels)))
	ax.set_yticklabels(labels)
	ax.set_ylim([-1, len(levels)])

