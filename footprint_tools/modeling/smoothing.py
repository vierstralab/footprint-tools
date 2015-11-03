import numpy as np

class moving_trimmed_mean(object):

	def __init__(self, half_window_width = 50, trim = [0.05, 0.95]):

		self.half_window_width = half_window_width

		lower = np.floor(trim[0] * self.half_window_width * 2)
		upper = np.ceil(trim[1] * self.half_window_width * 2)

		self.trimmed_mean_func = lambda x: np.rint(np.mean(np.sort(x)[lower:upper]))

	def smooth(self, counts, half_window_width):

		idx = np.vstack( np.arange(i-self.half_window_width, i+self.half_window_width+1) for i in np.arange(self.half_window_width, len(counts)-self.half_window_width) )

		return np.apply_along_axis(self.trimmed_mean_func, 1, counts[idx])
		
		#if i - self._smoothing_half_window_width < 0:
		#	res = 0
		#elif i + self._smoothing_half_window_width > len(self._window_counts):
		#	res = 0
		#else:
		#	limits = np.percentile(self._window_counts[i-self._smoothing_half_window_width:i+self._smoothing_half_window_width+1], self._trim)
		#	res = scipy.stats.tmean(self._window_counts[i-self._smoothing_half_window_width:i+self._smoothing_half_window_width+1], limits = limits)

		#return res
