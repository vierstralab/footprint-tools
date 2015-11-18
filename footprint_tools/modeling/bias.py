'''
JDV Aug.2015

footprint toolkit/bias_model

Class used to perform various shufflings of cleavage data
based on sequence within a local window. A hexamer (6mer)
model is currently hardcoded, but should change in the 
future to accomodate other sizes.
'''

import numpy as np
import itertools

import random

class bias_model(object):

	def __init__(self):
		self.model = {}
		self.offset = 0

	def get_value(self, key):
		try:
			return self.model[key]
		except:
			return 1e-6

	def shuffle(self):
		"""Randomly shuffle the bias model
		"""
		ret = bias_model()
		ret.model = { x: y for (x, y) in zip(self.model.keys(), sorted(self.model.values(), key = lambda k: random.random())) }
		ret.offset = self.offset
		return ret

	def predict(self, probs, n = 100):
		"""Compute cleavage propensities from sequence
		Args:

		Return:
		"""
		# Random sampling
		#x = np.random.choice(len(probs), size = n, p = probs/np.sum(probs))
		#return np.bincount(x, minlength = len(probs))
		
		# Deterministic
		return np.around(probs/np.sum(probs) * n)

class kmer_model(bias_model):

	def __init__(self, filepath):
		bias_model.__init__(self)
		self.offset = 3
		self.k = 6
		self.read_model(filepath)
	
	def read_model(self, filepath):
		try:
			for line in open(filepath, 'r'):
				(seq, prob) = line.strip().split('\t')
				self.model[seq.upper()] = float(prob)
		except IOError:
			raise IOError("Cannot open file: %s" % filepath)
	
	def probs(self, seq):
		return np.array([ self.get_value( seq[(i-self.offset):(i+(self.k-self.offset))] ) for i in range(self.offset, len(seq)-(self.k-self.offset)) ])

class uniform_model(bias_model):

	def __init__(self):
		bias_model.__init__(self)
		for seq in itertools.product('ATCG', repeat=6):
			self.model[''.join(seq)] = 1.0

	def probs(self, seq):
		return np.ones( len(seq) )
