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
		shuffled = sorted(self.model.values(), key=lambda k: random.random())
		
		ret = bias_model()
		ret.model = { x: y for (x, y) in zip(self.model.keys(), shuffled) }
		ret.offset = self.offset

	# Note: fw and rv are already oriented to the strand (see cutcounts.py)
	# Also assumed that the sequence is correctly oriented (reverse complement)
	# for negative strand
	
	# New more generalizable function

	def predict(self, probs, n = 100):
		
		#res = np.zeros(len(probs))
		#for i in np.random.choice(len(probs), size = n, p = probs/np.sum(probs)):
		#	res[i] += 1.0
		#return res

		#x = np.random.choice(len(probs), size = n, p = probs/np.sum(probs))
		#return np.bincount(x, minlength = len(probs))

		return np.around(probs/np.sum(probs) * n)

class kmer_model(bias_model):

	def __init__(self, filepath):

		bias_model.__init__(self)

		self.offset = 3
		self.k = 6

		self.read_model(filepath)
	
	def read_model(self, filepath):

		try:
			file = open(filepath, 'r')
		except IOError:
			raise IOError("Cannot open file: %s" % filepath)

		for line in file:
			(seq, prob) = line.strip().split('\t')
			self.model[seq.upper()] = float(prob)
	

	# returns an array of probabilities values that is smaller the the sequence
	def probs(self, seq):

		return np.array([ self.get_value( seq[(i-self.offset):(i+(self.k-self.offset))] ) for i in range(self.offset, len(seq)-(self.k-self.offset)) ])

		#p = np.zeros(len(seq))

		#for i in range(self._offset, len(seq) - self._k - self.offset + 1):
		#	try:
		#		p[i] = self.model[ seq[(i-self._offset):(i+self._k-self._offset)] ]
		#	except KeyError:
		#		p[i] = 1e-6

		#return p

class uniform_model(bias_model):

	def __init__(self):

		bias_model.__init__(self)

		for seq in itertools.product('ATCG', repeat=6):
			self.model[''.join(seq)] = 1.0

	def probs(self, seq):

		probs = np.ones( len(seq) )

		return probs
