
''' '''
import numpy as np

class genomic_interval(object):

	def __init__(self, chrom, start, end, name = 0, score = 0, strand = '+'):

		self.chrom = str(chrom)
		self.start = int(start)
		self.end = int(end)
		self.name = str(name)
		self.score = float(score)
		self.strand = str(strand)

	def __len__(self):

		return self.end - self.start
	
	def __str__(self):
	
		return '\t'.join( [ str(x) for x in [self.chrom, self.start, self.end, self.name, self.score, self.strand] ] )

	def widen(self, w):

		return genomic_interval(self.chrom, self.start - w, self.end + w)

''' '''

class genomic_interval_set(object):

	def __init__(self, filepath = None):

		self.intervals = []

		if filepath:

			self.__read_bed_file(filepath)
	
	def __len__(self):

		return len(self.intervals)

	def __iter__(self):

		for x in self.intervals:

			yield x

	def __getitem__(self, i):

		return self.intervals[i]

	def add(self, other):

		if type(other) == genomic_interval_set:
			self.intervals.extend(other.intervals)
		else:
			self.intervals.append(other)

	def __read_bed_file(self, filepath):

		try:

			file = open(filepath, "r")

		except:

			raise IOError("Cannot open file: %s" % filepath)

		for line in file:

			fields = line.strip().split('\t')

			if len(fields) < 6:

				raise IOError("Malformed BED entry!")

			self.add( genomic_interval(fields[0], fields[1], fields[2], fields[3], fields[4], fields[5]) )	
