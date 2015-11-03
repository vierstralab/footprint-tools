'''
JDV Aug.2015

footprint toolkit/cutcouts

Class used to get per-nucleotide cleavage counts directly from
a BAM alignment file using the pysam library. Largely inspired
by Piper et al. Wellington algorithm.
'''

import pysam
import numpy as np

import genomic_interval

class bamfile(object):

	def __init__(self, filepath, chunksize = 1000):

		try:
			self.samfile = pysam.Samfile(filepath)
		except:

			raise IOError("Cannot open BAM file: %s" % filepath)

		self.cache = { i: {"+": {}, "-": {}} for i in self.samfile.references }
		self.lookup = { i: [] for i in self.samfile.references }

		self.offset = -1 # a hack for the mis-aligned data from 2010
		self.CHUNK_SIZE = chunksize

	def __add(self, chrom, start, end):
        
		for alignedread in self.samfile.fetch(chrom, max(start, 0), end):

				if alignedread.mapq < 1:
					continue

				if alignedread.is_reverse:
					a = int(alignedread.aend) - 1 - self.offset
					if a < end:	
						self.cache[chrom]["-"][a] = self.cache[chrom]["-"].get(a, 0.0) + 1.0
				else:
					a = int(alignedread.pos) + self.offset # pysam is zero-based!
					if a >= start:
						self.cache[chrom]["+"][a] = self.cache[chrom]["+"].get(a, 0.0) + 1.0
		
		self.lookup[chrom].append(start)

	def __lookup(self, chrom, start, end):

		lbound = int(np.floor(start / float(self.CHUNK_SIZE)) * float(self.CHUNK_SIZE))
		ubound = int(np.ceil(end / float(self.CHUNK_SIZE)) * float(self.CHUNK_SIZE))

		for i in range(lbound, ubound, self.CHUNK_SIZE):
           		
			if i not in self.lookup[chrom]:
                		
				self.__add(chrom, i, i + self.CHUNK_SIZE)
        
		#fills in with zeroes where the hash table contains no information for each strand.
		fw_cutarray  = np.array([self.cache[chrom]["+"].get(i, 0.0) for i in range(start, end)])
		rev_cutarray = np.array([self.cache[chrom]["-"].get(i, 0.0) for i in range(start, end)])

		return {"+": fw_cutarray, "-": rev_cutarray}
	
	def __getitem__(self, interval):

		chrom = interval.chrom
		start = interval.start
		end = interval.end

		flip = True if interval.strand == '-' else False
		
		ret = self.__lookup(chrom, start, end)

		if flip:
			ret["+"], ret["-"] = ret["-"][::-1], ret["+"][::-1]

		return ret

