# Copyright 2015 Jeff Vierstra

import pysam
import numpy as np

from collections import defaultdict

class ReadError(Exception):
	pass

class GenotypeError(Exception):
	pass

class bamfile(object):

	"""Class to access a BAM file (largely inspired/copied from Piper et al.)"""
	
	def __init__(self, filepath, min_qual = 1, remove_dups = False, remove_qcfail = True, chunksize = 1000, offset = (0, -1)):

		try:
			self.samfile = pysam.Samfile(filepath, "rb")
		except:

			raise IOError("Cannot open BAM file: %s" % filepath)

		self.offset = offset #-1 # a hack for the mis-aligned data from 2010
		self.CHUNK_SIZE = chunksize

		self.min_qual = min_qual
		self.remove_dups = remove_dups
		self.remove_qcfail = remove_qcfail

	def __validate_read(self, read):
		if self.remove_qcfail and read.is_qcfail:
			raise ReadError()
		if self.remove_dups and read.is_duplicate:
			raise ReadError()
		if read.mapping_quality < self.min_qual:
			raise ReadError()

		return read

	def __get_read_mate(self, read):
		fpos = self.samfile.tell()
		try:
			mate = self.samfile.mate(read)
		except ValueError:
			mate = None
		finally:
			self.samfile.seek(fpos)
		return mate

	def __read_pair_generator(self, chrom, start, end):
		
		read_dict = defaultdict(lambda: [None, None])

		for read in self.samfile.fetch(chrom, max(start-10, 0), end+10):

			try:

				self.__validate_read(read)

				# if single-end just pass the read
				if not read.is_paired:
					yield read, None 
				
				# we think it is paired
				else:

					if not read.is_proper_pair or read.is_secondary or read.is_supplementary:
						continue

					qname = read.query_name
					if qname not in read_dict:
						if read.is_read1:
							read_dict[qname][0] = read
						else:
							read_dict[qname][1] = read
					else:
						if read.is_read1:
							yield read, read_dict[qname][1]
						else:
							yield read_dict[qname][0], read
						del read_dict[qname]

			except ReadError, e:
				continue


	def __add_read(self, read, fw, rev):
		if read.is_reverse:
			a = int(read.reference_end)+self.offset[1]
			rev[a] = rev.get(a, 0.0) + 1.0
		else:
			a = int(read.reference_start)+self.offset[0]
			fw[a] = fw.get(a, 0.0) + 1.0

	# DEPRECATED
	# def __lookup(self, chrom, start, end, flip=False):

	# 	tmp_fw = {}
	# 	tmp_rev = {}

	# 	for read in self.samfile.fetch(chrom, max(start-10, 0), end+10):
			
	# 		try:
	# 			self.__validate_read(read)
	# 			self.__add_read(read, tmp_fw, tmp_rev)

	# 			#if read.is_reverse:
	# 			#	a = int(read.reference_end) + self.offset[1]
	# 			#	if a < end:	
	# 			#		tmp_rev[a] = tmp_rev.get(a, 0.0) + 1.0
	# 			#else:
	# 			#	a = int(read.reference_start)  + self.offset[0] # pysam is zero-based!
	# 			#	if a >= start:
	# 			#		tmp_fw[a] = tmp_fw.get(a, 0.0) + 1.0

	# 		except:
	# 			pass

	# 	fw_cutarray = np.array([tmp_fw.get(i, 0.0) for i in range(start, end)])
	# 	rev_cutarray = np.array([tmp_rev.get(i, 0.0) for i in range(start, end)])

	# 	return {
	# 		"+": rev_cutarray[::-1] if flip else fw_cutarray, 
	# 		"-": fw_cutarray[::-1] if flip else rev_cutarray
	# 		}

	def __lookup(self, chrom, start, end, flip=False):

		tmp_fw = {}
		tmp_rev = {}

		for read1, read2 in self.__read_pair_generator(chrom, max(start-10, 0), end+10):

			self.__add_read(read1, tmp_fw, tmp_rev)
			if read2:
				self.__add_read(read2, tmp_fw, tmp_rev)


		fw_cutarray = np.array([tmp_fw.get(i, 0.0) for i in range(start, end)])
		rev_cutarray = np.array([tmp_rev.get(i, 0.0) for i in range(start, end)])


		return {
 			"+": rev_cutarray[::-1] if flip else fw_cutarray, 
			"-": fw_cutarray[::-1] if flip else rev_cutarray
			}

	def __getitem__(self, interval):

		chrom = interval.chrom
		start = interval.start
		end = interval.end

		flip = True if interval.strand == '-' else False

#		ret = self.__lookup(chrom, start, end) if not self.use_cache else self.__lookup_from_cache(chrom, start, end)
		ret = self.__lookup(chrom, start, end, flip)
		return ret

	def __lookup_allelic(self, chrom, start, end, pos, ref, alt, flip=False):

		#visited_read_pairs = set()

		tmp_ref_fw = {}
		tmp_ref_rev = {}
		tmp_alt_fw = {}
		tmp_alt_rev = {}

		for read1, read2 in self.__read_pair_generator(chrom, max(start-10, 0), end+10):

			#try read1 for genotype
			try:
				
				offset = pos-read1.reference_start
				s = read1.query_sequence[offset]

				mismatches = int(read1.get_tag("XM", with_value_type=False))

				if s==ref and mismatches<=1:
					fw=tmp_ref_fw
					rev=tmp_ref_rev
				elif s==alt and mismatches<=2:
					fw=tmp_alt_fw
					rev=tmp_alt_rev
				else:
					raise GenotypeError()

				self.__add_read(read1, fw, rev)
				if read2:
					self.__add_read(read2, fw, rev)
				continue

			except GenotypeError, e:
				continue
			except IndexError, e:
				pass


			# if we make if here and a paired read is available lets try it
			if not read2:
				continue

			#try read2 for genotype
			try:
				
				offset = pos-read2.reference_start
				s = read2.query_sequence[offset]

				mismatches = int(read2.get_tag("XM", with_value_type=False))

				if s==ref and mismatches<=1:
					fw=tmp_ref_fw
					rev=tmp_ref_rev
				elif s==alt and mismatches<=2:
					fw=tmp_alt_fw
					rev=tmp_alt_rev
				else:
					raise GenotypeError()

				self.__add_read(read1, fw, rev)
				self.__add_read(read2, fw, rev)

				continue

			except GenotypeError, e:
				continue
			except IndexError, e:
				pass

		ref_fw_cutarray = np.array([tmp_ref_fw.get(i, 0.0) for i in range(start, end)])
		ref_rev_cutarray = np.array([tmp_ref_rev.get(i, 0.0) for i in range(start, end)])
		alt_fw_cutarray = np.array([tmp_alt_fw.get(i, 0.0) for i in range(start, end)])
		alt_rev_cutarray = np.array([tmp_alt_rev.get(i, 0.0) for i in range(start, end)])

		return {
			ref: {
				"+": ref_rev_cutarray[::-1] if flip else ref_fw_cutarray, 
				"-": ref_fw_cutarray[::-1] if flip else ref_rev_cutarray
				},
			alt: {
				"+": alt_rev_cutarray[::-1] if flip else alt_fw_cutarray, 
				"-": alt_fw_cutarray[::-1] if flip else alt_rev_cutarray
			}
		}



	def get_allelic_reads(self, interval, pos, ref, alt, flip=False):

		chrom = interval.chrom
		start = interval.start
		end = interval.end

		flip = True if interval.strand == '-' else False

		ret = self.__lookup_allelic(chrom, start, end, pos, ref, alt, flip)
		return ret
