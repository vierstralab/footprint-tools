# Copyright 2015 Jeff Vierstra

import pysam
import numpy as np

from collections import defaultdict

class ReadError(Exception):
	pass

class GenotypeError(Exception):
	pass

class bamfile(object):
	"""Class to access a BAM file (largely inspired/copied from Piper et al.) 

		:param filepath: File path to BAM alignment file (must contain associated inddex)
		:type filepath: str
		:param min_qual: Filter reads by minimim mapping quality (MAPQ)
		:type min_qual: int
		:param remove_dups: Remove reads with duplicate flag (512) set
		:type remove_dups: bool
		:param remove_qcfail: Remove reads with QC fail flag (1024) set
		:type remove_qcfail: bool
		:param offset: Position offsets to apply to the `+` and `-` strands (default =(0, -1))
		:type filepath: tuple
	"""
	def __init__(self, filepath, min_qual = 1, remove_dups = False, remove_qcfail = True, offset = (0, -1)):
		"""Constructor"""
		try:
			self.samfile = pysam.Samfile(filepath, "rb")
		except:

			raise IOError("Cannot open BAM file: %s" % filepath)

		self.offset = offset #-1 # a hack for the mis-aligned data from 2010


		self.min_qual = min_qual
		self.remove_dups = remove_dups
		self.remove_qcfail = remove_qcfail

	def validate_read(self, read):
		"""
		Validate BAM tag

		:param read: Alignment record from BAM/SAM/CRAM file
		:type read: pysam.AlignedSegment 
		
		:return: Validated read 
		:rtype: :class:`pysam.AlignedSegment` or `NoneType`

		:raises ReadError: If read does not pass 1) QC flag, 2) duplicate, or 3) minimum mapping quality	
		"""
		if self.remove_qcfail and read.is_qcfail:
			raise ReadError()
		if self.remove_dups and read.is_duplicate:
			raise ReadError()
		if read.mapping_quality < self.min_qual:
			raise ReadError()

		return read

	def get_read_mate(self, read):
		"""
		Fetch the mate pair for paired-end reads

		:param read: Alignment record from BAM/SAM/CRAM file
		:type read: pysam.AlignedSegment
		:return: Mate pair read
		:rtype: :class:`pysam.AlignedSegment` or `NoneType`

		"""
		fpos = self.samfile.tell()
		try:
			mate = self.samfile.mate(read)
		except ValueError:
			mate = None
		finally:
			self.samfile.seek(fpos)
		return mate

	def read_pair_generator(self, chrom, start, end):
		"""
		Generator function that returns sequencing tags within a given region

		:param chrom: Chromosome (contig
		:type chrom:  str
		:param start: Start coordinate
		:type start: int
		:param end: End coordinate
		:type end: int
	
		:return: A tuple of :class:`pysam.AlignedSegments`. Elements may be `NoneType` if single-end sequencing or one of pairs falls outisde of query range.
		:rtype: tuple

		"""
		read_dict = defaultdict(lambda: [None, None])

		for read in self.samfile.fetch(chrom, max(start-10, 0), end+10):

			try:

				self.validate_read(read)

				# Single-end:,just pass the read
				if not read.is_paired:
					yield read, None 
				
				#  Pair-end; do some validation
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

			except ReadError as e:
				continue

		""" Flush out the rest of dictionary. (for example if one the 
		mates wasn't in the original region). It might be reasonable to 
		manually grab the remaining read using `__get_read_mate`"""
		for k, reads in read_dict.items():
			yield reads[0], reads[1]


	def __add_read(self, read, fw, rev):
		"""

		"""
		if read.is_reverse:
			a = int(read.reference_end)+self.offset[1]
			rev[a] = rev.get(a, 0.0) + 1.0
		else:
			a = int(read.reference_start)+self.offset[0]
			fw[a] = fw.get(a, 0.0) + 1.0

	def __lookup(self, interval):
		"""

		"""	
		chrom = interval.chrom
		start = interval.start
		end = interval.end
		flip = True if interval.strand == '-' else False

		tmp_fw = {}
		tmp_rev = {}

		for read1, read2 in self.read_pair_generator(chrom, max(start-10, 0), end+10):
			if read1:
				self.__add_read(read1, tmp_fw, tmp_rev)
			if read2:
				self.__add_read(read2, tmp_fw, tmp_rev)

		fw_cutarray = np.array([tmp_fw.get(i, 0.0) for i in range(start, end)])
		rev_cutarray = np.array([tmp_rev.get(i, 0.0) for i in range(start, end)])

		return {
 			"+": rev_cutarray[::-1] if flip else fw_cutarray, 
			"-": fw_cutarray[::-1] if flip else rev_cutarray
			}

	def __validate_genotype(self, read, pos, ref, alt):
		"""

		Returns
		-------
		bool
			True if read contains reference allele and passes filters
			False if read contains alternate allele and passes filters

		Raises
		------
		IndexError
		"""
		if not read:
			return 0

		var_offset=pos-read.reference_start

		try:
			base_call=read.query_sequence[var_offset]
		except IndexError as e:
			return 0

		mismatches=int(read.get_tag("XM", with_value_type=False))
		if base_call==ref and mismatches<=1:
			return ref
		if  base_call==alt and mismatches<=2:
			return alt
			raise GenotypeError()



		return read.flag & (1<<12)

	def __lookup_allelic(self, chrom, start, end, pos, ref, alt, flip=False):
		"""

		"""
		#visited_read_pairs = set()

		tmp_ref_fw = {}
		tmp_ref_rev = {}
		tmp_alt_fw = {}
		tmp_alt_rev = {}

		for read1, read2 in self.read_pair_generator(chrom, max(start-10, 0), end+10):

			#try read1 for genotype
			try:

				# Check reads for proper genotypes; raise Genotype error
				# if some thing is problematic, if a read pair doesn't exist or
				# doesn't overlap the SNV then returns 0
				read1_gt=self.__validate_genotype(read1)
				read2_gt=self.__validate_genotype(read2)

				if (not read1_gt) and (not read2_gt): # neither read overlaps SNV
					continue # got to get read pair
				elif (read1_gt and read2_gt) and read1_gt!=read2_gt: # discordant genotypes
					raise GenotypeError("Discordant genotypes between mate reads!")
				elif read1_gt==ref or read2_gt==ref: # Matches REF allele
					fw=tmp_ref_fw
					rev=tmp_ref_rev
				elif read1_gt==alt or read2_gt==alt: # Matches ALT allele
					fw=tmp_alt_fw
					rev=tmp_alt_rev
				else:
					raise GenotypeError()

				if read1:
					self.__add_read(read1, fw, rev)
				if read2:
					self.__add_read(read2, fw, rev)				

			except GenotypeError as e:
				continue
	
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


	def __getitem__(self, x):
		"""
		
		:param x: genomic positions to overlap reads
		:type x: genomic_interval or pysam.VariantRecord 
		:return: Dictionary of :class:`numpy.array` corresponding to  counts of the 5` ends of reads mapping to the `+` and `-` strands
		:rtype: dict
		"""

		if isinstance(x, genomic_interval):
			return self.__lookup(x)
		elif isinstance(x, VariantRecord):
			return self.__lookup_allelic(x)
		else:
			raise TypeError()
