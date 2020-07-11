# Copyright 2015 Jeff Vierstra

import genome_tools
import pysam

import numpy as np
from collections import defaultdict

class ReadError(Exception):
	pass

class GenotypeError(Exception):
	pass

class bamfile(object):
	"""Class to access BAM files
	
	Attributes
	----------
	min_qual : int
	    Filter reads by minimim mapping quality (MAPQ)
	offset : tuple
	    Position offsets to apply to the `+` and `-` strands (default =(0, -1))
	remove_dups : bool
	    Remove reads with duplicate flag (512) set
	remove_qcfail : bool
	    Remove reads with QC fail flag (1024) set
	samfile : pysam.Samfile
	    SAM/BAM file object
	"""
	
	def __init__(self, filepath, min_qual = 1, remove_dups = False, remove_qcfail = True, offset = (0, -1)):
		"""Constructor
		
		Parameters
		----------
		filepath : str
		    Fileath to SAM/BAM file
		min_qual : int, optional
		    Filter reads by minimim mapping quality (MAPQ)
		remove_dups : bool, optional
		    Remove reads with duplicate flag (512) set
		remove_qcfail : bool, optional
		    Remove reads with QC fail flag (1024) set
		offset : tuple, optional
		    Position offsets to apply to the `+` and `-` strands (default =(0, -1))
		
		Raises
		------
		IOError
		    If file is not found, has no index, or something else bad happened while trying to open BAM file
		"""

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
		
		Parameters
		----------
		read : pysam.AlignedSegment
		    Read from BAM/SAM file
		
		Returns
		-------
		pysam.AlignedSegment
		    Same read as input
		
		Raises
		------
		ReadError
		    Raises error if read fails QC flag, is a duplicate or MAPQ < minimum
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
		
		Parameters
		----------
		read : pysam.AlignedSegment
			One  of the read-pairs
		
		Returns
		-------
		pysam.AlignedSegment
		    The corresponding mate read
		
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
		
		Parameters
		----------
		chrom : str
		    Chromosome
		start : int
		    Start coordinate
		end : int
		    End coordinate
		
		Yields
		------
		tuple
		    A tuple of :class:`pysam.AlignedSegments`. Elements may be `NoneType` if single-end sequencing or one of pairs falls outisde of query range.
		
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

		"""Flush out the rest of dictionary. (for example if one the 
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
		"""Getter function to retrieve cutcounts
		
		Parameters
		----------
		x : genome_tools.genomic_interval or pysam.VariantRecord
		    Retrieve cleavages over a windowed region or resolve allelically over a known variant
		
		Returns
		-------
		dict
		    A dictionary of read counts resolved to each strand, and each allele (if approriate)
		
		Raises
		------
		TypeError
		    If input is neither a genome_tools.genomic_interval nor a pysam.VariantRecord
		"""
		if isinstance(x, genome_tools.genomic_interval):
			return self.__lookup(x)
		elif isinstance(x, pysam.VariantRecord):
			return self.__lookup_allelic(x)
		else:
			raise TypeError()
