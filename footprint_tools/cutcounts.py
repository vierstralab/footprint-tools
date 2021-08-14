"""
This modules contains classes and functions to compute cleavage counts
directly from an alignment file.
"""
import genome_tools
import pysam

import numpy as np
from collections import defaultdict

import logging
logger = logging.getLogger(__name__)

class ReadError(Exception):
    
    ERROR_ALIGNMENT = (0, "Read alignment problematic (QC fail, duplicate, or MAPQ < 1)")
    ERROR_5PROXIMITY = (1, "Variant too close to 5' end of tag") 
    ERROR_BASEQ = (2, "Base quality < 20")
    ERROR_GENOTYPE = (3, "Base does not match reference or expected alternate allele")
    ERROR_MISMATCH = (4, "Read contains too many mismatches")
    
    def __init__(self, e):
        self.value = e[0]
        self.message = e[1]


class GenotypeError(Exception):
    pass

class bamfile(object):
    """Class to access BAM files
    
    Attributes
    ----------
    min_qual : int
        Filter reads by minimim mapping quality (MAPQ)
    offset : tuple
        Position offsets to apply to the `+` and `-` strands,. DNase I data (0, -1).
        Tn5-derived data would use (4,-5). (default = (0, -1))
    remove_dups : bool
        Remove reads with duplicate flag (512) set
    remove_qcfail : bool
        Remove reads with QC fail flag (1024) set
    samfile : pysam.Samfile
        SAM/BAM file object
    """
    
    def __init__(self, 
                filepath, 
                min_qual = 1, 
                remove_dups = False,
                remove_qcfail = True,
                offset = (0, -1)):
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
            If file is not found
        ValueError
            If file has not associated index
        """

        try:
            self.samfile = pysam.Samfile(filepath, "rb")
        except:
            raise IOError("Cannot open BAM file: %s" % filepath)

        self.offset = offset #-1 # a hack for the mis-aligned data from 2010
        self.min_qual = min_qual
        self.remove_dups = remove_dups
        self.remove_qcfail = remove_qcfail
        
    def close(self):
        """Closes BAM file
        """
        if self.samfile:
            self.samfile.close()
        return True

    def __del__(self):
        return self.close()

    def validate_read(self, read):
        """Validate BAM tag
        
        Parameters
        ----------
        read : :class:`pysam.AlignedSegment`
            Read from BAM/SAM file
        
        Returns
        -------
        read: :class:`pysam.AlignedSegment`
            Same read as input
        
        Raises
        ------
        ReadError
            Raises error if read fails QC flag, is a duplicate or MAPQ < minimum
        """

        if self.remove_qcfail and read.is_qcfail:
            raise ReadError(ReadError.ERROR_ALIGNMENT)
        if self.remove_dups and read.is_duplicate:
            raise ReadError(ReadError.ERROR_ALIGNMENT)
        if read.mapping_quality < self.min_qual:
            raise ReadError(ReadError.ERROR_ALIGNMENT)

        return read

    def _get_read_mate(self, read):
        """Fetch the mate pair for paired-end reads
        
        Parameters
        ----------
        read : :class:`pysam.AlignedSegment`
            One  of the read-pairs
        
        Returns
        -------
        mate: :class:`pysam.AlignedSegment`
            The corresponding mate read two input read
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
        """Generator function that returns sequencing tags within a given region
        
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
        reads: tuple
            A tuple of :class:`pysam.AlignedSegment`. Elements may be NoneType if 
            single-end sequencing or one of pairs falls outisde of query range.
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
            
    def _add_read(self, read, fw, rev):
        """Add read to output dictionary

        Parameters
        ----------
        read: :class:`pysam.AlignedSegment`
            Aligned read segment to add
        fw: dict
            Dictionary of read counts on forward strand
        rev: dict
            Dictionary of read counts on reverse strand
        """
        if read.is_reverse:
            a = int(read.reference_end)+self.offset[1]
            rev[a] = rev.get(a, 0.0) + 1.0
        else:
            a = int(read.reference_start)+self.offset[0]
            fw[a] = fw.get(a, 0.0) + 1.0

    def _get_fragment(self, read):
        """Returns a fragment from a mapped read
        
        Parameters
        ----------
        read : class:`pysam.AlignedSegment`
            An aligned read

        Returns
        -------
        fragment: genomic_interval
            Fragment that corresponds to `read`. Start and end coordinates
            correspond to the 5' ends of the paired reads.
        """
        tlen = read.template_length
        if read.is_reverse:
            end = int(read.reference_end)+self.offset[1]
            start = end + tlen
        else:
            start = int(read.reference_start)+self.offset[0]
            end = start + tlen

        return genome_tools.genomic_interval(read.reference_name, start, end)   

    def lookup(self, interval):
        """Lookup reads in a defined genomic region

        Parameters
        ----------
        intervals: iterable (genomic_interval)

        Returns
        -------
        counts: dict
            Dictionary of read counts (keys: '+' or '-), which contain arrays
            with counts on each strand
        """

        chrom = interval.chrom
        start = interval.start
        end = interval.end
        flip = True if interval.strand == '-' else False

        tmp_fw = {}
        tmp_rev = {}
        reads = []

        for read1, read2 in self.read_pair_generator(chrom, max(start-10, 0), end+10):
            if read1:
                self._add_read(read1, tmp_fw, tmp_rev)
                reads.append(self._get_fragment(read1))
            if read2:
                self._add_read(read2, tmp_fw, tmp_rev)
            
        fw_cutarray = np.array([tmp_fw.get(i, 0.0) for i in range(start, end)])
        rev_cutarray = np.array([tmp_rev.get(i, 0.0) for i in range(start, end)])

        return {
            "+": rev_cutarray[::-1] if flip else fw_cutarray, 
            "-": fw_cutarray[::-1] if flip else rev_cutarray,
            "fragments": reads
        }

    def _validate_genotype(self, read, pos, ref, alt):
        """Validate read genotype

        Parameters
        ----------
        read : :class:`pysam.AlignedSegment`
            Read to validate
        pos: int
            Genomic postion of variant (0-based)
        ref : str
            Reference allele
        alt : str
            Alternate allele

        Returns
        -------
        out : str or bool 

        Raises
        ------
        ReadError
        """

        if not read:
            return None

        try:

            offset = pos-read.reference_start
            if offset<0:
                raise IndexError

            base_call = read.query_sequence[offset]
            qual = read.query_qualities[offset]

            if qual<20:
                raise ReadError(ReadError.ERROR_BASEQ)

            offset_5p = read.reference_end-1-pos if read.is_reverse else pos-read.reference_start
            if offset_5p<=3:
                raise ReadError(ReadError.ERROR_5PROXIMITY)

            mm = int(read.get_tag("XM", with_value_type = False))

            if base_call==ref: 
                if mm>1:
                    raise ReadError(ReadError.ERROR_MISMATCH)
                return ref
            elif base_call==alt:
                if mm>2:
                    raise ReadError(ReadError.ERROR_MISMATCH)
                return alt
            else:
                raise ReadError(ReadError.ERROR_GENOTYPE)
        
        # throws when variant pos is outside of read, handle these differently
        except IndexError as e:
            return None

    def lookup_allelic(self, chrom, start, end, pos, ref, alt, flip=False):
        """Lookup function for allelically resolved counts

        Parameters
        ----------

        Returns
        -------

        """
        #visited_read_pairs = set()

        tmp_ref_fw = {}
        tmp_ref_rev = {}
        tmp_alt_fw = {}
        tmp_alt_rev = {}
        tmp_non_fw = {}
        tmp_non_rev = {}

        reads_ref = []
        reads_alt = []
        reads_non = []
        
        for read1, read2 in self.read_pair_generator(chrom, max(start-10, 0), end+10):

            #try read1 for genotype
            try:
                # Check reads for proper genotypes; raise Genotype error
                # if some thing is problematic, if a read pair doesn't exist or
                # doesn't overlap the SNV then returns 0
                read1_gt = self._validate_genotype(read1, pos, ref, alt)
                read2_gt = self._validate_genotype(read2, pos, ref, alt)

                if (not read1_gt) and (not read2_gt): # neither read overlaps SNV
                    fw = tmp_non_fw
                    rev = tmp_non_rev
                    reads = reads_non
                elif (read1_gt and read2_gt) and read1_gt != read2_gt: # Both reads have a GT, but discordant genotypes
                    raise ReadError(ReadError.ERROR_GENOTYPE)
                elif read1_gt == ref or read2_gt == ref: # Matches REF allele
                    fw = tmp_ref_fw
                    rev = tmp_ref_rev
                    reads = reads_ref
                elif read1_gt == alt or read2_gt == alt: # Matches ALT allele
                    fw = tmp_alt_fw
                    rev = tmp_alt_rev
                    reads = reads_alt
                else:
                    raise ReadError()

                if read1:
                    self._add_read(read1, fw, rev)
                if read2:
                    self._add_read(read2, fw, rev)					

                if read1:
                    reads.append(self._get_fragment(read1))
                else:
                    reads.append(self._get_fragment(read2))
    
            except ReadError as e:
                continue
    
        ref_fw_cutarray = np.array([tmp_ref_fw.get(i, 0.0) for i in range(start, end)])
        ref_rev_cutarray = np.array([tmp_ref_rev.get(i, 0.0) for i in range(start, end)])

        alt_fw_cutarray = np.array([tmp_alt_fw.get(i, 0.0) for i in range(start, end)])
        alt_rev_cutarray = np.array([tmp_alt_rev.get(i, 0.0) for i in range(start, end)])
        
        non_fw_cutarray = np.array([tmp_non_fw.get(i, 0.0) for i in range(start, end)])
        non_rev_cutarray = np.array([tmp_non_rev.get(i, 0.0) for i in range(start, end)])

        return {
            ref: {
                "+": ref_rev_cutarray[::-1] if flip else ref_fw_cutarray, 
                "-": ref_fw_cutarray[::-1] if flip else ref_rev_cutarray,
                "fragments": reads_ref
                },
            alt: {
                "+": alt_rev_cutarray[::-1] if flip else alt_fw_cutarray, 
                "-": alt_fw_cutarray[::-1] if flip else alt_rev_cutarray,
                "fragments": reads_alt
            },
            "other": {
                "+": non_rev_cutarray[::-1] if flip else non_fw_cutarray, 
                "-": non_fw_cutarray[::-1] if flip else  non_rev_cutarray,
                "fragments": reads_non
            }
        }


    def __getitem__(self, x):
        """General function to retrieve cutcounts. Currently supports
        only genomic_intervals and pysam.VariantRecord as input.
        
        Parameters
        ----------
        x : :class:`genome_tools.genomic_interval` or :class:`pysam.VariantRecord`
            Retrieve cleavages over a windowed region or resolve allelically 
            over a known variant
        
        Returns
        -------
        counts : dict
            A dictionary of read counts resolved to each strand, and 
            each allele (if approriate)
        
        Raises
        ------
        TypeError
            Input is neither a :class:`genome_tools.genomic_interval` nor 
            :class:`pysam.VariantRecord`
        """
        if isinstance(x, genome_tools.genomic_interval):
            return self.lookup(x)
        elif isinstance(x, pysam.VariantRecord):
            return self.lookup_allelic(x)
        else:
            raise TypeError(f"Query type not supported: {type(x)}")
