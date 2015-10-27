import sys 

from footprint_tools import genomic_interval

def bed3_iterator(file_handle):

	for line in file_handle:
		
		fields = line.strip().split()
		
		chrom = fields[0]
		start = int(fields[1])
		end = int(fields[2])

		yield genomic_interval.genomic_interval(chrom, start, end)

def bed6_iterator(file_handle):

	for line in file_handle:
		
		fields = line.strip().split()
		
		chrom = fields[0]
		start = int(fields[1])
		end = int(fields[2])
		name = fields[3]
		strand = fields[4]
		score = fields[5]

		yield genomic_interval.genomic_interval(chrom, start, end, name, strand)