#!/usr/bin/env python

from __future__ import print_function, division

import sys
from argparse import ArgumentParser

import numpy as np

import pysam
import pyfaidx

from genome_tools import bed, genomic_interval, genomic_interval_set

def parse_options(args):

	parser = ArgumentParser(description = "Compute an aggregated DNase I cleavage profile from a collection of regions")

	parser.add_argument("tabix_file", metavar = "tabix_file", type = str,
                        help = "Path to TABIX-format footpringing file")

	parser.add_argument("interval_file", metavar = "interval_file", type = str, 
                        help = "File path to BED file")

	parser.add_argument("output_dir", metavar = "output_dir", type = str, 
                        help = "Path to output directory (must exist)")

	parser.add_argument("--padding", metavar = "N", type = int,
						dest = "padding", default = 25,
						help = "Padding regions by this amount of nuleotides"
						" (default: %(default)s)")

	parser.add_argument("--fasta", metavar = "PATH", type = str,
						dest = "fasta_file",
						help = "Return FASTA sequences")

	parser.add_argument("--phylop", metavar = "PATH", type = str,
						dest = "phylop_file",
						help = "Return phyloP scores")

	parser.add_argument("--posterior", metavar = "PATH", type = str,
					dest = "posterior_file",
					help = "Return phyloP scores")

	parser.add_argument("--posterior-col", metavar = "PATH", type = int,
					dest = "posterior_col", default = 1,
					help = "Return phyloP scores")


	return parser.parse_args(args)


def main(argv = sys.argv[1:]):

	args = parse_options(argv)

	intervals = genomic_interval_set(bed.bed6_iterator(open(args.interval_file)))
	tabix = pysam.TabixFile(args.tabix_file)

	fasta = pyfaidx.Fasta(args.fasta_file) if args.fasta_file else None
	phylop_file = pysam.TabixFile(args.phylop_file) if args.phylop_file else None

	posterior_file = pysam.TabixFile(args.posterior_file) if args.posterior_file else None

	L = len(intervals[0]) + (2*args.padding) + 1
	N = len(intervals)

	exp = np.zeros((N, L), dtype = np.int)
	obs = np.zeros((N, L), dtype = np.int)
	e = np.zeros(L)
	o = np.zeros(L)

	seq = np.chararray((N, L-1))

	posterior = np.zeros((N, L), dtype = np.float)
	posterior_mu = np.zeros(N, dtype = np.float)

	phylop = np.zeros((N, L-1), dtype = np.float)
	p = np.zeros(L)
	
	
	for i, interval in enumerate(intervals):

		start = interval.start - args.padding
		end = interval.end + args.padding + 1

		e.fill(0)
		o.fill(0)

		for row in tabix.fetch(interval.chrom, start, end, parser = pysam.asTuple()):
			j = np.int(row[1]) - start
			e[j] = np.float(row[3])
			o[j] = np.float(row[4])

		if interval.strand == "+":
			exp[i,:] = e
			obs[i,:] = o
		else:
			exp[i,:] = e[::-1]
			obs[i,:] = o[::-1]

		if fasta:
			if interval.strand == "+":
				seq[i,:] = list(fasta[interval.chrom][start:end].seq.upper())[:-1]
			else:
				seq[i,:] = list(fasta[interval.chrom][start:end].complement.seq.upper())[::-1][1:]

		if phylop_file:
			p.fill(0)
			for row in phylop_file.fetch(interval.chrom, start, end, parser = pysam.asTuple()):
				j = np.int(row[1]) - start
				p[j] = np.float(row[4])

			if interval.strand == "+":
				phylop[i,:] = p[:-1]
			else:
				phylop[i,:] = p[::-1][1:]
		
		if posterior_file:
			p.fill(0)
			for row in posterior_file.fetch(interval.chrom, start, end, parser = pysam.asTuple()):
				j = np.int(row[1]) - start
				p[j] = np.float(row[args.posterior_col-1])

			if interval.strand == "+":
				posterior[i,:] = p
			else:
				posterior[i,:] = p[::-1]

			posterior_mu[i] = np.mean(posterior[i,args.padding:-args.padding])

	np.savetxt("%s/exp" % args.output_dir, exp, delimiter = "\t", fmt = "%d")
	np.savetxt("%s/obs" % args.output_dir, obs, delimiter = "\t", fmt = "%d")
	
	if fasta:
		np.savetxt("%s/seq" % args.output_dir, seq, delimiter = "\t", fmt = "%s")
	
	if phylop_file:
		np.savetxt("%s/phylop" % args.output_dir, phylop, delimiter = "\t", fmt = "%0.4f")
	
	if posterior_file:
		np.savetxt("%s/posterior" % args.output_dir, posterior, delimiter = "\t", fmt = "%0.4f")
		np.savetxt("%s/posterior_mu" % args.output_dir, posterior_mu, delimiter = "\t", fmt = "%0.4f")

if __name__ == "__main__":
    sys.exit(main())
		
