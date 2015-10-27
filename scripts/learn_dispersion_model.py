import sys
import argparse

# import footprint_tools
sys.path.append('/home/jvierstra/proj/code/footprint-tools')
from footprint_tools import bed, genomic_interval
from footprint_tools.modeling import bias, dispersion

# numpy
import numpy as np

# fasta index
import pyfaidx

class kmer_action(argparse.Action):
	def __call__(self, parser, namespace, values, option_string = None):
		try:
			setattr(namespace, self.dest, bias.kmer_model(values[0]))
		except IOError, e:
			 raise argparse.ArgumentError(self, str(e))

parser = argparse.ArgumentParser(description = "Learn a dispersion model")
parser.add_argument("bam_file", metavar = "bam_file", type = str, help = "File path to BAM-format tag sequence file")
parser.add_argument("fasta_file", metavar = "fasta_file", type = str, help = "File path to genome FASTA file (requires associated index in same folder)")
parser.add_argument("interval_file", metavar = "interval_file", type = str, help = "File path to BED file")

grp = parser.add_mutually_exclusive_group(required = True)
grp.add_argument("--kmer", metavar = "MODEL_FILE", dest = "bias_model", nargs = 1, action = kmer_action, help = "Use a k-mer model for local bias")
grp.add_argument("--uniform", dest = "bias_model", action = "store_const", const = bias.uniform_model(), help = "Use a uniform model for local bias")

parser.add_argument("--half_win_width", metavar = "N", type = int, help = "Half window width (nt) to apply bias model (default: %(default)s)", default = 10)

parser.add_argument("--procs", metavar = "N", type = int, help = "Number of available processors (default: %(default)s)", default = 8)
parser.add_argument("--chunksize", metavar = "N", type = int, help ="Number of intervals to process per chunk (default: %(default)s", default = 500)

args = parser.parse_args()

# read the intervals BED file
intervals = genomic_interval.genomic_interval_set(bed.bed3_iterator(open(args.interval_file)))

# make expected vs. observed histogram
h = dispersion.build_histogram_multicore(args.bam_file, args.fasta_file, intervals, args.bias_model, half_window_width = args.half_win_width, size = (200, 1000), processes = args.procs, chunksize = args.chunksize)

# learn a dispersion model
dm = dispersion.learn_dispersion_model(h)

# write out model in JSON format
print dispersion.write_dispersion_model(dm)