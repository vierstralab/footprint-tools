import sys
import argparse

#import footprint_tools
sys.path.append('/home/jvierstra/proj/code/footprint-tools')
from footprint_tools import genomic_interval
from footprint_tools.modeling import bias, dispersion

#numpy
import numpy as np

#fasta index
import pyfaidx

parser = argparse.ArgumentParser(description = "Learn a dispersion model")
parser.add_argument("bam_file", metavar = "bam_file", type = str, help = "File path to BAM-format tag sequence file")
parser.add_argument("fasta_file", metavar = "fasta_file", type = str, help = "File path to genome FASTA file (requires associated index in same folder)")
parser.add_argument("bias_model_file", metavar = "bias_model_file", type = str, help = "File path to a sequence bias model")
parser.add_argument("interval_file", metavar = "interval_file", type = str, help = "File path to BED file")

parser.add_argument("--half_win_width", metavar = "N", type = int, help = "Half window width (nt) to apply bias model (default: %(default)s)", default = 10)
parser.add_argument("--procs", metavar = "N", type = int, help = "Number of available processors (default: %(default)s)", default = 8)

args = parser.parse_args()


# read the bias model
bm = bias.kmer_model(args.bias_model_file)

# read the intervals BED file
intervals = genomic_interval.genomic_interval_set(args.interval_file)

# make expected vs. observed histogram
h = dispersion.build_histogram_multicore(args.bam_file, args.fasta_file, intervals, bm, half_window_width = args.half_win_width, size = (200, 1000), processes = args.procs, chunksize = 500)

# learn a dispersion model
dm = dispersion.learn_dispersion_model(h)

# write out model in JSON format
print dispersion.write_dispersion_model(dm)