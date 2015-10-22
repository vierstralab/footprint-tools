import sys
import argparse

#import footprint_tools
sys.path.append('/home/jvierstra/proj/code/footprint-tools')
from footprint_tools import genomic_interval, cutcounts, modeling
from footprint_tools.modeling import bias

#numpy
import numpy as np

#fasta index
import pyfaidx

parser = argparse.ArgumentParser(description = "Compute expected DNase I per-nucleotide cleavages")
parser.add_argument("bam_file", metavar = "bam_file", type = str, help = "File path to BAM-format tag sequence file")
parser.add_argument("fasta_file", metavar = "fasta_file", type = str, help = "File path to genome FASTA file (requires associated index in same folder)")
parser.add_argument("bias_model_file", metavar = "bias_model_file", type = str, help = "File path to a sequence bias model")
parser.add_argument("interval_file", metavar = "interval_file", type = str, help = "File path to BED file")

parser.add_argument("--half_win_width", metavar = "N", type = int, help = "Half window width (nt) to apply bias model (default: %(default)s)", default = 10)

grp = parser.add_mutually_exclusive_group()
grp.add_argument("--smooth_tmean", dest = "smoothing_class", action = "store_const", const = modeling.smoothing.moving_trimmed_mean(), help = "Trimmed mean smoothing of expected window counts", default = None)

args = parser.parse_args()

# open tags file
reads = cutcounts.bamfile(args.bam_file)

# open fasta file
faidx = pyfaidx.Fasta(args.fa_file)

# read the bias model
bm = bias.kmer_model(args.bias_model_file)

# read the intervals BED file
intervals = genomic_interval.genomic_interval_set(intervals_file)

#def windowed_chi_squared(x, w = 3):
#    ret = np.zeros(len(x))
#    for i in np.arange(w, len(x)-w+1):
#        s = sorted(x[i-w:i+w+1])
#    	ret[i] = -2 * np.sum(s[1:-2])
#    return ret

#lnpvals = np.log( [ dm.p_value(e, o) for e, o in zip(exp, obs) ] )
#chisqvals = windowed_chi_squared(lnpvals)

for interval in intervals:

    res = res = modeling.predict_interval(reads, faidx, interval, bm, half_window_width = args.half_win_width, smoothing_class = modeling.smoothing.moving_trimmed_mean())
    
    obs = res["obs"]['+'][1:] + res["obs"]['-'][:-1]
    exp = res["exp"]['+'][1:] + res["exp"]['-'][:-1]
    
    for i in range(len(obs)):
	   sys.stdout.write("%s\t%d\t%d\t%d\t%d\n" % (interval.chrom, interval.start + i + 1, interval.start + i + 2, obs[i], exp[i]))