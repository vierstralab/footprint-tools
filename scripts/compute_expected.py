import sys
import argparse

# import footprint_tools
sys.path.append('/home/jvierstra/proj/code/footprint-tools')
from footprint_tools import bed, genomic_interval, cutcounts, modeling
from footprint_tools.modeling import bias, smoothing, dispersion

# numpy
import numpy as np

#scipy
import scipy.stats

# fasta index
import pyfaidx

class kmer_action(argparse.Action):
    def __call__(self, parser, namespace, values, option_string = None):
        try:
            setattr(namespace, self.dest, bias.kmer_model(values[0]))
        except IOError, e:
             raise argparse.ArgumentError(self, str(e))

class dispersion_model_action(argparse.Action):
    def __call__(self, parser, namespace, values, option_string = None):
        try:
            setattr(namespace, self.dest, dispersion.load_dispersion_model(values[0]))
        except IOError, e:
             raise argparse.ArgumentError(self, str(e))

parser = argparse.ArgumentParser(description = "Compute expected DNase I per-nucleotide cleavages")
parser.add_argument("bam_file", metavar = "bam_file", type = str, help = "File path to BAM-format tag sequence file")
parser.add_argument("fasta_file", metavar = "fasta_file", type = str, help = "File path to genome FASTA file (requires associated index in same folder)")
parser.add_argument("interval_file", metavar = "interval_file", type = str, help = "File path to BED file")

grp_bm = parser.add_mutually_exclusive_group(required = True)
grp_bm.add_argument("--kmer", metavar = "MODEL_FILE", dest = "bias_model", nargs = 1, action = kmer_action, help = "Use a k-mer model for local bias")
grp_bm.add_argument("--uniform", dest = "bias_model", action = "store_const", const = bias.uniform_model(), help = "Use a uniform model for local bias")

grp_smooth = parser.add_mutually_exclusive_group()
grp_smooth.add_argument("--smooth_tmean", dest = "smoothing_class", action = "store_const", const = smoothing.moving_trimmed_mean(), help = "Trimmed mean smoothing of expected window counts", default = None)

parser.add_argument("--half_win_width", metavar = "N", type = int, help = "Half window width (nt) to apply bias model (default: %(default)s)", default = 5)

parser.add_argument("--disp_model", nargs = 1, metavar = "MODEL_FILE", dest = "dispersion_model", action = dispersion_model_action, help = "Compute p-values using a custom dispersion model (JSON format)", default = None)

args = parser.parse_args()

# open tags file
reads = cutcounts.bamfile(args.bam_file)

# open fasta file
faidx = pyfaidx.Fasta(args.fasta_file)

# read the intervals BED file
intervals = genomic_interval.genomic_interval_set(bed.bed3_iterator(open(args.interval_file)))

def windowed_chi_squared(x, w = 3):
    """
    Compute chi-squared values from groups of ln p-values
    """

    chi = np.zeros(len(x))
    p = np.zeros(len(x))
    for i in np.arange(w, len(x)-w+1):
        #s = sorted(x[i-w:i+w+1])
    	#ret[i] = -2 * np.sum(s[1:-2])
        chi[i] = -2 * np.sum(x[i-w:i+w+1])
        p[i] = scipy.stats.chi2.logsf(chi[i], 6)
    return (chi, p)

#
for interval in intervals:

    res = modeling.predict_interval(reads, faidx, interval, args.bias_model, half_window_width = args.half_win_width, smoothing_class = args.smoothing_class)
    
    exp = res["exp"]['+'][1:] + res["exp"]['-'][:-1]
    obs = res["obs"]['+'][1:] + res["obs"]['-'][:-1]
    
    if args.dispersion_model:
        """
        If user supplies a negative binomial fit, calculate 
        p-values and smoothed chi-squared scores
        """

        lnpvals_down = np.array([args.dispersion_model.log_p_value(e, o) for e, o in zip(exp, obs)]) 
        chisq_vals, chisq_lnpvals = windowed_chi_squared(lnpvals_down)

        for i in range(len(obs)):
            sys.stdout.write("%s\t%d\t%d\t%d\t%d\t%0.4f\t%0.4f\t%0.4f\n" % (interval.chrom, interval.start + i + 1, interval.start + i + 2, obs[i], exp[i], lnpvals_down[i], chisq_vals[i], chisq_lnpvals[i]))
    else:
        """
        If not, ouput just the observed and expected counts 
        as per bias model and smoothing strategy
        """

        for i in range(len(obs)):
            sys.stdout.write("%s\t%d\t%d\t%d\t%d\n" % (interval.chrom, interval.start + i + 1, interval.start + i + 2, obs[i], exp[i]))
