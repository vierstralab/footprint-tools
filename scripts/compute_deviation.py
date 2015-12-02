import sys

from argparse import ArgumentParser, Action, ArgumentError

# import footprint_tools
sys.path.append('/home/jvierstra/proj/code/footprint-tools')

from footprint_tools import bed, genomic_interval, cutcounts, stats
from footprint_tools.modeling import bias, predict
from footprint_tools.stats import fdr, dispersion

# fasta index
import pyfaidx

#numpy
import numpy as np

class kmer_action(Action):
    def __call__(self, parser, namespace, values, option_string = None):
        try:
            setattr(namespace, self.dest, bias.kmer_model(values[0]))
        except IOError, e:
             raise ArgumentError(self, str(e))

class dispersion_model_action(Action):
    def __call__(self, parser, namespace, values, option_string = None):
        try:
            setattr(namespace, self.dest, dispersion.load_dispersion_model(values[0]))
        except IOError, e:
             raise ArgumentError(self, str(e))

# TODO: implement
# class list_parse_action(Action):
#    def __call__(self, parser, namespace, values, option_string = None):
#        try:
#            pass
#        except IOError, e:
#             raise ArgumentError(self, str(e))

def parse_options(args):

    parser = ArgumentParser(description = "Compute expected DNase I per-nucleotide cleavages")

    parser.add_argument("bam_file", metavar = "bam_file", type = str,
                        help = "File path to BAM-format tag sequence file")

    parser.add_argument("fasta_file", metavar = "fasta_file", type = str, 
                        help = "File path to genome FASTA file (requires associated"
                        " FA index in same folder)")

    parser.add_argument("interval_file", metavar = "interval_file", type = str, 
                        help = "File path to BED file")

    grp_bm = parser.add_argument_group("bias modeling options")

    grp_bm.add_argument("--kmer", metavar = "MODEL_FILE", dest = "bias_model", 
                        nargs = 1, action = kmer_action, default = bias.uniform_model(),
                        help = "Use a k-mer model for local bias (supplied by file). If"
                        " argument is not provided the model defaults to uniform sequence"
                        " bias.")

    grp_bm.add_argument("--half_win_width", metavar = "N", type = int, default = 5,
                        help = "Half window width to apply bias model."
                        " (default: %(default)s)")

    grp_sm = parser.add_argument_group("smoothing options")
    
    grp_sm.add_argument("--smooth_half_win_width", metavar = "N", type = int, default = 50,
                        help = "Half window width to apply smoothing model. When set to"
                        " zero no smoothing is applied. (default: %(default)s)")

    grp_sm.add_argument("--smooth_clip", metavar = "N", type = float, default = 0.05,
                        help = "Fraction of signal to clip when computing trimmed mean."
                        " (default: %(default)s)")

    grp_st = parser.add_argument_group("statistics options")

    grp_st.add_argument("--disp_model", nargs = 1, metavar = "MODEL_FILE", 
                        dest = "dispersion_model", action = dispersion_model_action, default = None,
                        help = "Dispersion model for negative binomial tests. If argument"
                        " is not provided then no stastical output is provided. File is in"
                        " JSON format and generated using the 'learn_dispersion_model'"
                        " script included in the software package.")

    # TODO: implement
    # grp_st.add_argument("--fdr_cutoffs", nargs = 1, metavar = "[N, ...]", dest = "fdr_cutoffs",
    #                     action = list_parse_action, default = [],
    #                     help = "FDR cutoff at which to report footprints.")

    grp_st.add_argument("--fdr_shuffle", nargs = 1, metavar = "N", type = int,
                        dest = "fdr_shuffle", default = 25,
                        help = "Number of times to shuffle data for FDR calculation."
                        " (default: %(default)s)")

    return parser.parse_args(args)


def main(argv = sys.argv[1:]):

    args = parse_options(argv)

    reads = cutcounts.bamfile(args.bam_file)
    fasta = pyfaidx.Fasta(args.fasta_file)
    intervals = bed.bed3_iterator(open(args.interval_file))

    bm = args.bias_model
    dm = args.dispersion_model

    win_pvals_func = lambda x: stats.stouffers_z(np.ascontiguousarray(x), 3)

    formatter = fmt = "\t{:0.0f}\t{:0.0f}"
    if args.dispersion_model:
        formatter += "\t{:0.4g}\t{:0.4g}\t{:0.4g}"
    
    for interval in intervals:

        counts = predict.predict_interval(reads, fasta, interval, args.bias_model, args.half_win_width, args.smooth_half_win_width, args.smooth_clip)
        exp = counts["exp"]['+'][1:] + counts["exp"]['-'][:-1]
        obs = counts["obs"]['+'][1:] + counts["obs"]['-'][:-1]

        if args.dispersion_model:

            pvals = dm.p_values(exp, obs)
            win_pvals = win_pvals_func(pvals)

            # Re-shuffle data
            pvals_null = dm.resample_p_values(exp, args.fdr_shuffle)
            win_pvals_null = np.apply_along_axis(win_pvals_func, 0, pvals_null)

            # False-positive rate at each base
            fpr = fdr.emperical_fpr(win_pvals_null, win_pvals)

            # Segment FPR cutoffs
            # for cutoff in args.fdr_cutoffs:
            #    pass

            out = np.column_stack((exp, obs, pvals, win_pvals, fpr))

        else:

            out = np.column_stack((exp, obs))

        for i in range(out.shape[0]):
            print '{}\t{}\t{}'.format(interval.chrom, interval.start+i, interval.start+i+1) + formatter.format(*(out[i,:]))

    return 0
    
if __name__ == "__main__":
    sys.exit(main())
