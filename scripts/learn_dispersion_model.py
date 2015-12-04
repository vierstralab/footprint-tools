from __future__ import print_function, division

import sys

from argparse import ArgumentParser, Action, ArgumentError

# import footprint_tools
sys.path.append('/home/jvierstra/proj/code/footprint-tools')

from footprint_tools import bed, genomic_interval, cutcounts
from footprint_tools.modeling import bias, predict
from footprint_tools.stats import dispersion

# fasta index
import pyfaidx

#numpy
import numpy as np

#mp
import multiprocessing as mp

#from functools import partial

class kmer_action(Action):
    def __call__(self, parser, namespace, values, option_string = None):
        try:
            setattr(namespace, self.dest, bias.kmer_model(values[0]))
        except IOError, e:
             raise ArgumentError(self, str(e))

def parse_options(args):

    parser = ArgumentParser(description = "Learn a negative binomial dispersion model")

    parser.add_argument("bam_file", metavar = "bam_file", type = str,
                        help = "File path to BAM-format tag sequence file")

    parser.add_argument("fasta_file", metavar = "fasta_file", type = str, 
                        help = "File path to genome FASTA file (requires associated"
                        " FASTA index in same folder; see documentation on how"
                        " to create an index)")

    parser.add_argument("interval_file", metavar = "interval_file", type = str, 
                        help = "File path to BED file")

    grp_bm = parser.add_argument_group("bias modeling options")

    grp_bm.add_argument("--bm", metavar = "MODEL_FILE", dest = "bias_model", 
                        nargs = 1, action = kmer_action, default = bias.uniform_model(),
                        help = "Use a k-mer model for local bias (supplied by file). If"
                        " argument is not provided the model defaults to uniform sequence"
                        " bias.")

    grp_bm.add_argument("--half_win_width", metavar = "N", type = int, default = 5,
                        help = "Half window width to apply bias model."
                        " (default: %(default)s)")

    grp_ot = parser.add_argument_group("other options")

    grp_ot.add_argument("--processors", metavar = "N", type = int,
                        dest = "processors", default = 8,
                        help = "Number of processors to use."
                        " (default: all available processors)")

    return parser.parse_args(args)

def process_func(pred, size):
	
	(obs_counts, exp_counts, win_counts) = pred.compute()

	obs = obs_counts['+'][1:] + obs_counts['-'][:-1]
	exp = exp_counts['+'][1:] + exp_counts['-'][:-1]

	res = np.zeros(size)

	for o, e in zip(obs, exp):
		try:
			res[e, o] += 1.0
		except IndexError:
			pass

	return res

class process_callback(object):

	def __init__(self, size):

		self.res = np.zeros(size)

	def __call__(self, res):

		self.res += res

def main(argv = sys.argv[1:]):

	args = parse_options(argv)

	reads = cutcounts.bamfile(args.bam_file)
	fasta = pyfaidx.Fasta(args.fasta_file)
	intervals = bed.bed3_iterator(open(args.interval_file))

	size = (200, 1000)

	hist_func = process_callback(size)

	pool = mp.Pool(args.processors)

	for interval in genomic_interval.genomic_interval_set(intervals):

		region = predict.prediction(reads, fasta, interval, args.bias_model, args.half_win_width, 0, 0)

		pool.apply_async(process_func, args = (region, size,), callback = hist_func)

		while pool._taskqueue.qsize() > 1000:
			pass

	pool.close()
	pool.join()

	# Learn model from histogram
	model = dispersion.learn_dispersion_model(hist_func.res)

	# Write model
	print(dispersion.write_dispersion_model(model), file = sys.stdout)

	# Success!
	return 0

if __name__ == "__main__":
    sys.exit(main())
