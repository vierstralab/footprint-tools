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
from functools import partial

class kmer_action(Action):
    def __call__(self, parser, namespace, values, option_string = None):
        try:
            setattr(namespace, self.dest, bias.kmer_model(values[0]))
        except IOError, e:
             raise ArgumentError(self, str(e))

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

    grp_st = parser.add_argument_group("other options")

    grp_st.add_argument("--processors", metavar = "N", type = int,
                        dest = "processors", default = 8,
                        help = "Number of processors to use."
                        " (default: all available processors)")

    return parser.parse_args(args)


def chunks_list(x, chunksize):

	"""Chunk data for parallelization"""

	n = max(1, chunksize)
	return [ x[i:i+n] for i in range(0, len(x), n) ]

def build_histogram(intervals, bam_file, fasta_file, bm, half_window_width, size):
	
	"""Build histogram"""

	reads = cutcounts.bamfile(bam_file)
	fasta = pyfaidx.Fasta(fasta_file)

	h = np.zeros(size)

	for interval in intervals:
		
		res = predict.predict_interval(reads, fasta, interval, bm, half_window_width, 0, 0.0)
		
		obs = res["obs"]['+'][1:] + res["obs"]['-'][:-1]
		exp = res["exp"]['+'][1:] + res["exp"]['-'][:-1]
		
		for o, e in zip(obs, exp):
			try:
				h[e, o] += 1.0
			except IndexError:
				pass

	return h

def main(argv = sys.argv[1:]):

	args = parse_options(argv)

	# Make a partial function to wrap arguments
	wrapper_func = partial(build_histogram, bam_file = args.bam_file, fasta_file = args.fasta_file, bm = args.bias_model, half_window_width = args.half_win_width, size = (200, 1000))

	# Read intervals into memory
	intervals = genomic_interval.genomic_interval_set(bed.bed3_iterator(open(args.interval_file)))

    # Run pool of threads
	pool = mp.Pool(processes = args.processors)
	res_chunks = pool.map(wrapper_func, chunks_list(intervals, len(intervals) / args.processors))

	# Combine the histograms
	res = res_chunks[0]
	for i in np.arange(1, len(res_chunks)):
		res += res_chunks[i]

	# Learn histogram
	dm = dispersion.learn_dispersion_model(res)

	# Save it
	print >> sys.stdout, dispersion.write_dispersion_model(dm)

	# Success!
	return 0

if __name__ == "__main__":
    sys.exit(main())
