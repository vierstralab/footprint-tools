import argh
from argh.decorators import named, arg
from argh.exceptions import CommandError

from tqdm import tqdm
import multiprocessing as mp

import numpy as np
import scipy as sp
import pysam

from genome_tools import bed

import footprint_tools
from footprint_tools import cutcounts
from footprint_tools.modeling import bias, dispersion, predict
from footprint_tools.cli.utils import chunkify

import logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

def tuple_ints(arg):
	"""
	Function to parser a tuple of integers from command line
	"""
	try:
		fw, rev = list(map(int, arg.split(',')))
		return (fw, rev)
	except:
		raise CommandError("Argument must be a tuple -- i.e., 0,-1")

class process_callback():
	"""
	Class that aggregates results from processor functions
	In this case, it just sums up the 2d histogram matricies
	"""
	def __init__(self, hist_size):
		self.x = np.zeros(hist_size)

	def __call__(self, other):
		self.x += other

def process_func(bam_file, fasta_file, bm, intervals, hist_size, proc_id, **kwargs):
	"""
	Function that can be used to parallize reading 
	the BAM files and computing expected counts

	Returns an np.array histogram of dimension hist_size
	"""

	bam_kwargs = { k:kwargs.pop(k) for k in ["min_qual", "remove_dups", "remove_qcfail", "offset"] }
	bam_reader = cutcounts.bamfile(bam_file, **bam_kwargs)
	fasta_reader = pysam.FastaFile(fasta_file)

	predictor = predict.prediction(bam_reader, fasta_reader, bm, **kwargs)

	res_hist = np.zeros(hist_size)

	progress_desc = "Chunk {}".format(proc_id)

	with tqdm(total=len(intervals), desc=progress_desc, position=proc_id+1, ncols=80) as progress_bar:

		for interval in intervals:

			obs, exp, win = predictor.compute(interval)
			
			obs = obs['+'][1:] + obs['-'][:-1]
			exp = exp['+'][1:] + exp['-'][:-1]
		
			for o, e in zip(obs, exp):
				try:
					res_hist[int(e), int(o)] += 1.0
				except IndexError:
					pass

			progress_bar.update(1)

	return res_hist

@named('learn_dm')
@arg('interval_file', 
	help="File path to BED file")
@arg('bam_file', 
	help='Path to BAM-format tag alignment file')
@arg("fasta_file", help = "Path to genome FASTA file (requires associated FASTA index in same folder; see documentation on how to create an index)")
@arg("--bias_model_file", help="Use a k-mer model for local bias (supplied by file). If argument is not provided the model defaults to uniform sequence bias.")
@arg("--min_qual", help="Ignore reads with mapping quality lower than this threshold.", default=1)
@arg("--remove_dups", help="Remove duplicate reads from analysis", default=False)
@arg("--remove_qcfail", help="Remove QC-failed reads from analysis", default=False)
@arg("--offset", help="BAM file offset (enables support for other datatypes -- e.g. Tn5/ATAC)", default=(0,-1), type=tuple_ints)
@arg("--half_win_width", help="Half window width to apply bias model.", default=5)
@arg("--n_threads", help="Number of processors to use. (default: all available processors)", default=mp.cpu_count())
def run(interval_file,
		bam_file,
		fasta_file,
		bias_model_file=None,
		min_qual=1,
		remove_dups=False,
		remove_qcfail=False,
		bam_offset=(0, -1),
		half_win_width=5,
		n_threads=mp.cpu_count()):
    """
    Learn a negative binomial dispersion model from data corrected for intrinsic sequence preference.
    """
	hist_size = (200, 1000) # hard coded histogram size -- for now...
	hist_agg = callback(hist_size)

	intervals = list(bed.bed3_iterator(open(interval_file)))

	# Load bias model (if specified)
	if bias_model_file:
		logger.info(f"Loading bias model from file {bias_model_file}")
		bm = bias.kmer_model(bias_model_file)
	else:
		logger.info(f"No bias model file specified -- using uniform model")
		bm = bias.uniform_model()

	logger.info(f"BED file contains {len(intervals):,} regions")

	proc_kwargs = {
		"min_qual": min_qual,
		"remove_dups": remove_dups,
		"remove_qcfail": remove_qcfail,
		"offset": bam_offset,
		"half_window_width": half_win_width,
	}

	pool = mp.Pool(n_threads)

	logger.info(f"Using {n_threads} threads to computed expected cleavage counts")

	for i, chunk in enumerate(chunkify(intervals, n_threads)):
		pool.apply_async(process_func, args=(bam_file, fasta_file, bm, chunk, hist_size, i), kwds=proc_kwargs, callback=hist_agg)

	pool.close()
	pool.join()

	logger.info("Finished computing expected cleavage counts!")

	# Learn model from histogram

	logger.info("Learning dispersion model")

	model = dispersion.learn_dispersion_model(hist_agg.x)

	# Write model

	model_file = os.path.abspath(os.path.join(os.getcwd(), "dm.json"))

	logger.info("Writing dispersion model to {}".format(model_file))

	with open(model_file, "w") as f:
		print(dispersion.write_dispersion_model(model), file = f)

	# Success!
	return 0