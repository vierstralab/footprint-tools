from os import cpu_count
import sys

import argh
from argh.decorators import named, arg

from tqdm import tqdm
import multiprocessing as mp
from functools import partial

import pysam

import numpy as np
import scipy.stats

import logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

def process_func(contig, tabix_filehandle, fdr_cutoff, exp_cutoff):
	"""
	Fetch nucleotides per contig
	"""
	res = []

	for row in tabix_filehandle.fetch(contig=contig, parser=pysam.asTuple(), multiple_iterators=True):

		exp = np.float64(row[3])
		obs = np.float64(row[4])
		fdr = np.float64(row[7])

		if fdr <= fdr_cutoff and exp >= exp_cutoff:
			res.append( (obs+1)/(exp+1) )
	return res

@named('learn_beta')
@arg('tabix_file',
	help='Path to TABIX-format file (ouput from "find_fps" command')
@arg('--fdr_cutoff',
	type=float,
	default=0.05,
	help='Only consider nucleotides with FDR <= this value')
@arg('--exp_cutoff',
	type=int,
	default=10,
	help='Only consider nucleotides with expected cleavages >= this value')
@arg('--n_threads',
	type=int,
	help='Number of processors to use (min=1)',
	default=max(1, mp.cpu_count()))
def run(tabix_file, fdr_cutoff=0.05, exp_cutoff=10, n_threads=max(1, mp.cpu_count())):
	"""Learn the parameters of a Beta distribution for a single sample.

	Output:

	Note: This step is required to compute posterior footprint probabilities.
	"""

	tabix_filehandle = pysam.TabixFile(tabix_file)
	contigs = tabix_filehandle.contigs
	ratios = []

	pool = mp.Pool(n_threads)
	func = partial(process_func, tabix_filehandle=tabix_filehandle, fdr_cutoff=fdr_cutoff, exp_cutoff=exp_cutoff)

	progress_desc="Contigs processed"
	for item in tqdm(pool.imap_unordered(func, contigs), desc=progress_desc, total=len(contigs), ncols=80):
		ratios.extend(item)

	ratios = np.array(ratios)
	logger.info(f"{len(ratios):,} positions analyzed")

	i = np.isfinite(ratios) & (ratios>0) & (ratios<1)
	logger.info(f"Using {np.sum(i)} positions to to copmute Beta parameters")

	(a, b) = scipy.stats.beta.fit(ratios[np.isfinite(ratios) & (ratios>0) & (ratios<1)], floc = 0, fscale = 1)[0:2]

	print("f{a}\t{b}", file=sys.stdout)
