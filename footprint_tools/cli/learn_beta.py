import sys

import argh
from argh.decorators import named, arg

import pysam

import numpy as np
import scipy.stats

import logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

@named('learn_beta')
@arg('tabix_file',
    help='Path to TABIX-format file (e.g., ouput from "find_fps" command')
@arg('--fdr_cutoff',
    type=float,
    default=0.05,
    help='Only consider nucleotides with FDR <= this value.')
@arg('--exp_cutoff',
    type=int,
    default=10,
    help='Only consider nucleotides with expected cleavages >= this value.')
def run(tabix_file, fdr_cutoff=0.05, exp_cutoff=10):
    """Learn the parameters of a Beta distribution for a single sample.

    This step is required to compute posterior footprint probabilities.
    """
    handle = pysam.TabixFile(tabix_file)

	ratios = []

	for row in handle.fetch(parser = pysam.asTuple()):

		exp = np.float64(row[3])
		obs = np.float64(row[4])
		fdr = np.float64(row[7])

		if fdr <= fdr_cutoff and exp >= exp_cutoff:
			ratios.append( (obs+1)/(exp+1) )

	handle.close()

	ratios = np.array(ratios)
	(a, b) = scipy.stats.beta.fit(ratios[np.isfinite(ratios) & (ratios>0) & (ratios<1)], floc = 0, fscale = 1)[0:2]

	print("f{a}\t{b}", file=sys.stdout)

    return 0