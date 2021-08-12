import sys

from argh.decorators import named, arg

import pysam

import numpy as np
import scipy.stats

import logging
logger = logging.getLogger(__name__)

from yaspin import yaspin
from yaspin.spinners import Spinners

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
@arg('bedgraph_file',
    help='Path to bedgraph file (ouput from "detect" command')
@arg('--fdr_cutoff',
    type=float,
    help='Only consider nucleotides with FDR <= this value')
@arg('--exp_cutoff',
    type=int,
    help='Only consider nucleotides with expected cleavages >= this value')
def run(bedgraph_file, 
        fdr_cutoff=0.05, 
        exp_cutoff=10):
    """Learn the parameters of a Beta distribution for a single sample.

    Output:

    Note: This step is required to compute posterior footprint probabilities.
    """
    
    filehandle = open(bedgraph_file, 'r')

    total_lines = 0
    with yaspin(Spinners.hamburger) as sp:
        for line in filehandle:
            total_lines +=1

            if total_lines % 1000000:
                sp.text = total_lines

    return 0
