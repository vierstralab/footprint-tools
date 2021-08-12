import sys

from argh.decorators import named, arg

import numpy as np
import scipy.stats

import logging
logger = logging.getLogger(__name__)

from yaspin import yaspin
from yaspin.spinners import Spinners

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

    total_lines = total_passed = 0
    obs_over_exp = []

    logger.info(f"Filtering nucleotides with expected count <= {exp_cutoff} and within a FDR {fdr_cutoff} footprint")

    with yaspin(Spinners.bouncingBar, text='Reading nucleotides') as sp:
        
        filehandle = open(bedgraph_file, 'r')

        for line in filehandle:
            total_lines +=1
            if total_lines % 500000 == 0:
                sp.text = "Reading nucleotides -- {:,} / {:,} passed filters".format(total_passed, total_lines)

            fields = line.strip().split('\t')
            exp = float(fields[3])
            obs = float(fields[4])
            fdr = float(fields[7])

            if fdr <= fdr_cutoff and exp >= exp_cutoff:
                obs_over_exp.append( (obs+1)/(exp+1) )
                total_passed += 1
        
        filehandle.close()
        
        sp.text = "Fitting Beta distribution"

        obs_over_exp = np.array(obs_over_exp)

        i = np.isfinite(obs_over_exp) & (obs_over_exp>0) & (obs_over_exp<1)
        
        prior = scipy.stats.beta.fit(obs_over_exp[i], floc = 0, fscale = 1)[0:2]

    logger.info(f"{total_passed,:} positions used for fitting distribution")
    logger.info(f"Beta distribution parameters - \u03b1 = {prior[0]}, \u03b2 = {prior[1]} ")

    print("%0.4f\t%0.4f" % (prior[0], prior[1]), file = sys.stdout)

    return 0
