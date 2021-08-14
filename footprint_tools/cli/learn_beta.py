import click
from click_option_group import optgroup

from yaspin import yaspin
from yaspin.spinners import Spinners

import numpy as np
import scipy.stats

import logging
logger = logging.getLogger(__name__)

@click.command(name='learn_beta')
@click.argument('bedgraph_file')
@click.option('--fdr_cutoff',
    type=click.FLOAT, default=0.05,
    help='Only consider nucleotides with FDR <= this value')
@click.option('--exp_cutoff',
    type=click.INT, default=10,
    help='Only consider nucleotides with expected cleavages >= this value')
@optgroup.group('Output options')
@optgroup.option('--outfile',
    type=click.STRING, default='beta.txt',
    help='Output file path')
def run(bedgraph_file, 
        fdr_cutoff=0.05, 
        exp_cutoff=10,
        outfile='beta.txt'):
    """Learn the parameters of a Beta distribution

    BEDGRAPH_FILE is a path to bedgraph file containg the per-nucleotide 
    statistics generated from the 'detect' command.
    
    Outputs the fit Beta-distribution parameters to <outfile>

    Note: This step is required to compute posterior footprint probabilities.
    """

    total_lines = total_passed = 0
    obs_over_exp = []

    logger.info(f"Filtering nucleotides with expected count >= {exp_cutoff} and within a FDR {fdr_cutoff} footprint")

    with yaspin(Spinners.bouncingBar, text='Reading nucleotides') as sp:
        
        try:
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
            
        except IOError as e:
            logger.critical(e)
            click.Abort()
        
        sp.text = "Fitting Beta distribution"

        obs_over_exp = np.array(obs_over_exp)

        i = np.isfinite(obs_over_exp) & (obs_over_exp>0) & (obs_over_exp<1)
        
        prior = scipy.stats.beta.fit(obs_over_exp[i], floc = 0, fscale = 1)[0:2]

    logger.info(f"{total_passed:,} positions used for fitting distribution")
    logger.info(f"Beta distribution parameters: \u03b1 = {prior[0]}, \u03b2 = {prior[1]}")

    logger.info(f"Writing parameters to file {outfile}")

    try:
        with open(outfile, 'w') as f:
            print("%0.4f\t%0.4f" % (prior[0], prior[1]), file = f)
    except IOError as e:
        logger.critical(e)
        click.Abort()

