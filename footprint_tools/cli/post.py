from footprint_tools.data import processor
import click

from multiprocessing import cpu_count

from pandas.core.indexes import interval

import numpy as np
import scipy as sp
import pandas as pd
import pysam

from genome_tools import genomic_interval
from footprint_tools.modeling import dispersion
from footprint_tools.stats import posterior
from footprint_tools.cli.utils import tuple_args, get_kwargs
from footprint_tools.data.process import process
from footprint_tools.data.utils import numpy_collate_concat

from tqdm import tqdm
from tqdm.contrib.logging import logging_redirect_tqdm

import logging
logger = logging.getLogger(__name__)

# kill numpy warnings
np.seterr(all="ignore")

class posterior_stats(process):
    def __init__(self, interval_file, samples_data, fdr_cutoff):

        self.intervals = pd.read_table(interval_file)
        self.samples_data = samples_data
        self.fdr_cutoff = fdr_cutoff

        self.tabix_files = None # these get loaded on first call of __getitem__
        self.disp_models = [dispersion.load_dispersion_model(fn) for fn in self.samples_data["dm_file"]]
        self.betas = self.samples_data[["beta_a", "beta_b"]].to_numpy()

    def _open_tabix_files(self):
        self.tabix_files = [pysam.TabixFile(fn) for fn in self.samples_data["tabix_file"]]

    def _load_data(self, interval):
        assert self.tabix_files or len(self.tabix_files) > 0
        
        n = len(self.tabix_files)
        m = len(interval)

        obs = exp = fdr = w = np.zeros((n, m), dtype=np.float)
        i = j = 0

        for i, tbf in enumerate(self.tabix_files):
            try:
                for row in tbf.fetch(interval.chrom, interval.start,interval.end, parser=pysam.asTuple()):
                    j = int(row[1])-interval.start
                    obs[i, j] = np.float(row[3])
                    exp[i, j] = np.float(row[4])
                    fdr[i, j] = np.float(row[7])
                    w[i, j] = 1
            except:
                pass
        
        return obs, exp, fdr, w

    def cleanup(self):
        """Clean-up filehandlers"""
        if self.tabix_files and len(self.tabix_files) > 0:
            [tbf.close() for tbf in self.tabix_files]
        self.tabix_files = None

    def __len__(self):
        return len(self.intervals)
    
    def __getitem__(self, index):
        """Process data for a single interval"""

        if len(self.tabix_files)==0:
            self._open_tabix_files()

        chrom, start, end = (self.intervals.iat[index, 0], 
                        self.intervals.iat[index, 1], 
                        self.intervals.iat[index, 2])

        interval = genomic_interval(chrom, start, end)
        
        obs, exp, fdr, w = self._load_data(interval)

        prior = posterior.compute_prior_weighted(fdr, w, cutoff=self.fdr_cutoff)   
        delta = posterior.compute_delta_prior(obs, exp, fdr, self.beta, cutoff=self.fdr_cutoff)
        ll_on = posterior.log_likelihood(obs, exp, self.dms, delta=delta, w=3) 
        ll_off = posterior.log_likelihood(obs, exp, self.disp_models, w=3)

        post = -posterior.posterior(prior, ll_on, ll_off)
        post[post<=0] = 0.0

        return {
            'interval': interval,
            'post': post,
        }

@click.command(name='posterior')
@click.argument('sample_data_file')
@click.argument('interval_file')
@click.option('--fdr_cutoff', type=click.FLOAT,
    default=0.05, show_default=True,
    help='FDR cutoff to use when computing priors')
@click.option('--post_cutoff', type=click.FLOAT,
    default=0.2, show_default=True,
    help='Print only positions where that maximum posterior across'
        'all samples meets this threshold. Used to control for'
        'file size.')
@click.option('--n_threads', type=click.IntRange(1, cpu_count()),
    default=cpu_count(), show_default=True,
    help='Number of processors to use')
@click.option('--batch_size', type=click.INT,
    default=100, show_default=True,
    help='Batch size of intervals to process')
@click.option('--outprefix', type=click.STRING,
    default='out', show_default=True,
    help='Output prefix')
def run(sample_data_file,
        interval_file, 
        fdr_cutoff=0.05,
        post_cutoff=0.2,
        n_threads=cpu_count()):
    """Compute footprint posterior probabilities

    \b
    Inputs:
    interval_file       Path to BED-formatted file contained genomic regions to be analyzed
    sample_data_file    Path to file that contains sample data. File is tab-delimited with the
                        following columns:

    \b
                        id          Sample identifier (unique)
                        tabix_file  Path to TABIX-format cleavage statistics file
                        dm_file     Path to dataset JSON-encoded dispersion model file
                        beta_a      \u03b1 parameter (see 'learn_beta' command)
                        beta_b      \u03b2 parameter

                        Note: File must contain a header row. Lines ignored when '#' is first
                        character.

    \b
    Output:

    """

    logger.info(f"Loading sample file {sample_data_file}")

    sample_data = pd.read_table(sample_data_file, header=0)

    logger.info("Verifying input files")
    
    raise click.UsageError("Command not yet implemented")
    