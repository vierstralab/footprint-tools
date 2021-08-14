import click
from click_option_group import optgroup

from multiprocessing import cpu_count

import numpy as np
import scipy as sp
import pandas as pd

import pysam
pysam.set_verbosity(0)

from genome_tools import genomic_interval
from genome_tools.data.dataset import dataset

from footprint_tools.modeling import dispersion
from footprint_tools.stats import posterior

from footprint_tools.cli.utils import (verify_tabix_file, write_stats_to_output, 
                                        write_output_header)

from tqdm import tqdm
from tqdm.contrib.logging import logging_redirect_tqdm

import logging
logger = logging.getLogger(__name__)

# kill numpy warnings
np.seterr(all="ignore")

# columns required in sample data file
required_sample_data_cols = ["id", "tabix_file", "dm_file", "beta_a", "beta_b"]

class posterior_stats(dataset):
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
        """TODO: Use a genome_tools.data.extractors.tabix_extractor"""
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
@optgroup.group('Output options')
@optgroup.option('--post_cutoff', type=click.FLOAT,
    default=0.2, show_default=True,
    help='Print only positions where that maximum posterior across '
        'all samples meets this threshold. Used to control for '
        'file size.')
@optgroup.option('--outprefix', type=click.STRING,
    default='out', show_default=True,
    help='Output prefix')
@optgroup.group('Other options')
@optgroup.option('--n_threads', type=click.IntRange(1, cpu_count()),
    default=cpu_count(), show_default=True,
    help='Number of processors to use')
@optgroup.option('--batch_size', type=click.INT,
    default=100, show_default=True,
    help='Batch size of intervals to process')
def run(sample_data_file,
        interval_file, 
        fdr_cutoff=0.05,
        post_cutoff=0.2,
        batch_size=100,
        n_threads=cpu_count(),
        outprefix='out'):
    """Compute footprint posterior probabilities

    Applies an emperical Bayesian approach to compute the posterior probability a
    nucleotide is protected by jointly analyzing many datasets.

    INTERVAL_FILE is a BED-formatted file containing genomic regions to be 
    analyzed. SAMPLE_DATA_FILE file that specifying sample metadata. 
    
    SAMPLE_DATA_FILE is tab-delimited with the following columns::
    
    \b
        id          Sample identifier (unique)
        tabix_file  Path to TABIX-format cleavage statistics file
        dm_file     Path to dataset JSON-encoded dispersion model file
        beta_a      \u03b1 parameter (see 'learn_beta' command)
        beta_b      \u03b2 parameter
    
    Note: File must contain a header row. Lines ignored when '#' 
    is first character.

    A bedGraph-like file ({N_samples}+3 columns) with the folowing
    annotations::

    \b
        contig start start+1 -log(1-p)_1 ... -log(1-p)_N
    
    where N is the total number samples. Note that values are
    1-posterior probability. Columns (samples) are in the same 
    order as the sample data file.
    """

    logger.info(f"Loading sample data file {sample_data_file}")
    try:
        sample_data = pd.read_table(sample_data_file, header=0, comment='#')
        if not all([col in sample_data.columns for col in required_sample_data_cols]):
            raise ValueError(f"Sample data file does not contain the required columns: {required_sample_data_cols}")
    except Exception as e:
        logger.critical(e)
        raise click.Abort()

    # Validate and load inputs
    logger.info("Validating input files")
    try:
        [verify_tabix_file(fn) for fn in sample_data["tabix_file"]]
        [open(fn).close() for fn in sample_data["dm_file"]]

    except IOError as e:
        logger.critical(e)
        raise click.Abort()
    
    # Open output files
    try:
        # Output stats file
        output_bedgraph_file = outprefix + '.bedgraph'
        logger.info(f"Writing per-nucleotide stats to {output_bedgraph_file}")
        output_bedgraph_filehandle = open(output_bedgraph_file , 'w')

        #write header lines
        write_output_header(sample_data["id"], output_bedgraph_filehandle)

    except IOError as e:
        logger.critical(e)
        raise click.Abort()

    # Create data processor and iterator
    dl = posterior_stats(interval_file, sample_data, fdr_cutoff)
    dl_iter = dl.batch_iter(batch_size=batch_size, num_workers=n_threads)

    # filter function to apply when writing posteriors to file
    filter_fn = lambda x: np.min(x, axis=1) > post_cutoff

    with logging_redirect_tqdm():

        for batch in tqdm(dl_iter, colour='#cc951d'):
            
            for interval, stats in zip(batch["interval"], batch["stats"]):
                write_stats_to_output(interval, stats, output_bedgraph_filehandle, filter_fn=filter_fn)
    
    output_bedgraph_filehandle.close()
