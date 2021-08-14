import click
from click_option_group import optgroup

from multiprocessing import cpu_count

import numpy as np
import scipy as sp
import pandas as pd
import pysam

from genome_tools import genomic_interval
from genome_tools.data.dataset import dataset
from genome_tools.data.utils import numpy_collate_concat

from footprint_tools import cutcounts
from footprint_tools.modeling import bias, predict, dispersion

from footprint_tools.cli.utils import (tuple_args, get_kwargs, verify_bam_file, verify_fasta_file)

from tqdm import tqdm
from tqdm.contrib.logging import logging_redirect_tqdm

import logging
logger = logging.getLogger(__name__)

# kill numpy warnings
np.seterr(all="ignore")

class expected_counts(dataset):
    """Computes observed and expected cleavage counts
    """
    def __init__(self, interval_file, bam_file, fasta_file, bm, **kwargs):
       
        self.intervals = pd.read_table(interval_file, header=None)

        logger.info(f"BED file contains {len(self.intervals):,} regions")

        self.bam_file = bam_file
        self.fasta_file = fasta_file
        self.bm = bm

        # kwargs for cutcounts.bamfile
        self.counts_reader_kwargs = get_kwargs(
            [
                'min_qual',
                'remove_dups',
                'remove_qcfail',
                'offset'
            ], 
            kwargs)

        self.fasta_reader_kwargs = {}

        self.counts_predictor_kwargs = get_kwargs(
            [
                'half_win_width',
            ],
            kwargs)

        self.counts_extractor = None
        self.fasta_extractor = None
        self.count_predictor = None

    def __len__(self):
        return len(self.intervals)

    def __getitem__(self, index):
        """Process data for a single interval"""
        
        # Open file handlers on first call. This avoids problems when
        # parallel processing data with non-thread safe code (i.e., pysam)
        if not self.counts_extractor:
            self.counts_extractor = cutcounts.bamfile(self.bam_file, **self.counts_reader_kwargs)
            self.fasta_extractor = pysam.FastaFile(self.fasta_file, **self.fasta_reader_kwargs)
            self.count_predictor = predict.prediction(self.counts_extractor, self.fasta_extractor, 
                                                        self.bm, **self.counts_predictor_kwargs)

        chrom, start, end = (self.intervals.iat[index, 0], 
                             self.intervals.iat[index, 1], 
                             self.intervals.iat[index, 2])

        interval = genomic_interval(chrom, start, end)

        obs, exp, _ = self.count_predictor.compute(interval)
        obs = obs['+'][1:] + obs['-'][:-1]
        exp = exp['+'][1:] + exp['-'][:-1]

        return np.column_stack((exp, obs))

@click.command(name='learn_dm')
@click.argument('interval_file')
@click.argument('bam_file')
@click.argument('fasta_file')
@optgroup.group('Bias modeling options')
@optgroup.option('--bias_model_file', type=click.STRING,
    help='Use a k-mer model for sequence bias (supplied by file). '
        'If argument is not provided the model defaults to uniform '
        'sequence bias.')
@optgroup.option('--half_win_width', type=click.INT,
    default=5, show_default=True,
    help='Half window width to apply bias model')
@optgroup.group('Read filtering options')
@optgroup.option('--min_qual', type=click.INT, 
    default=1, show_default=True,
    help='Ignore reads with mapping quality lower than this threshold')
@optgroup.option('--keep_dups', type=click.BOOL,
    default=False, show_default=True,
    help='Keep duplicate reads')
@optgroup.option('--keep_qcfail', type=click.BOOL,
    default=False, show_default=True,
    help='Keep QC-failed reads')
@optgroup.group('Output options')
@optgroup.option('--outfile', type=click.STRING,
    default='dm.json', show_default=True,
    help='Output prefix')
@optgroup.group('Other options')
@optgroup.option('--bam_offset', type=click.STRING,
    default="0,-1", show_default=True, callback=tuple_args(int),
    help='BAM file offset (enables support for other datatypes -- e.g. Tn5/ATAC)')
@optgroup.option('--n_threads', type=click.IntRange(1, cpu_count()),
    default=cpu_count(), show_default=True,
    help='Number of processors to use')
@optgroup.option('--batch_size', type=click.INT,
    default=100, show_default=True,
    help='Batch size of intervals to process')
def run(interval_file,
        bam_file,
        fasta_file,
        bias_model_file=None,
        min_qual=1,
        keep_dups=False,
        keep_qcfail=False,
        bam_offset=(0, -1),
        half_win_width=5,
        n_threads=cpu_count(),
        batch_size=100,
        outfile='dm.json'):
    """Learn a negative binomial dispersion model 
    
    This creates expected cleavage counts from data corrected for intrinsic 
    sequence preferences and then builds a dispersion model used for footprint
    detection and analysis.

    INTERVAL_FILE is a BED-formatted file contained genomic regions to be analyzed.
    BAM_FILE is the path to a BAM-format tag alignment file. FASTA_FILE is the path
    to genome FASTA file (requires associated FASTA index in same folder (see 
    documentation on how to create an index).

    Outputs a JSON-formated dispersion model
    """

    proc_kwargs = {
        "min_qual": min_qual,
        "remove_dups": ~keep_dups,
        "remove_qcfail": ~keep_qcfail,
        "offset": bam_offset,
        "half_win_width": half_win_width,
    }

    # Validate and load input files
    logger.info("Validating input files")
    try:
        verify_bam_file(bam_file)
        verify_fasta_file(fasta_file)

        # Load bias model (if specified), otherwise use the uniform model
        if bias_model_file:
            logger.info(f"Loading bias model from file {bias_model_file}")
            bm = bias.kmer_model(bias_model_file)
        else:
            logger.info(f"No bias model file specified -- using uniform model")
            bm = bias.uniform_model()

    except IOError as e:
        logger.critical(e)
        raise click.Abort()

    # Initiate data processor
    dl = expected_counts(interval_file, bam_file, fasta_file, bm, **proc_kwargs)
    dl_iter = dl.batch_iter(batch_size=batch_size, collate_fn=numpy_collate_concat, num_workers=n_threads)

    hist_size = (200, 1000)
    hist = np.zeros(hist_size, dtype=int)

    with logging_redirect_tqdm():

        for cnts in tqdm(dl_iter, colour='#cc951d'):
            for i in range(cnts.shape[0]):
                try:
                    # expected, observed
                    hist[int(cnts[i,0]),int(cnts[i, 1])] += 1
                # ignore counts bigger than histogram bounds
                except IndexError:
                    pass

    logger.info("Learning dispersion model")
    
    model = dispersion.learn_dispersion_model(hist)

    logger.info("Writing dispersion model to {}".format(outfile))

    try:
        with open(outfile, "w") as output_filehandle:
            print(dispersion.write_dispersion_model(model), file=output_filehandle)

    except IOError as e:
        logger.critical(e)
        raise click.Abort()