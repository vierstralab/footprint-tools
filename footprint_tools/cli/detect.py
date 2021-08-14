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

from footprint_tools import cutcounts
from footprint_tools.modeling import bias, predict, dispersion
from footprint_tools.stats import fdr, windowing

from footprint_tools.cli.utils import (list_args, tuple_args, get_kwargs, 
                                        verify_bam_file, verify_fasta_file, write_output_header,
                                        write_stats_to_output, write_segments_to_output)

from tqdm import tqdm
from tqdm.contrib.logging import logging_redirect_tqdm

import logging
logger = logging.getLogger(__name__)

# kill numpy warnings
np.seterr(all="ignore")

class deviation_stats(dataset):
    """Class that computes per-nucleotide cleavage
    deviation statistics
    """
    def __init__(self, interval_file, bam_file, fasta_file, bm, dm, **kwargs):

        self.intervals = pd.read_table(interval_file, header=None)

        logger.info(f"BED file contains {len(self.intervals):,} regions")

        self.bam_file = bam_file
        self.fasta_file = fasta_file
        self.bm = bm
        self.dm = dm

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
                'smoothing_half_win_width',
                'smoothing_clip'
            ],
            kwargs)

        self.fdr_shuffle_n = kwargs['fdr_shuffle_n']
        self.seed = kwargs['seed']

        self.counts_extractor = None
        self.fasta_extractor = None
        self.count_predictor = None

        self.win_pval_fn = lambda z: windowing.stouffers_z(np.ascontiguousarray(z), 3)

    def cleanup(self):
        """Needs to be implemented"""
        pass

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

        assert len(obs) == len(exp)
        n = len(obs)

        if self.dm:
            try:
                pvals = self.dm.p_values(exp, obs)
                win_pvals = self.win_pval_fn(pvals)

                _, pvals_null = self.dm.sample(exp, self.fdr_shuffle_n)
                win_pvals_null = np.apply_along_axis(self.win_pval_fn, 0, pvals_null)

                efdr = fdr.emperical_fdr(win_pvals_null, win_pvals)
            except Exception as e:
                logger.warning(f"Error computing stats for '{interval.chrom}:{interval.start}-{interval.end}'")
                pvals = win_pvals = efdr = np.ones(n) # should change to return 'nan'
            finally:
                stats = np.column_stack((exp, obs, -np.log(pvals), -np.log(win_pvals), efdr))
        else:
            stats = np.column_stack((exp, obs))

        return {
            'interval': interval,
            'stats': stats
        }

@click.command(name='detect')
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
@optgroup.group('Smoothing options')
@optgroup.option('--smooth_half_win_width', type=click.INT,
    default=50, show_default=True,
    help='Half window width to apply smoothing model. When set to '
        '0, no smoothing is applied.')
@optgroup.option('--smooth_clip', type=click.FLOAT,
    default=0.01, show_default=True,
    help='Fraction of bases to clip when computing trimmed mean in the smoothing window')
@optgroup.group("Statistics options")
@optgroup.option('--dispersion_model_file', type=click.STRING,
    help='Dispersion model for negative binomial tests. If argument '
        'is not provided then no stastical output is provided. File is in '
        'JSON format and generated using the command ``learn_dm``')
@optgroup.option('--fdr_shuffle_n', type=click.INT,
    default=100, show_default=True,
    help='Number of times to shuffle data for FDR calculation')
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
@optgroup.option('--outprefix', type=click.STRING,
    default='out', show_default=True,
    help='Output prefix')
@optgroup.option('--write_footprints', type=click.STRING, 
    default="0.001,0.01,0.05", show_default=True, callback=list_args(float),
    help='Output footprints at specified FDRs')
@optgroup.group('Other options')
@optgroup.option('--bam_offset', type=click.STRING,
    default="0,-1", show_default=True, callback=tuple_args(int),
    help='BAM file offset (enables support for other datatypes -- e.g. Tn5/ATAC)')
@optgroup.option('--seed', type=click.INT,
    help='Seed for random number generation (not currently used)')
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
        dispersion_model_file=None,
        min_qual=1,
        keep_dups=False,
        keep_qcfail=False,
        bam_offset=(0, -1),
        half_win_width=5,
        smooth_half_win_width=50,
        smooth_clip=0.01,
        fdr_shuffle_n=50,
        seed=None,
        n_threads=cpu_count(),
        batch_size=100,
        write_footprints=[0.001, 0.01, 0.05],
        outprefix='out'):
    """Compute per-nucleotide cleavage deviation statistics	

    INTERVAL_FILE is a BED-formatted file contained genomic regions to be analyzed.
    BAM_FILE is the path to a BAM-format tag alignment file. FASTA_FILE is the path
    to genome FASTA file (requires associated FASTA index in same folder (see 
    documentation on how to create an index).

    Per-nucleotde statistics are written to a bedGraph file: <outprefix>.bedgraph).
    Footprints are *optionally* written to a BED3-format file: <outprefix>.fdr<thresh>.bed).
    Footprint FDR threshold(s) can be specified using the ``--write_footprints`` option.
    Output file column definitions are written to a header within each file.
    """
    proc_kwargs = {
        "min_qual": min_qual,
        "remove_dups": ~keep_dups,
        "remove_qcfail": ~keep_qcfail,
        "offset": bam_offset,
        "half_win_width": half_win_width,
        "smoothing_half_win_width": smooth_half_win_width,
        "smoothing_clip": smooth_clip,
        "fdr_shuffle_n": fdr_shuffle_n,
        "seed": seed, # not used...yet
    }

    # Validate and load inputs
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

        # Load dispersion model (if specified)
        if dispersion_model_file:
            logger.info(f"Loading dispersion model from file {dispersion_model_file}")
            dm = dispersion.load_dispersion_model(dispersion_model_file)
        else:
            logger.info(f"No dispersion model file specified -- not be reporting base-level statistics and footprints")
            dm = None
            write_footprints = []

    except IOError as e:
        logger.critical(e)
        raise click.Abort()

    # Open output files
    try:
        # Output stats file
        output_bedgraph_file = outprefix + '.bedgraph'
        logger.info(f"Writing per-nucleotide stats to {output_bedgraph_file}")
        output_bedgraph_filehandle = open(output_bedgraph_file , 'w')
        write_output_header(["exp", "obs", "-log(pval)", "-log(winpval)", "fdr"], file=output_bedgraph_filehandle)

        # Output footprints filex
        output_bed_file_template = outprefix + '.fdr{0}.bed'
        output_bed_filehandles = {}

        if len(write_footprints) > 0:
            logger.info(f"Writing FDR thresholded footprints to {output_bed_file_template.format('{threshold}')} for threshold \u22f2 {write_footprints}")
            for t in write_footprints:
                fh =  open(output_bed_file_template.format(t), 'w')
                write_output_header(["name", "fdr"], file=fh, extra=f"thresholded @ FDR {t}")
                output_bed_filehandles.update({t:fh})

    except IOError as e:
        logger.critical(e)
        raise click.Abort()
    
    # Create data processor and iterator
    dl = deviation_stats(interval_file, bam_file, fasta_file, bm, dm, **proc_kwargs)
    dl_iter = dl.batch_iter(batch_size=batch_size, num_workers=n_threads)

    with logging_redirect_tqdm():
        
        for batch in tqdm(dl_iter, colour='#cc951d'):

            for interval, stats in zip(batch["interval"], batch["stats"]):
                # write stats
                write_stats_to_output(interval, stats, output_bedgraph_filehandle)

                # write footprints
                for thresh, fh in output_bed_filehandles.items():
                    # fdr is last column in stats array

                    write_segments_to_output(interval, stats[:,-1], thresh, file=fh, decreasing=True)

    output_bedgraph_filehandle.close()
    [f.close() for f in output_bed_filehandles.values()]