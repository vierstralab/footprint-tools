import sys

import argh
from argh.decorators import named, arg

from tqdm import tqdm
import multiprocessing as mp

import numpy as np
import scipy as sp
import pysam

from genome_tools import bed

from footprint_tools import cutcounts
from footprint_tools.modeling import bias, predict, dispersion
from footprint_tools.stats import fdr, windowing

from footprint_tools.cli.utils import chunkify, tuple_ints

import logging
logger = logging.getLogger(__name__)

def read_func(bam_file, fasta_file, bm, dm, intervals, q, **kwargs):
    """Reads BAM file, computes expected cleavages and associated
        statistics and outputs to a multiprocessing pool queue
    """
    bam_kwargs = { k:kwargs.pop(k) for k in ["min_qual", "remove_dups", "remove_qcfail", "offset"] }
    bam_reader = cutcounts.bamfile(bam_file, **bam_kwargs)
    fasta_reader = pysam.FastaFile(fasta_file)

    predict_kwargs = { k:kwargs.pop(k) for k in ["half_window_width", "smoothing_half_window_width", "smoothing_clip"] }
    predictor = predict.prediction(bam_reader, fasta_reader, bm, **predict_kwargs)

    # Args used for FDR sampling procedure
    seed = kwargs.pop("seed")
    if seed:
        random.seed(seed)
        np.random.seed(seed)
    fdr_shuffle_n = kwargs.pop("fdr_shuffle_n") 

    win_pvals_func = lambda z: windowing.stouffers_z(np.ascontiguousarray(z), 3)
    
    # Read and process each region
    for interval in intervals:

        obs, exp, win = predictor.compute(interval)

        obs = obs['+'][1:] + obs['-'][:-1]
        exp = exp['+'][1:] + exp['-'][:-1]

        n = len(obs)

        if dm:
            try:
                pvals = dm.p_values(exp, obs)
                _, pvals_null = dm.sample(exp, fdr_shuffle_n)

                win_pvals = win_pvals_func(pvals)
                win_pvals_null = np.apply_along_axis(win_pvals_func, 0, pvals_null)
                
                efdr = fdr.emperical_fdr(win_pvals_null, win_pvals)

            except Exception as e:
                pvals = win_pvals = efdr = np.ones(n)

            finally:
                stats = np.column_stack((exp, obs, -np.log(pvals), -np.log(win_pvals), efdr))
        else:
            stats = np.column_stack((exp, obs))

        q.put( (interval, stats) )

        while q.qsize() > 100:
            pass

def write_func(q, outfile, total):
    """Function to write output from multiple threads to a single file
    """    
    f = outfile

    progress_desc="Regions processed"
    with tqdm(total=total, desc=progress_desc, ncols=80) as progress_bar:

        while 1:
            data = q.get()

            # If sentinel value is detected, return thread
            if data == None:
                q.task_done()
                break

            interval, stats = data

            for i in range(stats.shape[0]):
                coords = "{}\t{:d}\t{:d}".format(interval.chrom, interval.start+i, interval.start+i+1)
                val_string = "\t".join( ["{:0.4f}".format(val) for val in stats[i,:]])
                print(coords + "\t" + val_string, file=f)

            q.task_done()

            progress_bar.update(1)

@named('detect')
@arg('interval_file', 
    help='File path to BED file contain regions to analyzed')
@arg('bam_file', 
    help='Path to BAM-format tag alignment file')
@arg('fasta_file',
    help='Path to genome FASTA file (requires associated FASTA '
        'index in same folder; see documentation on how to create '
        'an index)')
@arg('--bias_model_file',
    type=str,
    default=None,
    help='Use a k-mer model for local bias (supplied by file). '
        'If argument is not provided the model defaults to uniform '
        'sequence bias.')
@arg('--dispersion_model_file',
    type=str,
    default=None,
    help='Dispersion model for negative binomial tests. If argument '
        'is not provided then no stastical output is provided. File is in '
        'JSON format and generated using the command learn_dm')
@arg('--min_qual',
    type=int,
    default=1,
    help='Ignore reads with mapping quality lower than this threshold')
@arg('--remove_dups',
    help='Remove duplicate reads from consideration',
    default=False)
@arg('--remove_qcfail',
    help='Remove QC-failed reads from consideration',
    default=False)
@arg('--bam_offset',
    type=tuple_ints,
    default=(0, -1),
    help='BAM file offset (enables support for other datatypes -- e.g., Tn5/ATAC)')
@arg('--half_win_width',
    help='Half window width to apply bias model',
    default=5)
@arg('--smooth_half_win_width',
    type=int,
    help='Half window width to apply smoothing model. When set to '
        '0, no smoothing is applied.',
    default=50)
@arg('--smooth_clip',
    type=float,
    help='Fraction of bases to clip when computing trimmed mean in the smoothing window',
    default=0.01)
@arg('--fdr_shuffle_n',
    type=int,
    default=50,
    help='Number of times to shuffle data for FDR calculation')
@arg('--seed',
    default=None,
    help='Seed for random number generation')
@arg('--n_threads',
    help='Number of processors to use (min=2)',
    default=max(2, mp.cpu_count()))
def run(interval_file,
        bam_file,
        fasta_file,
        bias_model_file=None,
        dispersion_model_file=None,
        min_qual=1,
        remove_dups=False,
        remove_qcfail=False,
        bam_offset=(0, -1),
        half_win_width=5,
        smooth_half_win_width=50,
        smooth_clip=0.01,
        fdr_shuffle_n=50,
        seed=None,
        n_threads=max(2, mp.cpu_count())):
    """Compute per-nucleotide cleavage deviation statistics	

    Output:
        bedGraph file written to `stdout`:
            contig start start+1 obs exp -log(pval) -log(winpval) fdr

        Note that output data is not sorted -- pipe to `sort -k1,1 -k2,2n` for sorted output
    """
    intervals = list(bed.bed3_iterator(open(interval_file)))
    
    logger.info(f"BED file contains {len(intervals):,} regions")

    proc_kwargs = {
        "min_qual": min_qual,
        "remove_dups": remove_dups,
        "remove_qcfail": remove_qcfail,
        "offset": bam_offset,
        "half_window_width": half_win_width,
        "smoothing_half_window_width": smooth_half_win_width,
        "smoothing_clip": smooth_clip,
        "fdr_shuffle_n": fdr_shuffle_n,
        "seed": seed,
    }

    # Load bias model (if specified), otherwise use the uniform model
    if bias_model_file:
        logger.info(f"Loading bias model from file {bias_model_file}")
        bm = bias.kmer_model(bias_model_file)
    else:
        logger.info(f"No bias model file specified -- using uniform model")
        bm = bias.uniform_model()

    # Load dispersion model (if specified)
    if dispersion_model_file:
        dm = dispersion.load_dispersion_model(dispersion_model_file)
        logger.info(f"Loaded dispersion model from file {dispersion_model_file}")

    else:
        logger.info(f"No dispersion model file specified -- will not be reporting base-level statistics")
        dm = None

    q = mp.JoinableQueue()

    read_procs = []
    for i, chunk in enumerate(chunkify(intervals, max(1, n_threads-1))):
           p = mp.Process(target=read_func, args=(bam_file, fasta_file, bm, dm, chunk, q), kwargs=proc_kwargs)
           read_procs.append(p)

    write_proc = mp.Process(target=write_func, args=(q, sys.stdout, len(intervals)))

    logger.info(f"Using {len(read_procs)} threads to compute footprint statistics")
    logger.info("Using 1 thread to write footprint statistics")

    [p.start() for p in read_procs]
    write_proc.start()

    try:
        [p.join() for p in read_procs]
        q.join() # block until queue is empty after processing is done
        q.put(None) # sends sentinal signal
        write_proc.join() # wait for writer proc to return
        logger.info("Finished computing and writing footprint statistics!")
    except KeyboardInterrupt:
        [p.terminate() for p in read_procs]
        write_proc.terminate()
    
    return 0
