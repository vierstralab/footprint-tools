import sys

import argh
from argh.decorators import named, arg

from tqdm import tqdm
import multiprocessing as mp

import numpy as np
import pandas as pd
import scipy.stats
import pysam

from genome_tools import bed

from footprint_tools.modeling import dispersion
from footprint_tools.stats import bayesian, segment

from footprint_tools.cli.utils import chunkify

import logging
logger = logging.getLogger(__name__)

def read_func(tabix_files, intervals, queue):
    """Reads TABIX files and outputs to a multiprocessing pool queue
    """

    tabix_handles = [pysam.TabixFile(f) for f in tabix_files]
    n_datasets = len(tabix_handles)

     # Write to input queue
    for interval in intervals:

        l = len(interval)

        obs = np.zeros((n_datasets, l), dtype = np.float64)
        exp = np.zeros((n_datasets, l), dtype = np.float64)
        fdr = np.ones((n_datasets, l), dtype = np.float64)
        w = np.zeros((n_datasets, l), dtype = np.float64)

        for i, tabix in enumerate(tabix_handles):

            try:
                for row in tabix.fetch(interval.chrom, interval.start, interval.end, parser=pysam.asTuple()):
                    j = int(row[1])-interval.start
                    exp[i, j] = np.float64(row[3])
                    obs[i, j] = np.float64(row[4])
                    fdr[i, j] = np.float64(row[7])
                    w[i, j] = 1
            except:
                pass

        queue.put( (interval, exp, obs, fdr, w) )

        # Stop memory from getting out of control in the processing
        while queue.qsize() > 100:
            pass

    [handle.close() for handle in tabix_handles]

def process_func(disp_models, beta_priors, read_q, write_q, fdr_cutoff=0):
    """Listens to read queue and computes posteriors
    """
    disp_models=[]
    beta_priors

    while 1:
        data = read_q.get()
        # If sentinel value is detected, return thread
        if data == None:
            read_q.task_done()
            break

        (interval, exp, obs, fdr, w) = data

        prior = bayesian.compute_prior_weighted(fdr, w, cutoff=fdr_cutoff)   
        scale = bayesian.compute_delta_prior(obs, exp, fdr, beta_priors, cutoff=fdr_cutoff)

        ll_on = bayesian.log_likelihood(obs, exp, disp_models, delta = scale, w=3) 
        ll_off = bayesian.log_likelihood(obs, exp, disp_models, w=3)

        # Compute posterior
        post = -bayesian.posterior(prior, ll_on, ll_off)
        post[post <= 0] = 0.0

        read_q.task_done()

        write_q.put( (interval, post) )

def write_func(q, outfile, total, log_post_cutoff=0):
    """Writer function
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

            interval, post = data
            
            # Write ouput; filter out positions by posterior cutoff
            for i in np.where(np.nanmax(post, axis=0) > log_post_cutoff)[0]:
                out = f'{interval.chrom}\t{interval.start:d}\t{interval.start+1:d}\t'
                out += '\t'.join(map(str, post[:,i]))
                print(out, file=outfile)

            q.task_done()

            progress_bar.update(1)

@named('posterior')
@arg('sample_data_file',
    type=str,
    help='')
@arg('interval_file',
    type=str,
    help='')
@arg('--fdr_cutoff',
    type=float,
    default=0.05,
    help='')
@arg('--post_cutoff',
    type=float,
    default=0.2,
    help='')
@arg('--n_threads',
    type=int,
    help='Number of processors to use (min=4)',
    default=max(4, mp.cpu_count()))
def run(sample_data_file, interval_file, fdr_cutoff=0.05, post_cutoff=0.2, n_threads=max(4, mp.cpu_count())):
    """Compute the posterior probability of cleavage data

    Output:
        A bedGraph-like file written to `stdout`
            contig start start+1 -log(posterior) ... N_samples
    
        Note that output data is not sorted -- pipe to `sort -k1,1 -k2,2n` for sorted output
    """

    logger.info(f"Reading sample data file ({sample_data_file}) and verifying inputs")

    sample_data = pd.read_table(sample_data_file, header=0, comment='#')

    tabix_files = []
    disp_models = []
    beta_priors = np.ones((len(sample_data), 2))

    # Load and parse input files
    try:

        for i, sample in enumerate(sample_data.itertuples()):

            # Test whether we can successfully open the TABIX file; throws exception on errors
            f = pysam.TabixFile(sample["tabix_file"]).close()
            tabix_files.append(sample["tabix_file"])

            # Load dispersion model
            dm = dispersion.load_dispersion_model(sample["dispersion_model"])
            disp_models.append(dm)

            # Read beta prior parameters
            with open(sample["beta_prior_file"], 'r') as f:
                params = f.readline().strip().split('\t')
                if len(params)!=2:
                    raise ValueError(f"Beta prior file malformed -- must contain 2 columns ({sample['beta_prior_file']}")
                beta_priors[i,:] = np.array(params, dtype = np.float64)

    except (IOError, ValueError) as e:
        logger.critical(e)
        return 1

    # Load intervals file
    intervals = list(bed.bed3_iterator(open(interval_file)))
    
    logger.info(f"BED file contains {len(intervals):,} regions")
    
    # Processing queues
    read_q = mp.JoinableQueue()
    write_q = mp.JoinableQueue()

    #
    read_procs = []
    for i, chunk in enumerate(chunkify(intervals, 2)):
        p = mp.Process(target=read_func, args=(sample_data["tabix_file"].tolist(), chunk, read_q))
        read_procs.append(p)

    logger.info(f"Using {len(read_procs)} threads to read footprint statistics")
    [p.start() for p in read_procs]

    #
    process_procs = []
    process_kwargs = {
        "fdr_cutoff": fdr_cutoff,
    }

    for i in range(n_threads-3):
        p = mp.Process(target=process_func, args=(disp_models, beta_priors, read_q, write_q), kwargs=process_kwargs)
        process_procs.append(p)

    logger.info(f"Using {len(process_procs)} threads to compute posteriors")
    [p.start() for p in process_procs]

    #
    write_kwargs = {
        "log_post_cutoff": -np.log(post_cutoff)
    }
    write_proc = mp.Process(target=write_func, args=(write_q, len(intervals)), kwargs=write_kwargs)
    write_proc.start()

    try:
        # Wait for readers to finish
        [p.join() for p in read_procs]

        # Block until all remaining regions are processed
        read_q.join() # wait till read queue is empty
        [read_q.put(None) for i in range(len(process_procs))] # sends sentinal signal
        [p.join() for p in process_procs] # block until thread exits
        
        # Block until all remaining regions are written to output and thread exits
        write_q.join() # wait till write queue is empty
        write_q.put(None) # sends sentinal signal
        write_proc.join() # block until thread exits

        logger.info("Finished computing and writing footprint posteriors!")

    except KeyboardInterrupt:
        [p.terminate() for p in read_procs]
        [p.terminate() for p in process_procs]
        write_proc.terminate()

    return 0