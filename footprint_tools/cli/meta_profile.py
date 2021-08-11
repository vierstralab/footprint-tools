import argh
from argh.decorators import named, arg

from tqdm import tqdm
import multiprocessing as mp
from functools import partial

import numpy as np

import pysam
import pyBigWig

from footprint_tools.cli.utils import chunkify

import logging
logger = logging.getLogger(__name__)

def read_tabix_func(tabix_file, intervals, col=4, dtype=int):
    """Read data from tabix
    """
    
    m = len(intervals[0]) # width of intervals
    n = len(intervals) # total intervals

    res = np.array((n, m), q=None, dtype=dtype)

    tabix_handle = pysam.TabixFile(tabix_file)

    for i, interval in enumerate(intervals):

        for row in tabix_handle.fetch(interval.chrom, interval.start, interval.end, parser=pysam.asTuple()):
            j = int(row[1])-interval.start
            res[i, j] = row[col]

            # Update progress queue
            if q:
                q.put(1)
    
    tabix_handle.close()

    return res

def read_bigwig_func(bigwig_file, intervals, q, dtype=np.float):
    """Read data from bigwig
    """
    
    m = len(intervals[0]) # width of intervals
    n = len(intervals) # total intervals

    res = np.array((n, m), dtype=dtype)

    bigwig_handle = pyBigWig.open(bigwig_file)

    for i, interval in enumerate(intervals):
        res[i,:] = bigwig_file.values(interval.chrom, interval.start, interval.end, numpy=True).astype(dtype)
        
        # Update progress queue	
        q.put(1)

    bigwig_handle.close()

    return res

def listener(q, total):
    """Listens to a queue and updates a progess bar
    """
    progress_desc="Elements processed"
    with tqdm(total=total, desc=progress_desc, ncols=80) as progress_bar:
        while 1:
            try:
                progress_bar.update(q.get(timeout=1))
                if progress_bar.n >= total:
                    break
            except:
                continue

@named('meta_profile')
@arg('interval_file',
    type=str,
    help='')
@arg('signal_file',
    type=str,
    help='')
@arg('file_type',
    help='')
@arg('n_threads',
    type=int,
    default=mp.cpu_count(),
    help='')
def run(interval_file, tabix_file, n_threads):
    """Generate metaprofile data
    """

    # Load intervals file
    intervals = list(bed.bed3_iterator(open(interval_file)))
    
    logger.info(f"BED file contains {len(intervals):,} elements")

    pool = mp.Pool(n_threads)
    q = mp.Queue()

    # func = partial()
    return 0
