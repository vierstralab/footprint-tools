import argh
from argh.decorators import named, arg

import numpy as np
import pandas as pd

from genome_tools import genomic_interval

from footprint_tools import cutcounts
from footprint_tools.data.dataset import dataset

from tqdm import tqdm

import logging
logger = logging.getLogger(__name__)

class profile_loader(dataset):
    """Class that enables iterable and parellelized loading and
    processing of data
    """
    
    def __init__(self, interval_file, bam_file):
        self.bam_file = bam_file
        self.intervals = pd.read_table(interval_file, header=None)

        self.bam_reader = None

    def __len__(self):
        return len(self.intervals)

    def __getitem__(self, index):
        
        # Open file handlers on first call. This avoids problems when
        # parallel processing data with non-thread safe code (i.e., pysam)
        if not self.bam_reader:
            self.bam_reader = cutcounts.bamfile(self.bam_file)

        chrom, start, end = self.intervals.iat[index, 0], self.intervals.iat[index, 1], self.intervals.iat[index, 2]

        interval = genomic_interval(chrom, start, end)

        counts = self.bam_reader[interval]
        counts.pop('fragments')

        return counts

@named('meta_profile')
@arg('interval_file',
    type=str,
    help='')
@arg('bam_file',
    type=str,
    help='')
@arg('n_threads',
    type=int,
    default=8,
    help='')
def run(interval_file, bam_file, n_threads):
    """Generate metaprofile data
    """

    ds = profile_loader(interval_file, bam_file)
    #print(ds[0])
    
    pos = [x["+"] for x in tqdm(ds.batch_iter(batch_size=1, num_workers=n_threads))]
    