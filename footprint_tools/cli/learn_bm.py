import argh
from argh.decorators import named, arg

import pysam

from footprint_tools.cli.utils import tuple_ints

import logging
logger = logging.getLogger(__name__)

@named('learn_bm')
@arg('bam_file',
    type=str,
    help='')
@arg('fasta_file',
    type=str,
    help='')
@arg('mappability_file',
    type=str,
    help='')
@arg('k',
    type=int,
    default=6
    help='')
@arg('mid',
    type=int,
    default=3
    help='')
@arg('bam_offset',
    type=tuple_ints,
    default=(0,-1),
    help='')
def run(bam_file, fasta_file, mappability_file, k=6, mid=3, offset=(0, -1)):
    """Learn a bias model
    """
    raise NotImplementedError
