import argh
from argh.decorators import named, arg

from footprint_tools.cli.utils import tuple_ints

import logging
logger = logging.getLogger(__name__)

@named('autocorr')
@arg('tabix_file',
    type=str,
    help='')
@arg('--range',
    type=tuple_ints,
    default=(-25, 25),
    help='')
def run(tabix_file, range):
    """Calculate nucleotide-level signal autocorrelation for genomic intervals
    """
    raise NotImplementedError
