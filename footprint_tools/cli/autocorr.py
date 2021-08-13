import click
from click_option_group import optgroup

from footprint_tools.cli.utils import tuple_args

import logging
logger = logging.getLogger(__name__)

@click.command('autocorr')
@click.argument('tabix_file')
@click.option('--range', type=str, 
    default="-25,25" callback=tuple_args(int),
    help='')
def run(tabix_file, range):
    """Calculate nucleotide-level signal autocorrelation for genomic intervals
    """
    raise NotImplementedError
