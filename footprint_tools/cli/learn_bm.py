import click

from footprint_tools.cli.utils import tuple_args

from tqdm import tqdm
from tqdm.contrib.logging import logging_redirect_tqdm

import logging
logger = logging.getLogger(__name__)

@click.command(name='learn_bm')
@click.argument('bam_file')
@click.argument('fasta_file')
@click.argument('mappability_file')
@click.option('-k', type=click.INT,
    default=6, show_default=True,
    help='Width of k-mer model')
@click.option('--mid', type=click.INT,
    default=3, show_default=True,
    help='Cleavage position of k-mer model relative to cut site')
@click.option('--bam_offset', type=click.STRING,
    default="0,-1", show_default=True, callback=tuple_args(int),
    help='BAM file offset (enables support for other datatypes -- e.g. Tn5/ATAC)')
def run(bam_file, fasta_file, mappability_file, k=6, mid=3, offset=(0, -1)):
    """Learn a sequence bias model
    """
    raise click.UsageError("Command not yet implemented")
