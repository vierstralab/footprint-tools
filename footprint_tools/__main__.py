import click

# Set up console logging from package file
# This configure root logger which is inherited 
# by all modules/submodules
import pkg_resources
import logging, logging.config
logging.config.fileConfig(pkg_resources.resource_filename(__name__, "logging.conf"))

logger = logging.getLogger(__name__)

import footprint_tools
import footprint_tools.cli.learn_bm as learn_bm
import footprint_tools.cli.learn_dm as learn_dm
import footprint_tools.cli.detect as detect
import footprint_tools.cli.learn_beta as learn_beta
import footprint_tools.cli.plot_dm as plot_dm
import footprint_tools.cli.post as post

epilog = """See http://github.com/jvierstra/footprint-tools for extended documentation.

Citation: Vierstra, J., Lazar, J., Sandstrom, R. et al. Global reference mapping of human transcription factor footprints. Nature 583, 729–736 (2020)

Written by Jeff Vierstra (jvierstra@altius.org) (2015-2021). Software licensed under GNU General Public License version 3."""

@click.group(epilog=epilog)
@click.version_option(version=footprint_tools.__version__)
def main():
    """footprint_tools: analysis of digital genomic footprints 
    
    This software packages enables de novo detection and analysis of genomic footprints from DNase 
    I data. The underying model simulates expected cleavage rates using a 6-mer DNase I cleavage
    preference model combined with density smoothing. Statistical significance of per-nucleotide 
    cleavages are computed from a series of emperically fit negative binomial distributions. In addition
    to footprint detection in a single, isolated dataset, this package has a statistical framework
    to jointly analyze 100s to 1000s of datasets in unison.
    """
    pass

main.add_command(learn_bm.run)
main.add_command(learn_dm.run)
main.add_command(detect.run)
main.add_command(learn_beta.run)
main.add_command(plot_dm.run)
main.add_command(post.run)