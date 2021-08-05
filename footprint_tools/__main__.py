import argh
import logging

import footprint_tools
import footprint_tools.cli.learn_dispersion as learn_dm

epilog = """Written by Jeff Vierstra (jvierstra@altius.org) (2015-2021). 

See http://github.com/jvierstra/footprint-tools for extended documentation.

Licensed under GNU General Public License version 3."""

logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p', level=logging.INFO)
logger = logging.getLogger(__name__)

def main()
    parser = argh.ArghParser(epilog=epilog)
    parser.add_commands([
        learn_dm.run,
    ])
    parser.add_argument('--version',
                    action='version',
                    version='%(prog)s ' + footprint_tools.__version__)
    argh.dispatch(parser)