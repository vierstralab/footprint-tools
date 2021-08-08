import argh
import logging

import footprint_tools
import footprint_tools.cli.learn_dispersion as learn_dm
import footprint_tools.cli.find_fps as find_fps
import footprint_tools.cli.learn_beta as learn_beta

epilog = """Written by Jeff Vierstra (jvierstra@altius.org) (2015-2021). 

See http://github.com/jvierstra/footprint-tools for extended documentation.

Software licensed under GNU General Public License version 3."""

# TODO: Color formating; load from a .conf file
logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p', level=logging.INFO)
logger = logging.getLogger(__name__)

def main()
    parser = argh.ArghParser(epilog=epilog)
    parser.add_commands([
        learn_dm.run,
        find_fps.run,
        learn_beta.run
    ])
    parser.add_argument('--version',
                    action='version',
                    version='%(prog)s ' + footprint_tools.__version__)
    argh.dispatch(parser)