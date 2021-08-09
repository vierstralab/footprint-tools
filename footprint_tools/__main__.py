import argh

import functools

import footprint_tools
import footprint_tools.cli.learn_dispersion as learn_dm
import footprint_tools.cli.find_fps as find_fps
import footprint_tools.cli.learn_beta as learn_beta
import footprint_tools.cli.plot_dm as plot_dm
import footprint_tools.cli.posterior as posterior

import pkg_resources
import logging
import logging.config

# Load logging config from "logging.conf" file
logging.config.fileConfig(pkg_resources.resource_filename(__name__, "logging.conf"))
logger = logging.getLogger(__name__)

epilog = """Written by Jeff Vierstra (jvierstra@altius.org) (2015-2021). 

See http://github.com/jvierstra/footprint-tools for extended documentation.

Software licensed under GNU General Public License version 3."""

def main():
	loggers = [logging.getLogger(name) for name in logging.root.manager.loggerDict]
	for l in loggers:
		print(l.handlers)

	parser = argh.ArghParser(epilog=epilog)
	parser.add_commands([
		learn_dm.run,
		find_fps.run,
		learn_beta.run,
		plot_dm.run,
		posterior.run,
	])
	parser.add_argument('--version',
					action='version',
					version='%(prog)s ' + footprint_tools.__version__)
	argh.dispatch(parser)