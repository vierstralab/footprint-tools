import argh
from argh.decorators import named, arg

import logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

@named('autocorr')
@arg('tabix_file',
	type=str,
	help='')
@arg('--range',
	type=tuple_int,
	default=(-25, 25),
	help='')
def run(tabix_file, range):
	"""Calculate signal autocorrelation for genomic intervals"""
	return 0