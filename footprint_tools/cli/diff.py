import argh
from argh.decorators import named, arg

import logging
logger = logging.getLogger(__name__)

@named('diff_test')
def run():
    """Base-pair resolution test for differential DNase I cleavage
    """
    raise NotImplementedError
