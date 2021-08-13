import click
from click_option_group import optgroup

logger = logging.getLogger(__name__)    

@click.command('diff_test')
def run():
    """Base-pair resolution test for differential DNase I cleavage
    """
    raise NotImplementedError
