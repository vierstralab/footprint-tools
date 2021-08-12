from multiprocessing import Value
import click
from click.exceptions import BadOptionUsage
from argh.exceptions import CommandError


def tuple_args(value_type=int):
    def _parse_tuple(ctx, params, value):
        """Function to parse a tuple of integers from command line"""
        try:
            items = tuple(map(value_type, value.split(',')))
            assert len(items) == 2
        except:
            raise click.BadOptionUsage(f'needs to be a comma-delimited tuple of type {value_type.__name__}')
        return items
    return _parse_tuple

def list_args(value_type=int):
    def _parse_list(ctx, params, value):
        """Function to parse a list of integers from command line"""
        try:
            items = list(map(value_type, value.split(',')))
        except:
            raise click.BadOptionUsage(f'needs to be a comma-delimited list of type {value_type.__name__}')
        return items
    return _parse_list

def get_kwargs(keys, kwargs):
    return {k:kwargs[k] for k in keys if k in kwargs}
