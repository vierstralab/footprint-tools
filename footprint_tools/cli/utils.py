import click
from argh.exceptions import CommandError

def tuple_args(ctx, param, value):
    """
    Function to parser a tuple of integers from command line
    """
    items = value.split(',')
    if len(items) > 2:
        raise click.BadOptionUsage('needs to be a comma-delimited tuple')

    return items

def list_args(arg_type=int):
    """
    Function to parser a list of integers from command line
    """
    def parse_list(s):
        try:
            return list(map(arg_type, s.split(',')))
        except Exception as e:
            raise CommandError("Argument must be a comma-delimited list -- i.e., 0,1,2,3")
    return parse_list

def fstr(template):
    return eval(f"f'{template}'")

def get_kwargs(keys, kwargs):
    return {k:kwargs[k] for k in keys if k in kwargs}

def chunkify(l, n):
    """Splits an iterable list into n chunks"""
    return [l[i::n] for i in range(n)]
