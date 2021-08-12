from typing import IO
from footprint_tools.data.processor import exception_wrapper
import click
import pysam

def tuple_args(value_type=int):
    def _parse_tuple(ctx, params, value):
        """Function to parse a tuple from command line"""
        try:
            items = tuple(map(value_type, value.split(',')))
            assert len(items) == 2
        except:
            raise click.BadOptionUsage(f'needs to be a comma-delimited tuple of type {value_type.__name__}')
        return items
    return _parse_tuple

def list_args(value_type=int):
    def _parse_list(ctx, params, value):
        """Function to parse a list from command line"""
        try:
            items = list(map(value_type, value.split(',')))
        except:
            raise click.BadOptionUsage(f'needs to be a comma-delimited list of type {value_type.__name__}')
        return items
    return _parse_list

def get_kwargs(keys, kwargs):
    return {k:kwargs[k] for k in keys if k in kwargs}

def verify_bam_file(fn):
    """Tries to open a file, raises IOError with problems"""
    try:
        pysam.AlignmentFile(fn).close()
    except IOError:
        raise IOError(f"No such file: {fn}")
    except ValueError:
        raise IOError(f"BAM-index not found for {fn}")

def verify_tabix_file(fn):
    """Tries to open a file, raises IOError with problems"""
    try:
        pysam.TabixFile(fn).close()
    except IOError:
        raise IOError(f"Not such file: {fn}")
    except ValueError:
        raise IOError(f"TABIX-index not found for {fn}")

def verify_fasta_file(fn):
    """Tries to open a file, raises IOError with problems"""
    try:
        pysam.FastaFile(fn).close()
    except IOError:
        raise IOError(f"Not such file: {fn}")
    except ValueError:
        raise IOError(f"FASTA-index not found for {fn}")