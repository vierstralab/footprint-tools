from pkg_resources import get_distribution

__version__ = get_distribution('footprint_tools').version
__all__ = ["cutcounts"]