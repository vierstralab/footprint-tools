import os
import sys

import sphinx_rtd_theme

extensions = [
    "sphinx_rtd_theme",
    "sphinx.ext.intersphinx",
    "sphinx.ext.autodoc",
    "sphinx.ext.viewcode",
]


intersphinx_mapping = {
    'pysam': ('http://pysam.readthedocs.org/en/latest', None),
}


sys.path.insert(0, os.path.abspath('../..'))

project = 'footprint-tools'
copyright = '2020, Jeff Vierstra'
author = 'Jeff Vierstra'
release = '1.0.0'
templates_path = ['_templates']
source_suffix = '.rst'
master_doc = 'index'
pygments_style = 'sphinx'
html_theme = "sphinx_rtd_theme"
html_static_path = ['_static']
