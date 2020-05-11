import os
import sys

import sphinx_rtd_theme

extensions = [
    "sphinx_rtd_theme",
]


sys.path.insert(0, os.path.abspath('../..'))
extensions = [
'sphinx.ext.autodoc',]

project = 'footprint-tools'
copyright = '2020, Jeff Vierstra'
author = 'Jeff Vierstra'
version = '1'
release = '0.0'
templates_path = ['_templates']
source_suffix = '.rst'
master_doc = 'index'
pygments_style = 'sphinx'
html_theme = "sphinx_rtd_theme"
html_static_path = ['_static']
