import os
import sys

from importlib.metadata import version as _pkg_version
from datetime import datetime
extensions = [
    "sphinx.ext.intersphinx",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.extlinks",
    "sphinx.ext.todo",
    "sphinx.ext.mathjax",
    "IPython.sphinxext.ipython_directive",
    "IPython.sphinxext.ipython_console_highlighting",
    "matplotlib.sphinxext.plot_directive",
    "numpydoc",
    "sphinx_panels",
    "sphinx_click",
    "nbsphinx"
]

autosectionlabel_prefix_document = True
numpydoc_show_class_members = False

intersphinx_mapping = {
    'pysam': ('http://pysam.readthedocs.org/en/latest', None),
    'numpy': ('http://numpy.readthedocs.org/en/latest', None),
}


sys.path.insert(0, os.path.abspath('../..'))

project = 'footprint-tools'
copyright = f'2015-{datetime.now().year}, Jeff Vierstra'
author = 'Jeff Vierstra'
version = _pkg_version('footprint_tools')  # docs build installs the package
templates_path = ['_templates']
source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}
master_doc = 'index'
pygments_style = 'sphinx'

html_theme = 'pydata_sphinx_theme'
html_theme_options = {
    "external_links": [],
    "github_url": "https://github.com/vierstralab/footprint-tools",
}

html_static_path = ['_static']
html_css_files = ['extra.css']



