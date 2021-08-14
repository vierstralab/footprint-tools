import os
import sys

import distutils.util

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
version = '1.3.2'
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
    "github_url": "https://github.com/jvierstra/footprint-tools",
}

html_static_path = ['_static']
html_css_files = ['extra.css']


# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
_libdir = "../../build/lib.%s-%s.%s" % (distutils.util.get_platform(),
                                        sys.version_info[0], sys.version_info[1])
print(_libdir)
if os.path.exists(_libdir):
    sys.path.insert(0, os.path.abspath(_libdir))


mathjax_path="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML"
