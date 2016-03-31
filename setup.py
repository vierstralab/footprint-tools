# Copyright 2015 Jeff Vierstra

from __future__ import absolute_import, division, print_function

import os
import sys

from setuptools import find_packages, setup
from distutils.command.build_clib import build_clib
from distutils.extension import Extension
from Cython.Distutils import build_ext

import numpy as np

from glob import glob

__version__ = "1.0.0"

if sys.version_info[0] != 2 or sys.version_info[1] < 7:
    print("Package requires Python version 2.7+")
    sys.exit(1)

cephes_include = "cephes"
cephes_src = glob("cephes/*.c")
cehpes_lib = ('cephes', { 'sources': cephes_src })

ext_modules = [
	Extension("footprint_tools.modeling.predict", 
		sources = ["footprint_tools/modeling/predict.pyx"]),
	Extension("footprint_tools.modeling.dispersion", 
		sources = ["footprint_tools/modeling/dispersion.pyx"]),
	Extension("footprint_tools.stats.windowing", 
		sources = ["footprint_tools/stats/windowing.pyx"]),
	Extension("footprint_tools.stats.distributions",
		 sources = ["footprint_tools/stats/distributions.pyx"]),
	Extension("footprint_tools.stats.segment",
		 sources = ["footprint_tools/stats/segment.pyx"]),
	Extension("footprint_tools.stats.fdr.bisect",
		 sources = ["footprint_tools/stats/fdr/bisect.pyx"])
]

scripts = ["scripts/ftd-learn-dispersion-model", "scripts/ftd-compute-deviation", "scripts/ftd-compute-posterior", "scripts/ftd-learn-beta-prior"]

install_requires = ["numpy>1.10", "scipy", "pysam>=0.8.2", "pyfaidx>=0.4.2", "statsmodels", "multiprocessing"]

setup(
	name = "footprint-tools",
	version = __version__,
	description = "",
	long_description = "",
	author = "Jeff Vierstra",
	author_email = "jvierstra@altiusinstitute.org",
	zip_safe = False,
	packages =  find_packages(),
	libraries = [cehpes_lib],
    ext_modules = ext_modules,
    include_dirs=[np.get_include(), cephes_include],
    cmdclass = {'build_clib': build_clib, 'build_ext': build_ext},
    install_requires = install_requires,
    scripts = scripts
)
