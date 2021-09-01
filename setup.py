# Copyright (C) 2015-2021 Jeff Vierstra (jvierstra@altius.org)

import os
import sys
import pathlib
from glob import glob

from setuptools import find_packages, setup
from distutils.command.build_clib import build_clib
from distutils.extension import Extension
from distutils import log
from Cython.Distutils import build_ext

import numpy as np

if sys.version_info[0] != 3 or sys.version_info[1] < 6:
    print("Package requires Python version 3.6+")
    sys.exit(1)

hcephes_include = "hcephes/include"
hcephes_src = glob('hcephes/src/**/*.c', recursive=True)
hcehpes_lib = ('hcephes', { 'sources': hcephes_src })

# see MANIFEST.in -- a slight hack to include all of the header files in the source distrubution

modules = [
    dict(name="footprint_tools.modeling.predict", sources=["footprint_tools/modeling/predict.pyx"]),
    dict(name="footprint_tools.modeling.dispersion", sources=["footprint_tools/modeling/dispersion.pyx"]),
    dict(name="footprint_tools.stats.utils", sources=["footprint_tools/stats/utils.pyx"]),
    dict(name="footprint_tools.stats.windowing", sources = ["footprint_tools/stats/windowing.pyx"]),
    dict(name="footprint_tools.stats.distributions.nbinom", sources=["footprint_tools/stats/distributions/nbinom.pyx"]),
    dict(name="footprint_tools.stats.distributions.invchi2", sources=["footprint_tools/stats/distributions/invchi2.pyx"]),
    dict(name="footprint_tools.stats.distributions.t", sources=["footprint_tools/stats/distributions/t.pyx"]),
    dict(name="footprint_tools.stats.differential", sources=["footprint_tools/stats/differential.pyx"]),
]

install_requires = [
    "cython",
    "numpy>=1.10",
    "scipy>=0.17",
    "pandas",
    "pysam>=0.15",
    "statsmodels",
    "genome-tools>=1.0.3",
    "pwlf",
    "simplejson",
    "tqdm",
    "click",
    "click-option-group",
    "matplotlib",
    "yaspin"
]

__version__ = "1.3.5"

setup(
    name = "footprint_tools",
    version = __version__,
    license = "GPL-3.0-or-later",
    description = "Genomic footprint detection",
    long_description=(pathlib.Path(__file__).parent / "README.md").read_text(),
    long_description_content_type="text/markdown",
    author = "Jeff Vierstra",
    author_email = "jvierstra@altius.org",
    url = "https://github.com/jvierstra/footprint-tools",
    download_url = "https://github.com/jvierstra/footprint-tools/archive/v{}.tar.gz".format(__version__),
    keywords = ["genomic footprints", "bioinformatics", "chromatin", "transcription factors", "dnase"],
    zip_safe = False,
    packages =  find_packages(),
    libraries = [hcehpes_lib],
    ext_modules = [Extension(**opts) for opts in modules],
    include_dirs=[np.get_include(), hcephes_include],
    package_data={'footprint_tools':['*.pxd', 'logging.conf']},
    cmdclass = {'build_clib': build_clib, 'build_ext': build_ext},
    install_requires = install_requires,
    include_package_data = True,
    #scripts = scripts,
    entry_points = {"console_scripts": ["ftd = footprint_tools.__main__:main"]},
    classifiers=[
        'Development Status :: 5 - Production/Stable', 
        'Intended Audience :: Science/Research', 
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],
)
