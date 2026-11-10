# Copyright (C) 2015-2026 Jeff Vierstra (jvierstra@altius.org)

# Package metadata and dependencies are in pyproject.toml; this file only builds
# the Cython extensions and the bundled hcephes C library.

from glob import glob

from setuptools import Extension, setup
from setuptools.command.build_clib import build_clib
from Cython.Distutils import build_ext

import numpy as np

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
    dict(name="footprint_tools.stats.distributions.normal", sources=["footprint_tools/stats/distributions/normal.pyx"]),
    dict(name="footprint_tools.stats.differential", sources=["footprint_tools/stats/differential.pyx"]),
]

ext_modules = [Extension(**opts) for opts in modules]

for e in ext_modules:
    e.cython_directives = {'language_level': "3"} #all are Python-3)

setup(
    libraries = [hcehpes_lib],
    ext_modules = ext_modules,
    include_dirs=[np.get_include(), hcephes_include],
    cmdclass = {'build_clib': build_clib, 'build_ext': build_ext},
)
