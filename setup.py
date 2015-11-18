from distutils.core import setup
from distutils.command.build_clib import build_clib
from distutils.extension import Extension
from Cython.Distutils import build_ext

from glob import glob

import numpy as np

cephes_include = "cephes"
cephes_src = glob("cephes/*.c")
cehpes_lib = ('cephes', { 'sources': cephes_src })
ext_modules = [
	Extension("footprint_tools.modeling.dispersion", 
		sources = ["footprint_tools/modeling/dispersion.pyx"]),
	Extension("footprint_tools.stats.windowing", 
		sources = ["footprint_tools/stats/windowing.pyx"]),
	Extension("footprint_tools.stats.distributions",
		 sources = ["footprint_tools/stats/distributions.pyx"]),
	Extension("footprint_tools.stats.bisect",
		 sources = ["footprint_tools/stats/bisect.pyx"])
]

setup(
	name = "ftd-tools",
	description = "",
	long_description = "",
	author = "Jeff Vierstra",
	author_email = "jeffrv@uw.edu",    
	packages = ["footprint_tools"],
	libraries = [cehpes_lib],
    ext_modules = ext_modules,
    include_dirs=[np.get_include(), cephes_include],
    cmdclass = {'build_clib': build_clib, 'build_ext': build_ext}
)
