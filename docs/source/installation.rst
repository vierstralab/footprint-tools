
Quickstart
==========

Requirements
------------

``footprint-tools`` requires Python (≥3.5) and depends on the following additional
packages:

- cython
- numpy
- scipy
- pysam 
- statsmodels
- genome_tools (http://www.github.com/jvierstra/genome_tools)
- pyfaidx (https://github.com/mdshw5/pyfaidx)
- pwlf (https://github.com/cjekel/piecewise_linear_fit_py)


We also recommend these non-python analysis tools:

- `samtools <http://www.htslib.org/>`_
- `BEDOPS <http://bedops.readthedocs.io>`_
- `kentUtils <https://github.com/ucscGenomeBrowser/kent/tree/master/src/utils>`_

Installation
------------

To install the latest stable release, type:

.. code::  bash

	pip install footprint-tools --user

Alternatively, you can get the latest build from GitHub:

.. code:: bash
	
	git clone https://github.com/jvierstra/footprint-tools.git
	cd footprint-tools
	python3 setup.py install --user
