
Quickstart
==========

footprint-tools requires Python is compatible with Python 3.6 or greater


.. panels::
    :card: + install-card
    :column: :column: col-12 p-3

    Using pip?
    ^^^^^^^^^^

	footprint-tools can be installed via pip from `PyPI <https://pypi.org/project/footprint-tools>`__.

    ++++

	Install using pip::
	
		pip install footprint-tools


	or for the lastest release::

		pip install git+https://github.com/jvierstra/footprint-tools.git
    
	---
    :column: col-12 p-3

    Source from GitHub
    ^^^^^^^^^^^^^^^^^^

    .. code-block:: bash

	git clone https://github.com/jvierstra/footprint-tools.git
	cd footprint-tools
	pip install .		


Dependencies
^^^^^^^^^^^^
- numpy
- scipy
- pandas
- pysam
- statsmodels
- pwlf
- simplejson
- cython
- genome_tools
- tqdm
- click
- click-option-group
- matplotlib
- yaspin

We also highly recommend these non-Python analysis tools:

- `samtools <http://www.htslib.org/>`_
- `BEDOPS <http://bedops.readthedocs.io>`_
- `kentUtils <https://github.com/ucscGenomeBrowser/kent/tree/master/src/utils>`_
