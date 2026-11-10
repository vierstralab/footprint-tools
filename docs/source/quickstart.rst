
Quickstart
==========

footprint-tools requires Python 3.10 or greater


.. panels::
    :card: + install-card
    :column: :column: col-12 p-3

    Using pip?
    ^^^^^^^^^^

	footprint-tools can be installed via pip from GitHub.

    ++++

	Install using pip::

		pip install git+https://github.com/vierstralab/footprint-tools.git
    
	---
    :column: col-12 p-3

    Source from GitHub
    ^^^^^^^^^^^^^^^^^^

    .. code-block:: bash

	git clone https://github.com/vierstralab/footprint-tools.git
	cd footprint-tools
	pip install .		


Dependencies
^^^^^^^^^^^^
- numpy
- scipy
- pandas
- pysam
- pwlf
- simplejson
- genome-tools (>= 3.0.0)
- tqdm
- click
- click-option-group
- matplotlib
- yaspin
- colorlog

We also highly recommend these non-Python analysis tools:

- `samtools <http://www.htslib.org/>`_
- `BEDOPS <http://bedops.readthedocs.io>`_
- `kentUtils <https://github.com/ucscGenomeBrowser/kent/tree/master/src/utils>`_
