footprint-tools: *de novo* genomic footprint detection 
======================================================

footprint-tools is a python module for *de novo* detection of genomic footprints from DNase I data.

**Features:**

- De novo footprint detection
- Consensus footprint detection (emperical Bayes)
- Differential footprinting
- API for programmatic access to cleavage data directly from sequence alignment files


.. warning::

	Please note that is documentation is a work in progress. Please contact us with any questions.


Contents
--------

.. toctree::
        :maxdepth: 2
        
        detect.rst
        posterior.rst
        differential.rst
        api.rst
        examples.rst

Requirements
------------

``footprint-tools`` requires Python (>=3.5) and depends on the following additional
packages:

- cython (>=0.22)
- numpy (>=1.10)
- scipy (>=0.16)
- pysam 
- statsmodels
- genome_tools (>=1.0.1) (http://www.github.com/jvierstra/genome_tools)
- pyfaidx (https://github.com/mdshw5/pyfaidx)
- pwlf (https://github.com/cjekel/piecewise_linear_fit_py)


We also recommend these non-python analysis tools:

- `samtools <http://www.htslib.org/>`_
- `BEDOPS <http://bedops.readthedocs.io>`_: fast and simpler ``bedtools`` alternative
- `kentUtils <https://github.com/ucscGenomeBrowser/kent/tree/master/src/utils>`_

Installation
------------

To install the latest release, type:

.. code::  bash

	pip install footprint-tools

Citation
--------

.. [Vierstra2020] Global reference mapping and dynamics of human transcription factor footprints. Vierstra J *et al.* (2020)
   	    bioRxiv. 
   	    `<https://doi.org/10.1101/2020.01.31.927798>`_


