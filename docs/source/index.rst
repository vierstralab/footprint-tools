footprint-tools: *de novo* genomic footprint detection 
======================================================

footprint-tools is a python module for *de novo* detection of genomic footprints from DNase I data.

Footprint-tools detects footprints by simulating expected
cleavage rates using a 6-mer DNase I cleavage preference model combined
with density smoothing. Statistical significance of per-nucleotide cleavages
are computed from a series emperically fit negative binomial distribution.


Contents
--------

.. toctree::
        :maxdepth: 2
        
        detect.rst
        posterior.rst
        differential.rst
        api.rst

Requirements
------------

``footprint-tools`` requires Python (>=3.5) and depends on the following additional
packages:

-  cython (>=0.22)
-  numpy (>=1.10)
-  scipy (>=0.16)
-  genome_tools (http://www.github.com/jvierstra/genome_tools)
-  pysam
-  pyfaidx
-  statsmodels
-  pwlf


We also recomended these analysis tools (non-python):

-  bedops
-  kentUtils

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


