# footprint-tools: *de novo* genomic footprint detection 

footprint-tools is a python module for *de novo* detection of genomic footprints from DNase I data.

Footprint-tools detects footprints by simulating expected
cleavage rates using a 6-mer DNase I cleavage preference model combined
with density smoothing. Statistical significance of per-nucleotide cleavages
are computed from a series emperically fit negative binomial distribution.

## Requirements

``footprint-tools`` requires Python (>=3.5) and depends on the following additional
packages:

* cython (>=0.22)
* numpy (>=1.10)
* scipy (>=0.16)
* genome_tools (>=1.01) (http://www.github.com/jvierstra/genome_tools)
* pysam
* pyfaidx
* statsmodels
* pwlf


We also recomended these analysis tools (non-python):

* bedops (https://bedops.readthedocs.io/en/latest/)
* kentUtils (https://github.com/ucscGenomeBrowser/kent)


## Installation

To install the latest release, type:
```
	pip install footprint-tools
```

footprint-tools documentation is available [here](http://footprint-tools.readthedoc.org/en/latest)

Question and comments should be sent to jvierstra (at) altius (dot) org
