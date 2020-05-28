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
* statsmodels
* pysam
* genome_tools (1.0.1) (http://www.github.com/jvierstra/genome_tools)
* pyfaidx (https://github.com/mdshw5/pyfaidx)
* pwlf (https://github.com/cjekel/piecewise_linear_fit_py)

We also recommend these non-python analysis tools:

* [samtools](http://www.htslib.org/)
* [BEDOPS](http://bedops.readthedocs.io)
* [kentUtils](https://github.com/ucscGenomeBrowser/kent/tree/master/src/utils)

## Installation

To install the latest release, type:
```
pip install footprint-tools
```

## Documentation & usage

User manual, API and examples can be found [here](http://footprint-tools.readthedocs.io)

## Citation

[Vierstra2020](https://doi.org/10.1101/2020.01.31.927798) Global reference mapping and dynamics of human transcription factor footprints. Vierstra J *et al.* (2020) bioRxiv.
