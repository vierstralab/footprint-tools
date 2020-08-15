# footprint-tools: *de novo* genomic footprint detection 

footprint-tools is a python module for *de novo* detection of genomic footprints from DNase I data.

Footprint-tools detects footprints by simulating expected
cleavage rates using a 6-mer DNase I cleavage preference model combined
with density smoothing. Statistical significance of per-nucleotide cleavages
are computed from a series emperically fit negative binomial distribution.

## Requirements

``footprint-tools`` requires Python (>=3.5) and depends on the following additional
packages:

* cython
* numpy
* scipy
* statsmodels
* pysam
* genome_tools (http://www.github.com/jvierstra/genome_tools)
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

[Vierstra2020](https://doi.org/10.1038/s41586-020-2528-x) Vierstra, J., Lazar, J., Sandstrom, R. *et al.* Global reference mapping of human transcription factor footprints. *Nature* **583**, 729–736 (2020)
