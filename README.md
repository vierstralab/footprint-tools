# footprint-tools: *de novo* genomic footprint detection 

footprint-tools is a python module for *de novo* detection of genomic footprints from DNase I data 
by simulating expected cleavage rates using a 6-mer DNase I cleavage preference model combined
with density smoothing. Statistical significance of per-nucleotide cleavages
are computed from a series emperically fit negative binomial distribution.

## Requirements

``footprint-tools`` requires Python 3.6+

We also recommend these non-Python analysis tools:

* [samtools](http://www.htslib.org/)
* [BEDOPS](http://bedops.readthedocs.io)
* [kentUtils](https://github.com/ucscGenomeBrowser/kent/tree/master/src/utils)

## Installation

To install the latest release, type:
```
pip install footprint-tools
```

If you run into errors, try installing footprint-tools in a conda environment (using the YAML file provided):
```
# Clone repository
git clone https://github.com/jvierstra/footprint-tools.git

# Switch to devel branch; note that “master” branch is still “old” code
git checkout devel

# Create conda enviroment from config YAML file
cd footprint-tools
conda env create -f conda-env.yml

# Activate conda environment
conda activate footprint-tools

# Run commands
ftd --version
ftd {commands}
```

## Documentation & usage

User manual, API and examples can be found [here](http://footprint-tools.readthedocs.io)

## Citation

[Vierstra2020](https://doi.org/10.1038/s41586-020-2528-x) Vierstra, J., Lazar, J., Sandstrom, R. *et al.* Global reference mapping of human transcription factor footprints. *Nature* **583**, 729–736 (2020)
