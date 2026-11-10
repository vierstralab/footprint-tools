# footprint-tools: *de novo* genomic footprint detection 

footprint-tools is a python module for *de novo* detection of genomic footprints from DNase I data 
by simulating expected cleavage rates using a 6-mer DNase I cleavage preference model combined
with density smoothing. Statistical significance of per-nucleotide cleavages
are computed from a series of emperically fit negative binomial distribution.


> [!CAUTION]
> There is a massive bug in the `posterior` footprint caller in versions 1.2.0 to 1.3.7. Please pull the lastest code from the master branch. 


Coding agents: read [AGENTS.md](AGENTS.md) first.

## Requirements

``footprint-tools`` requires Python 3.10+

We also recommend these non-Python analysis tools:

* [samtools](http://www.htslib.org/)
* [BEDOPS](http://bedops.readthedocs.io)
* [kentUtils](https://github.com/ucscGenomeBrowser/kent/tree/master/src/utils)

## Installation

To install the latest version from GitHub, type:
```
pip install git+https://github.com/vierstralab/footprint-tools.git
```

Or install footprint-tools in a conda environment (using the YAML file provided):
```
# Clone repository
git clone https://github.com/vierstralab/footprint-tools.git

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

To run footprint detection on many datasets, use the Nextflow pipeline [nf-footprinting](https://github.com/vierstralab/nf-footprinting).

Differential footprinting is available only from Python for now (`footprint_tools.stats.differential_bayesian`, see [methods](footprint_tools/stats/differential_bayesian/methods.md)); for an example, see [`DiffTest.ipynb`](https://github.com/vierstralab/regulome/blob/HEAD/Figure5.footprints/DiffTest.ipynb) in regulome.

## Citation

[Vierstra2020](https://doi.org/10.1038/s41586-020-2528-x) Vierstra, J., Lazar, J., Sandstrom, R. *et al.* Global reference mapping of human transcription factor footprints. *Nature* **583**, 729–736 (2020)
