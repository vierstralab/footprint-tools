# FTD (footprint detector)
A simple python package to detect footprints by simulating
expected cleavage rates using a 6-mer DNase I cleavage preference
model combined with tag smoothing. Statistical significance of
per-nucleotide cleavages are computed from an emperically fit 
negative binomial distribution.

## Requirements

FTD requires Python (>=2.7) and depends on the following additional packages:

* numpy
* scipy
* mpmath
* pysam
* pyfaidx

## Usage

Detecting footprints with FTD is easy and requires the execution of two scripts.

### Step 1: Align sequenced DNase I cleavages

### Step 2: Create an index of the reference genome FASTA file

### Step 3: Call DNase I hotspots

### Step 4: Download or create a 6-mer cleavage bias model

### Step 5: Create a dispersion (error) model

FTD using a negative binomial to compute the significance of per-nucleotide cleavage devations from the expected. The negative binomial has two parameters, mu and r. The script `learn_dispersion_model.py` emperically fits mu and r from the observed cleavage data and then interpolates all values using linear regression.

	[jvierstra@rotini footprint-tools]$ python scripts/learn_dispersion_model.py -h
	usage: learn_dispersion_model.py [-h] (--kmer MODEL_FILE | --uniform)
    	                             [--half_win_width N] [--procs N]
        	                         [--chunksize N]
            	                     bam_file fasta_file interval_file

	Learn a dispersion model

	positional arguments:
	  bam_file            File path to BAM-format tag sequence file
	  fasta_file          File path to genome FASTA file (requires associated
	                      index in same folder)
	  interval_file       File path to BED file

	optional arguments:
	  -h, --help          show this help message and exit
	  --kmer MODEL_FILE   Use a k-mer model for local bias
	  --uniform           Use a uniform model for local bias
	  --half_win_width N  Half window width (nt) to apply bias model (default: 5)
	  --procs N           Number of available processors (default: 8)
	  --chunksize N       Number of intervals to process per chunk (default: 500)

### Step 6: Compute per-nucleotide expected cleavages

	[jvierstra@rotini footprint-tools]$ python scripts/compute_expected.py  -h
	usage: compute_expected.py [-h] (--kmer MODEL_FILE | --uniform)
	                           [--smooth_tmean] [--half_win_width N]
	                           [--disp_model MODEL_FILE]
	                           bam_file fasta_file interval_file

	Compute expected DNase I per-nucleotide cleavages

	positional arguments:
	  bam_file              File path to BAM-format tag sequence file
	  fasta_file            File path to genome FASTA file (requires associated
	                        index in same folder)
	  interval_file         File path to BED file

	optional arguments:
	  -h, --help            show this help message and exit
	  --kmer MODEL_FILE     Use a k-mer model for local bias
	  --uniform             Use a uniform model for local bias
	  --smooth_tmean        Trimmed mean smoothing of expected window counts
	  --half_win_width N    Half window width (nt) to apply bias model (default:
	                        5)
	  --disp_model MODEL_FILE
	                        Compute p-values using a custom dispersion model (JSON
	                        format)

## SGE parallelization
