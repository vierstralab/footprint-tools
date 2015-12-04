# FTD (footprint detector)
A simple python package to detect footprints by simulating
expected cleavage rates using a 6-mer DNase I cleavage preference
model combined with tag smoothing. Statistical significance of
per-nucleotide cleavages are computed from an emperically fit 
negative binomial distribution.

## Introduction

## Requirements

FTD requires Python (>=2.7) and depends on the following additional packages:

* numpy
* scipy
* pysam
* pyfaidx

## Usage

Detecting footprints with FTD is easy and requires the execution of two scripts.

### Step 1: Align sequenced DNase I cleavages

FTD requires an alignment file in BAM format which can be made using any sequence alignment tool. FTD uses all reads with a MAPQ > 0.

### Step 2: Create an index of the reference genome FASTA file

The software uses an indexed FASTA file to enable rapid lookups of genomic sequences utilized by the sequence bias model. A FASTA file can be indexed using `samtools`.

	[jvierstra@rotini footprint-tools]$ samtools faidx hg19.all.fasta

### Step 4: Download or create a 6-mer cleavage bias model

The sequence bias model is the basis of FTD. A model file contains 2 columns that contain a sequence k-mer and a relative preference value. While the bias model can be of any k-mer size, we typically use 6mers with the cleavage ocurring betwenn the 3 and 4 base.
	
	ACTTGC	0.22214673145056082
	ACTTAC	0.21531706520159422
	TCTTGC	0.20841281881342563
	ACTTGT	0.2031757090292924
	TCTCGC	0.19681882102542664
	TCTTAC	0.19519174298158992
	ACTCGC	0.1917931347689513
	ACTCGT	0.18406049938472563
	TCTTGT	0.18256420745577184
	TCTCGA	0.17989100314058748

### Step 5: Create a dispersion (error) model

FTD using a negative binomial to compute the significance of per-nucleotide cleavage devations from the expected. The negative binomial has two parameters, mu and r. The script `learn_dispersion_model.py` emperically fits mu and r from the observed cleavage data and then interpolates all values using linear regression. `learn_dispersion_model.py` writes a dispersion model in JSON format to standard out which can then be used with all FTD analyses.

	[jvierstra@rotini footprint-tools]$ python scripts/learn_dispersion_model.py -h
	usage: learn_dispersion_model.py [-h] [--kmer MODEL_FILE] [--half_win_width N]
	                                 [--processors N]
	                                 bam_file fasta_file interval_file

	Learn a negative binomial dispersion model

	positional arguments:
	  bam_file            File path to BAM-format tag sequence file
	  fasta_file          File path to genome FASTA file (requires associated
	                      FASTA index in same folder; see documentation on how to
	                      create an index)
	  interval_file       File path to BED file

	optional arguments:
	  -h, --help          show this help message and exit

	bias modeling options:
	  --bm MODEL_FILE	  Use a k-mer model for local bias (supplied by file). If
	                      argument is not provided the model defaults to uniform
	                      sequence bias.
	  --half_win_width N  Half window width to apply bias model. (default: 5)

	other options:
	  --processors N      Number of processors to use. (default: all available
	                      processors)

The dispersion model is typically generated from a random subset of the accessible regions wthin the genome.

### Step 6: Compute per-nucleotide expected cleavages

	[jvierstra@rotini footprint-tools]$ python scripts/compute_deviation.py -h
	usage: compute_deviation.py [-h] [--bm MODEL_FILE] [--half_win_width N]
	                            [--smooth_half_win_width N] [--smooth_clip N]
	                            [--dm MODEL_FILE] [--fdr_shuffle_n N]
	                            [--processors N]
	                            bam_file fasta_file interval_file

	Compute expected DNase I per-nucleotide cleavages

	positional arguments:
	  bam_file              File path to BAM-format tag sequence file
	  fasta_file            File path to genome FASTA file (requires associated
	                        FASTA index in same folder; see documentation on how
	                        to create an index)
	  interval_file         File path to BED file

	optional arguments:
	  -h, --help            show this help message and exit

	bias modeling options:
	  --bm MODEL_FILE       Use a k-mer model for local bias (supplied by file).
	                        If argument is not provided the model defaults to
	                        uniform sequence bias.
	  --half_win_width N    Half window width to apply bias model. (default: 5)

	smoothing options:
	  --smooth_half_win_width N
	                        Half window width to apply smoothing model. When set
	                        to zero no smoothing is applied. (default: 50)
	  --smooth_clip N       Fraction of signal to clip when computing trimmed
	                        mean. (default: 0.05)

	statistics options:
	  --dm MODEL_FILE		Dispersion model for negative binomial tests. If
	                        argument is not provided then no stastical output is
	                        provided. File is in JSON format and generated using
	                        the 'learn_dispersion_model' script included in the
	                        software package.
	  --fdr_shuffle_n N     Number of times to shuffle data for FDR calculation.
	                        (default: 25)

	other options:
	  --processors N        Number of processors to use. (default: all available
	                        processors)

The `compute_deviation.py` script writes to standard out. The ouptput format is quasi-bedGraph such that the columns contain information about (3) expected cleavages, (4) observed cleavages, (5) -log p-value of the per-nucleotide deviation from expected, (6) -log of the combined p-values using Stouffers Z-score method, and (7) the  calibrated FDR of column 6. 

	[jvierstra@rotini footprint-tools]$ python scripts/compute_deviation.py --bm vierstra_et_al.txt --dm model.json
		reads.bam genome.fa dhs.bed
	chr1    39585441        39585442        0       0       0.1719  0.0017  1.1579
	chr1    39585441        39585442        0       0       0.1719  0.0017  1.1579
	chr1    39585441        39585442        0       0       0.1719  0.0017  1.1579
	chr1    39585441        39585442        0       0       0.1719  0.0017  1.1579
	chr1    39585441        39585442        0       0       0.1719  0.0017  1.1579
	...

### Step 7: Retrieving footprints

#### 7.1: P-value thresholded

#### 7.2: FDR thresholded

## SGE parallelization
