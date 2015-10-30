footprint-tools
===============
A simple python package to detect footprints by simulating
expected cleavage rates using a 6-mer DNase I cleavage preference
model combined with tag smoothing. Statistical significance of
per-nucleotide cleavages are computed from an emperically fit 
negative binomial distribution.

Requires:
---------
+ pysam
+ pyfaidx
+ numpy
+ scipy
