Differential footprint testing
==============================

To identify differences in footprint occupancy at true nucleotide resolution, we have developed a full Bayesian DFO (differential footprinting occupancy) test. The full documentation for this feature is currently under-construction; please see the Methods section of `Vierstra *et al.*, 2020 <https://doi.org/10.1038/s41586-020-2528-x>`_ for details.


Step 1: Call footprints on each dataset to be included in test
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
See :ref:`denovo-footprint-detection` for details on how to do this.

Step 2: Create metadata file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. |br| raw:: html

    <br>

Next, we create a metadata file that contains the pertinent information for each dataset. The format of this file is **tab-delimited**.

=== =========================  ==========================================
#   Column                     Description
=== =========================  ==========================================
1   ``group``                  Group label of dataset (either 'A' or 'B')
1   ``id``                     Dataset identifier
2   ``dispersion_model_file``  Dispersion model filepath  (JSON file) 
3   ``tabix_file``             Output file from ``ftd-compute-deviation`` |br|
                               Note: must be gzipped with tabix index 
=== =========================  ==========================================

.. note::

	Datasets must be listed in group order (i.e., All 'A' datasets occur before 'B' datasets).


Step 3: Run ``ftd-diff-test`` to identity differentially occupied nucleotides
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: bash

	[jvierstra@test0 ~]$ ftd-diff-test -h
	usage: ftd-diff-test [-h] [--processors N] [--tmpdir] metadata_file interval_file

	Perform a cleavage rate differential test between two groups

	positional arguments:
	  metadata_file   Path to configuration file containing metadata for samples
	  interval_file   File path to BED file

	optional arguments:
	  -h, --help      show this help message and exit

	Other options:
	  --processors N  Number of processors to use. Note that value excludes the minimum 2 threads that are dedicated to data I/O. (default: all available processors)
	  --tmpdir        Temporary directory to use (default: uses unix 'mkdtemp')

The ``ftd-diff-test`` script writes to standard out. 

**Ouptput file format:**

=== ============ ===========
#   Column       Description
=== ============ ===========
1   ``contig``   Chromosome
2   ``start``    Position (0-based)
3   ``end``      Position+1
4   ``lrt``      Likelihood ratio test -log *p*-value
5   ``winlp``    –log combined *p*-value (Stouffer's Z-score)
6   ``lf2c``     Estimated log2 difference between groups
=== ============ ===========
