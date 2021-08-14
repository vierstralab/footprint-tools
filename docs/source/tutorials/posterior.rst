.. _tutorial_posterior:

Emperical Bayes footprint detection
===================================

Overview
~~~~~~~~

To call footprints jointly considering hundreds-to-thousands of datasets, we developed an Empirical Bayes approach that computes the posterior probability that a particular region of the genome is occupied by a motif.


.. math::
	p(\theta_+|X, Y) = \frac{P(\theta_+) \mathcal{L}(X,Y|\theta_+)}{	p(\theta_+) \mathcal{L}(X,Y|\theta_+) + (1-P(\theta_+)) \mathcal{L}(X,Y|\theta_-)}

where :math:`X` corresponds to the observed cleavage rates and :math:`Y` represents the expected cleavage rates. The footprint prior, :math:`P(\theta_+)`, is the number of datasets that a nucleotide is found within a footprint (FDR<0.05; user-defined) divided by the number datasets in which that nucleotide lays within a DNase I hypersensitive site (as defined by `hotspot2 <https://github.com/Altius/hotspot2>`_). The likelihood function corresponding an **unoccupied nucleotide**, :math:`\mathcal{L}(X,Y|\theta_-)`, is the product of individual negative binomial probabilities corresponding to the dispersion model parameterized by the expected cleavage rate. 

.. math::
	\mathcal{L}(X,Y|\theta_-) = \sum_{j=i-3}^{i+3} \left( \begin{array}{c} x_j + r_j -1 \\ x_j \end{array} \right) (1-p_j)^{r_j} {p_j}^{x_j}

.. math::
	:nowrap:

	\begin{gather*}
	\theta_- = \Phi_{NB}(y_j) = (\mu_j, r_j)\\
	p_j = \frac{r_j}{r_j + \mu_j}
	\end{gather*}

Here, :math:`\Phi_{NB}(y_j)` is a function that returns the parameters :math:`(\mu, r)` of the fitted negative binomial at the expected cleavage rate :math:`y_j` (see :ref:`tutorial_detect_dispersion-model`).


The likelihood function corresponding to an **occupied nucleotide**, :math:`\mathcal{L}(X,Y │θ_+)`, is determined similar to the unoccupied case after scaling the expected cleavage rate by the expected depletion of cleavage at occupied nucleotides:

.. math::

	\theta_+ =  \Phi_{NB}(y_j \lambda_j) = (\mu^{+}_j, r^{+}_j) 

where :math:`\lambda_j` is the expected depletion of cleavage at an occupied nucleotide. We determine :math:`\lambda_j` by considering all datasets with an FDR 5% footprint at position :math:`j`. First, for each dataset we fit a Beta distribution to the ratio of observed over expected cleavages (depletion ratio) at all FDR 5% footprints identified within individual datasets (capping the ratio values at 1.0). Then, for each nucleotide we re-estimate the depletion ratio by updating the Beta distribution (:math:`\alpha' = \alpha + x_i`, :math:`\beta’ = \beta + (y_i–x_i)`). These updated parameters are used to generate maximum a posteriori (MAP) estimates of the depletion ratio (:math:`\mu_{MAP}`) and expected variation of this ratio (:math:`\sigma^2_{MAP}`) at each nucleotide. A per-nucleotide footprint depletion estimate is then finally calculated from the average of the MAP mean estimates weighted by the inverse of the MAP standard deviation considering all datasets with an identified footprint at that nucleotide. 

Step-by-step guide
~~~~~~~~~~~~~~~~~~~

Step 1: Call footprints in individual datasets
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To build the priors and likelihood functions, footprints need to first be identified 
in each individual dataset. Please see  :ref:`command_detect` for this task.

Step 2: Compute the per-dataset footprint priors
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Used the command :ref:`command_learn_bm` to compute the hyperparameters of the Beta prior. 

.. code:: bash

	ftd learn_beta \
		--fdr-cutoff 0.05 --exp-cutoff 10 \
		interval.bedgraph

.. note:: 

	This command is run for every dataset before continuing to Step 2.

Step 2: Create a metadata file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Next, we create a metadata file that contains  the pertinent 
information for each dataset. The format of this file is **tab-delimited** and 
must contain the the following columns in a header (additional columns are
permissable). Lines are ignore when ``#`` is the first character.

+---+-------------------------+------------------------------------------+
| # | Column                  | Description                              |
+===+=========================+==========================================+
| 1 | ``id``                  |  Dataset identifier                      |
+---+-------------------------+------------------------------------------+
| 2 | ``dm_file``             |  Dispersion model filepath  (JSON file)  |
+---+-------------------------+------------------------------------------+
| 3 | ``tabix_file``          |  Output file from :ref:`command_detect`  |
|   |                         |  Note: must be gzipped with tabix index  |
+---+-------------------------+------------------------------------------+
| 4 |  ``beta_a``             |  Beta distribution parameters            |
+---+-------------------------+  Output from :ref:`command_learn_beta`   |
| 5 | ``beta_b``              |                                          |
+---+-------------------------+------------------------------------------+

Step 3: Compute posterior probabilites
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The posterior footprint probabilities are called using the command :ref:`command_posterior`. 
This command takes both the metadata file created above and a BED-formated file containing 
the genomic regions where footprint detection will occur. Typically, the input regions are 
defined by merging the DNase I hotspots across all samples.

.. code:: bash

	ftd posterior sample_data.txt intervals.bed

**Example output:**

This script writes a bedgrah-like file. Each row consists of an individual nucleotide and 
columns correspond to datasets (in the same order as the input metdata file).


.. note::

	Because this is a potentially huge operation (millions of DHS vs. hundreds of samples), we 
	typicall split the input file (DHSs) into chunks and the parallel process the chunks on
	a high-performance computing cluster.

	.. code:: bash

		cat regions.bed | split -l 5000 -a 4 -d - regions.chunk.

		regions.chunk.0000
		regions.chunk.0001
		regions.chunk.0002
		...

Step 4: Retrieve footprints
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Footprints (per dataset) can be retrieved by thresholding on posterior probabilities

.. code:: bash

   cat per-nucleotide.posterior.bedgraph \
       | awk -v OFS="\t" -v col=45 -v thresh=0.01 \ # set column to dataset column
       		'$(col) >= -log(thresh) { print $1, $2-3, $3+3; }' \
       | sort-bed --max-mem 8G - \
       | bedops -m - \
   > footprints.bed




