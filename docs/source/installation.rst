
Quickstart
==========

Requirements
------------

``footprint-tools`` requires Python (≥3.6) and depends on the following additional
packages:

- numpy>=1.10
- scipy>=0.17
- pysam>=0.15
- pandas
- statsmodels
- pwlf
- simplejson
- cython
- genome_tools>=1.0.2
- tqdm
- click
- click-option-group
- matplotlib
- yaspin

We also recommend these non-python analysis tools:

- `samtools <http://www.htslib.org/>`_
- `BEDOPS <http://bedops.readthedocs.io>`_
- `kentUtils <https://github.com/ucscGenomeBrowser/kent/tree/master/src/utils>`_

Installation
------------

To install the latest stable release, type:

.. code::  bash

	pip install git+git://github.com/jvierstra/footprint-tools.git --user

Alternatively, you can get the latest build from GitHub:

.. code:: bash
	
	git clone https://github.com/jvierstra/footprint-tools.git
	cd footprint-tools
	pip install .

Finally, you can setup a ``conda`` enviroment using the provided YAML file:

.. code:: bash

	conda env create -f conda-env.yml

Basic usage
-----------

Footprint identification requries two steps: (1) learning a dispserion model 
and (2) apply dispersion model to data. This is steps are perform by distinct
commands:

.. code:: bash
	ftd learn_dm [intervals] [bam] [fasta] --bias_model_file [bm]
	ftd detect [intervals] [bam] [fasta] --bias_model_file [bm] --dispersion_model_file [dm]

```learn_dm``` returns a JSON-formatted dispersion model that is used for all 
downstream steps of footprint detect and analysis.