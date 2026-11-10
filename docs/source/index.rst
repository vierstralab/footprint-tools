footprint-tools: digital genomic footprint detection and analysis
=================================================================

A Python package for *de novo* detection of genomic footprints from DNase I data.

**Features:**

- De novo footprint detection
- Consensus footprint detection (emperical Bayes)
- Differential footprinting (Python API only for now, ``footprint_tools.stats.differential_bayesian``; see `methods <https://github.com/vierstralab/footprint-tools/blob/HEAD/footprint_tools/stats/differential_bayesian/methods.md>`_ and an example in `regulome <https://github.com/vierstralab/regulome/blob/HEAD/Figure5.footprints/DiffTest.ipynb>`_)
- API for programmatic access to cleavage data directly from sequence alignment files

.. warning::

	Please note that is documentation is a work in progress. Please `contact us <mailto:jvierstra@altius.org>`_ with any questions.

Contents
--------

.. toctree::
        :maxdepth: 1
        
        quickstart.rst
        commands.rst
        tutorials.rst
        api.rst

Pre-processed data
-------------------

As part of the `ENCODE Consortium <http://encodeproject.org>`_ we have generated digital genomic footprinting data for >240 cell types and tissues. These data and their metadata are hosted at both `<http://vierstra.org/resources/dgf>`_ and `ZENDODO <https://doi.org/10.5281/zenodo.3603548>`_.

Citation
--------

.. [Vierstra2020] Vierstra, J., Lazar, J., Sandstrom, R. *et al.* Global reference mapping of human transcription factor footprints. *Nature* **583**, 729–736 (2020) `<https://doi.org/10.1101/2020.01.31.927798>`_
