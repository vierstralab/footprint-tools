# footprint-tools — notes for coding agents

De novo detection of DNase I footprints: expected cleavage per base from a 6-mer bias model
and density smoothing, negative binomial statistics of observed vs. expected cleavage, Bayesian
posteriors across many samples, and differential footprints between groups of samples.
User docs: `README.md`, https://footprint-tools.readthedocs.io (`docs/source/`).

Downstream code imports it, so keep public names and signatures stable:
- genome-tools: footprint plot loaders (`loaders/footprint.py`: `stats.windowing.stouffers_z`,
  `stats.posterior`, `stats.differential`, `modeling.dispersion`) and the reads loader
  (`loaders/variant.py`: `cutcounts.BamFileExtractor`)
- nf-footprinting: the `ftd` CLI (`main.nf`, `consensus_footprints.nf`) and
  `stats.differential_bayesian.api` (`bin/extract_fp_data.py`)
- regulome Figure 5 notebooks (`differential_bayesian.api.load_fit`, `load_data_results`)

## When to use it

| Task | Use |
|---|---|
| Footprints in one sample from a BAM/CRAM | `ftd learn_dm` → `ftd detect` (CLI, below) |
| Footprints in many samples, on a cluster | the nf-footprinting pipeline (https://github.com/vierstralab/nf-footprinting), which runs the same `ftd` commands per sample |
| Footprints called jointly across samples | `ftd learn_beta` per sample, then `ftd posterior` (or nf-footprinting `consensus_footprints.nf`) |
| Differential footprints between groups of samples | Python only: `stats.differential_bayesian` (nf-footprinting `differential_footprints.nf` runs it per DHS) |
| Observed/expected cleavage of one region in Python (plots, custom statistics) | `cutcounts.BamFileExtractor` + `modeling.predict.ExpectedPredictor` |
| Plotting footprints and posteriors in multi-track figures | genome-tools `modular_plot` footprint components, not this package |

Not for: peak/hotspot calling (hotspot3), ATAC-seq (the bias model and offsets are for DNase I).

## Where things are

| Need | Look in |
|---|---|
| Cleavage counts from BAM/CRAM per strand (`BamFileExtractor.lookup`; `BamFile`, `bamfile` are older names) | `footprint_tools/cutcounts.py` |
| Sequence bias models (`KmerModel`, `UniformModel`) | `footprint_tools/modeling/bias.py` |
| Expected cleavage (`ExpectedPredictor.compute`) | `footprint_tools/modeling/predict.pyx` |
| Dispersion model (`DispersionModel`, `learn_dispersion_model`, `load_dispersion_model`, `write_dispersion_model`) | `footprint_tools/modeling/dispersion.pyx` |
| Windowed p-values (`stouffers_z`, `fishers_combined`, ...) | `footprint_tools/stats/windowing.pyx` |
| Empirical FDR, q-values | `footprint_tools/stats/fdr/__init__.py` |
| Posterior footprint probabilities across samples | `footprint_tools/stats/posterior.py` |
| Differential footprints (method write-up: `methods.md`) | `footprint_tools/stats/differential_bayesian/` |
| CLI (`ftd`) | `footprint_tools/__main__.py`, `footprint_tools/cli/*.py` |
| Diagnostic plots of a dispersion model | `footprint_tools/plotting/__init__.py` |
| Default 6-mer bias model | `data/vierstra_et_al.6mer-model.txt` |
| Tutorials (notebooks on public data) | `docs/source/tutorials/` |

`stats/differential.pyx` is the older differential test, still imported by genome-tools' footprint
loaders. `cli/learn_bm.py`, `cli/diff.py` and `cli/autocorr.py` are unfinished and not registered in `ftd`.

## How to use it

### CLI, one sample

```bash
ftd learn_dm --bias_model_file data/vierstra_et_al.6mer-model.txt --outfile dm.json \
    hotspots.bed reads.bam genome.fa                  # dispersion model; use all hotspots of the sample
ftd plot_dm dm.json --outfile dm.pdf                  # check the fit
ftd detect --bias_model_file data/vierstra_et_al.6mer-model.txt --dispersion_model_file dm.json \
    --outprefix out hotspots.bed reads.bam genome.fa  # out.bedgraph + out.fdr{0.001,0.01,0.05}.bed
ftd learn_beta --outfile beta.txt out.bedgraph        # Beta prior, only needed for `ftd posterior`
```

- BAM or CRAM is detected from the file (`--is_cram` is no longer needed). For CRAM, the FASTA
  must be the reference the CRAM was written with (htslib checks each contig's MD5).
- `out.bedgraph` columns: `chrom start end exp obs -log(pval) -log(winpval) fdr` (only `exp obs`
  without `--dispersion_model_file`). Footprints are runs of positions with `fdr` below the threshold.
- `--dedup` removes reads flagged as duplicates (default: keep them).

### CLI, joint posteriors across samples

`ftd posterior samples.tsv regions.bed --outprefix out`. `samples.tsv` has columns `id`, `dm_file`,
`tabix_file` (the sample's `out.bedgraph`, sorted, bgzipped, tabix-indexed with `tabix -0 -p bed`),
`beta_a`, `beta_b` (from `beta.txt`). Output: per-base `-log(1 - posterior)` for each sample.
`--post_cutoff 0.2` writes only positions where `1 - posterior < 0.2` in at least one sample.

### Python, one region

```python
from genome_tools import GenomicInterval
from genome_tools.data.extractors import FastaExtractor
from footprint_tools import cutcounts
from footprint_tools.modeling import bias, predict, dispersion
from footprint_tools.stats import windowing

reads = cutcounts.BamFileExtractor(bam_file)  # fasta_reference_filepath=... for CRAM
fasta = FastaExtractor(fasta_file)                # pysam.FastaFile also works
predictor = predict.ExpectedPredictor(reads, fasta, bias.KmerModel(bias_model_file),
                                      half_win_width=5, smoothing_half_win_width=50, smoothing_clip=0.01)
obs, exp, win = predictor.compute(GenomicInterval('chr19', 48363826, 48364602))
obs = obs['+'][1:] + obs['-'][:-1]           # dicts per strand; combine as the CLI does
exp = exp['+'][1:] + exp['-'][:-1]
dm = dispersion.load_dispersion_model('dm.json')
pvals = dm.p_values(exp, obs)
winpvals = windowing.stouffers_z(pvals, 3)   # 3 bp half window, as in `ftd detect`
```

`docs/source/tutorials/single_dataset.ipynb` runs this end to end on public data.

### Python, differential footprints

`fit_differential_footprints(groups, obs, exp, disp_models, config=DEFAULT_DF_CONFIG)`
(`differential_bayesian`) takes sample → group labels, samples x positions observed/expected
counts and sample → dispersion model, and returns a `DifferentialFootprintFit`
(`.calls(cutoff)`, `.region_table(cutoff)`). Every group needs at least two samples. Inside the
genome-tools plotting framework use `api.DifferentialFootprintLoader`; save and read per-DHS results
with `api.save_data_results` / `api.load_fit`; pool them genome-wide with
`python -m footprint_tools.stats.differential_bayesian.report <npz files> --prefix diff`.
Settings: `config.py` (`d_min`, the call threshold on the effect `d`, default 0.5).

## Gotchas

- Compiled Cython extensions (`.pyx`: `predict`, `dispersion`, `windowing`, `differential`,
  `distributions`, `stats.utils`). `No module named footprint_tools.stats.windowing` means Python
  found a `footprint_tools` folder without the built `.so` files:
  - running from inside a checkout that was never built (`pip install -e .` builds the `.so` in place);
  - a leftover `site-packages/footprint_tools/` folder without `__init__.py` (a namespace package
    shadowing the editable install; `footprint_tools.__file__` is `None`). Delete it and reinstall;
  - a kernel on another Python version than the one the `.so` was built for.
  After pulling changes to `.pyx` files, rerun `pip install -e .`.
- `GenomicInterval` has no strand unless one is given; `BamFileExtractor.lookup` treats a missing strand as `+`.
  Fragments from `lookup` carry `is_reverse`.
- `learn_dm` needs coverage, not many regions: the dispersion fit needs expected cleavages per base
  well beyond 25 (with defaults, at least 34 at >= 250 positions). Heavily duplicated libraries with
  `--dedup` can fall short; `learn_dispersion_model` then raises `DispersionModelFitError` and
  `ftd learn_dm` exits with a message. `r_breaks` must have exactly 4 values
  (`DispersionModel.fit_r` evaluates 5 segments).

## Working in this repo

- Package metadata and dependencies are in `pyproject.toml`; `setup.py` only builds the Cython
  extensions and the bundled `hcephes` C library. Bump the version in `pyproject.toml` (and
  `.bumpversion.cfg`).
- Build check: `pip install -e .` (builds in an isolated environment with Cython and numpy from
  `[build-system]`), then run `ftd --help` and each subcommand's `--help` from outside the checkout.
- No test suite. Quick check on real data: the tutorial notebook (`single_dataset.ipynb`), or
  `ftd learn_dm` + `ftd detect` on one chromosome of a sample's hotspots.
- Match the surrounding style: numpy-style docstrings, click options in `optgroup` groups, CLI
  datasets are `genome_tools.data.footprints.Dataset` subclasses with CamelCase names.
- Docs: the command pages are generated by sphinx-click from the click definitions; keep option
  help strings accurate. The README's "Differential footprinting" note and `methods.md` describe
  `differential_bayesian`; update them together with the code.
