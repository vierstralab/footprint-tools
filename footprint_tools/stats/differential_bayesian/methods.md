# differential_bayesian

Per DHS: footprints (bases where at least one group is depleted) and differential footprints (bases where one
group's samples are reliably more or less depleted than the samples of the other groups). All depths are log2
observed over expected cleavage.

| Step | Module | Output |
|---|---|---|
| Stage 1: each sample segmented on its own, no labels | `segmentation.py` | posterior mean and variance of the depth, samples × positions |
| Stage 2: group tracks, discriminability, footprint probability | `calls.py` | groups × positions |
| Container, calls, region table | `results.py` | `DifferentialFootprintFit` |
| DataBundle loader, NPZ storage | `api.py` | `DifferentialFootprintLoader`, `save_data_results`, `load_data_results` |

**Stage 1.** $y_{sj} \sim \mathrm{NB}(m_s(\lambda_{sj}2^{\theta_{sj}}), r_s(\lambda_{sj}2^{\theta_{sj}}))$ with each
sample's dispersion model (below the first breakpoint the mean is extended through the origin and $1/r$ is held at
its first-breakpoint value). $\theta_s$ is piecewise constant: block lengths ~ discretized gamma (mean 24.8, shape
1.9, 2–150 bp), block levels ~ $\mathcal N(0, 3^2)$ on a grid from −9 to 4 in steps of 0.2. Exact forward–backward
sums over all segmentations give $m_{sj}$, $v_{sj}$.

**Stage 2.** Group mean $\bar\theta_{gj}$ and standard error; effect $\Delta_{gj} = \bar\theta_{gj} - \mathrm{median}_{h\ne g}\bar\theta_{hj}$.
Pairwise separation $\pi_{gh,j} = \mathrm{mean}_{s\in g,t\in h}\,\Phi((m_{tj}-m_{sj})/\sqrt{v_{sj}+v_{tj}})$, the
probability that a sample of $g$ is more depleted than a sample of $h$; discriminability
$A_{gj} = \mathrm{mean}_{h\ne g}\,\pi_{gh,j}$ (more depleted) or of $1-\pi_{gh,j}$ (less), with a 95% lower bound from
resampling samples within groups.

**Calls.** Group $g$, ≥ 5 bp in a row: lower bound of $A$ ≥ 0.75 and $|\Delta| > \log_2 1.5$; ranked by the lower
bound. **Footprints**: ≥ 5 bp with $P(\min_g \bar\theta_{gj} < -\log_2 1.5) \ge 0.95$, groups as independent normals.

Nothing in stage 2 is fitted globally; settings are in `config.py`.
