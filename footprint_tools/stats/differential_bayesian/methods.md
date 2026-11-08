# differential_bayesian

Per DHS: footprints (positions where at least one group is depleted) and differential footprints (positions where a
set of groups is more or less depleted than the other groups). Depths are log2 observed over expected cleavage.

| Step | Module | Output |
|---|---|---|
| Stage 1: each sample segmented on its own, no labels | `segmentation.py` | posterior mean $m_{sj}$ and variance $v_{sj}$, samples × positions |
| Stage 2: effect d of every set, jackknife SD, permutation p-values | `pvalues.py` | per position: best set per direction, log10 p |
| Footprints, calls, multiple testing | `calls.py` | runs, BH, Storey's $\pi_0$ |
| Container, calls at a cutoff | `results.py` | `DifferentialFootprintFit` |
| DataBundle loader, NPZ storage | `api.py` | `DifferentialFootprintLoader`, `save_data_results`, `load_fit` |
| Genome-wide pooling | `report.py` | call, region (footprinted / not), DHS and cutoff tables, summary page |

**Stage 1.** $y_{sj} \sim \mathrm{NB}(m_s(\lambda_{sj}2^{\theta_{sj}}), r_s(\lambda_{sj}2^{\theta_{sj}}))$ with each
sample's dispersion model. $\theta_s$ is piecewise constant: block lengths ~ discretized gamma (mean 24.8, shape 1.9,
2–150 bp), block levels ~ $\mathcal N(0, 3^2)$ on a grid from −9 to 4 in steps of 0.2, integrated out per block. Exact
forward–backward over all segmentations gives $m_{sj}$ and $v_{sj}$.

**Effect.** For two samples $s, t$ the posterior depths differ by $D_{st} = \theta_t - \theta_s \sim \mathcal N(\mu_{st}, \sigma_{st}^2)$,
$\mu_{st} = m_t - m_s$, $\sigma_{st}^2 = v_s + v_t$: a difference $\mu_{st}$ and a distance
$\delta_{st} = E|D_{st}| = \mu_{st}(2\Phi(\mu_{st}/\sigma_{st}) - 1) + 2\sigma_{st}\phi(\mu_{st}/\sigma_{st})$.
For groups $g, h$, the ratio of their averages over all pairs $s\in g, t\in h$:
$$d_{gh} = \frac{\mathrm{mean}\,\mu_{st}}{\mathrm{mean}\,\delta_{st}} = \frac{\bar m_h - \bar m_g}{\mathrm{mean}\,\delta_{st}} \in [-1, 1],$$
the difference of group means over the mean absolute difference between the groups; $d_{gh} > 0$: $g$ more depleted,
$d_{hg} = -d_{gh}$. It is a magnitude-weighted Cliff's delta (with exact depths, the distance-weighted mean sign of the
pairwise differences); $|d_{gh}| = 1$ only when no pair goes against the means (Gini's transvariation is zero), and it
equals Dagum's relative economic affluence. For a set $S$ of $k$ groups, $d_S = \mathrm{mean}_{g\in S, h\notin S}\, d_{gh}
= \sum_{g\in S} R_g / (k(G-k))$ with $R_g = \sum_{h\ne g} d_{gh}$ (pairs inside $S$ cancel); $d_S > 0$: $S$ more depleted,
$d_{S^c} = -d_S$. Earlier versions used $W = (1 + d)/2$.

**Test.** $\mathrm{SD}_S$: grouped delete-one jackknife, $\mathrm{SD}_S^2 = \sum_g \frac{n_g-1}{n_g}\sum_{s\in g}(\Delta_s - \bar\Delta_g)^2$,
$\Delta_s$ the change in $d_S$ without sample $s$ (exact, from the row sums of $\delta$ and the group sums of $m$). Per
position and threshold $e$, $T = \max_{S,\pm} (|d_S| - e)/\mathrm{SD}_S$ over all sets of up to 5 groups (with both signs,
every split of the 11 groups). Reference: $T$ with $e = 0$ under 1000 label permutations (studentized, one permutation
for all positions); $p = (1 + \#\{T^{(\pi)} \ge T\})/(B+1)$, or a generalised Pareto tail fitted to the top 10% of the
permutation maxima when fewer than 10 reach $T$. Exact for "no difference"; for $e > 0$ approximate. Per DHS: the same
with the best window of 5 positions ($\max$ over windows of $\min T$).

**Calls.** Benjamini–Hochberg over all positions genome-wide at the call threshold $d_{\min}$ (default 0.6, i.e. $W \ge 0.8$) gives one
cutoff. Per direction, the best set at each position; a call is ≥ 5 consecutive positions passing in one direction; a
group belongs to it where it is in the best set for ≥ 5 consecutive positions. A more- and a less-depleted call may
overlap only where their sets share no group. **Footprints**: ≥ 5 positions with
$P(\min_g \bar\theta_{gj} < -\log_2 1.5) \ge 0.95$, groups as independent normals.

Settings are in `config.py`.
