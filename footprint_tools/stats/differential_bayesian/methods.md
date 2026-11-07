# differential_bayesian

Per DHS: footprints (positions where at least one group is depleted) and differential footprints (positions where a
set of groups is more or less depleted than the other groups). Depths are log2 observed over expected cleavage.

| Step | Module | Output |
|---|---|---|
| Stage 1: each sample segmented on its own, no labels | `segmentation.py` | posterior mean $m_{sj}$ and variance $v_{sj}$, samples × positions |
| Stage 2: W of every set, jackknife SD, permutation p-values | `pvalues.py` | per position: best set per direction, log10 p |
| Footprints, calls, multiple testing | `calls.py` | runs, BH, Storey's $\pi_0$ |
| Container, calls at a cutoff | `results.py` | `DifferentialFootprintFit` |
| DataBundle loader, NPZ storage | `api.py` | `DifferentialFootprintLoader`, `save_data_results`, `load_fit` |
| Genome-wide pooling | `report.py` | calls, DHS and cutoff tables, summary page |

**Stage 1.** $y_{sj} \sim \mathrm{NB}(m_s(\lambda_{sj}2^{\theta_{sj}}), r_s(\lambda_{sj}2^{\theta_{sj}}))$ with each
sample's dispersion model. $\theta_s$ is piecewise constant: block lengths ~ discretized gamma (mean 24.8, shape 1.9,
2–150 bp), block levels ~ $\mathcal N(0, 3^2)$ on a grid from −9 to 4 in steps of 0.2, integrated out per block. Exact
forward–backward over all segmentations gives $m_{sj}$ and $v_{sj}$.

**Effect.** For samples $s, t$: $D = m_t - m_s \sim \mathcal N(\mu, \sigma^2)$, $\sigma^2 = v_s + v_t$,
$E[D_+] = \mu\Phi(\mu/\sigma) + \sigma\phi(\mu/\sigma)$, $E|D| = 2E[D_+] - \mu$.
$W_{gh} = \sum_{s\in g,t\in h} E[D_+] / \sum E|D|$ (group $g$ more depleted than $h$), so $W_{gh} + W_{hg} = 1$.
For a set $S$ of $k$ groups, $W_S = \mathrm{mean}_{g\in S, h\notin S} W_{gh} = \big(\sum_{g\in S} R_g - \binom{k}{2}\big)/(k(G-k))$
with $R_g = \sum_{h\ne g} W_{gh}$; $W_S > \tfrac12$: $S$ more depleted, $< \tfrac12$: less.

**Test.** $\mathrm{SD}_S$: grouped delete-one jackknife, $\mathrm{SD}_S^2 = \sum_g \frac{n_g-1}{n_g}\sum_{s\in g}(\Delta_s - \bar\Delta_g)^2$,
$\Delta_s$ the change in $W_S$ without sample $s$ (exact, from row sums). Per position and threshold $w$,
$T = \max_{S,\pm} (|W_S - \tfrac12| - (w - \tfrac12))/\mathrm{SD}_S$ over all sets of up to 5 groups. Reference: $T$ with
$w = \tfrac12$ under 1000 label permutations (studentized, one permutation for all positions);
$p = (1 + \#\{T^{(\pi)} \ge T\})/(B+1)$, or a generalised Pareto tail fitted to the top 10% of the permutation maxima
when fewer than 10 reach $T$. Exact for "no difference"; for $w > \tfrac12$ approximate. Per DHS: the same with the best
window of 5 positions ($\max$ over windows of $\min T$).

**Calls.** Benjamini–Hochberg over all positions genome-wide at the call threshold $w_{\min}$ (default 0.8) gives one
cutoff. Per direction, the best set at each position; a call is ≥ 5 consecutive positions passing in one direction; a
group belongs to it where it is in the best set for ≥ 5 consecutive positions. A more- and a less-depleted call may
overlap only where their sets share no group. **Footprints**: ≥ 5 positions with
$P(\min_g \bar\theta_{gj} < -\log_2 1.5) \ge 0.95$, groups as independent normals.

Settings are in `config.py`.
