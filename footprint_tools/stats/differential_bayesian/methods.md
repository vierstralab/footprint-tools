# differential_bayesian

Per DHS: which bases are footprinted in any group, which groups are more or less depleted than
the others, and how much of the variation lies between groups (ICC). All quantities are log2.

## Pipeline

| Step | Module | Output |
|---|---|---|
| Group and sample summaries from counts | `summaries.py` | `Summaries`: x, v, sigma2 (G × n); theta_mean, theta_sd (S × n) |
| Two independent segmentations, exact posterior draws | `model.py` | depth draws, boundary probabilities, ICC |
| Tracks and calls | `results.py` | `FootprintFit` |
| Loaders, NPZ storage | `api.py` | `SummaryLoader`, `FootprintLoader`, `save_data_results`, `load_data_results` |
| Global constants by summed evidence | `hyper.py` | `ModelConfig` JSON |
| Sequence correction of the expected | `seqcorr.py` | coefficients JSON |

## Summaries

Counts follow a negative binomial around the 6-mer expectation scaled by 2^θ, with each sample's
dispersion model. Sample depletions θ ~ N(μ_g, σ²), with an inverse-gamma prior on σ² fitted per
group to per-base log-ratio variances. θ and σ² are integrated out on grids. Per group and base:
x, v = posterior mean and variance of μ_g under a flat prior (the segmentation supplies the prior
on μ); sigma2 = posterior mean of σ². Per-sample θ is summarized under N(x_g, σ² scale + v_g), for
display only. Bases with zero expected cleavage carry no information.

## Model

    x_gj = b'_B + δ_gA + ε_j + e_gj,   ε_j ~ N(0, τ²),   e_gj ~ N(0, v_gj)
    b'_B ~ N(0, β),   δ_A ~ N(0, ω² K) with Σ_g w_g δ_g = 0,   ω ~ half-Cauchy(s_ω) on a grid

B is a consensus block, A a group-pattern block; the two partitions are independent. With
w_g ∝ 1 / mean_j(v_gj) and P = I − 1wᵀ, the weighted mean m_j = Σ w_g x_gj depends only on b′ and
the contrasts P x_j only on δ, so the evidence factorizes: Z = Z_consensus · Z_pattern. Block
evidences are closed-form from prefix sums; each layer's sum over partitions is one forward
recursion, and draws are exact backward samples. Block lengths: discretized gammas, one per layer.

## Outputs (`FootprintFit`, per base)

- `depth_mean`, `depth_lo`, `depth_hi`: posterior mean and 5–95% interval of b′ + δ_g.
- `p_any`: P(min_g depth < 0). Consensus footprints: runs with `p_any` ≥ `prob` (0.95).
- `delta_mean`, `p_more`, `p_less`: Δ_g = depth_g − median of the other groups; P(Δ_g < −delta),
  P(Δ_g > delta), delta = 0.5. Calls: runs with probability ≥ `prob`.
- `p_pattern`, `p_consensus`: boundary probabilities of each partition between j and j + 1.
- `icc`: E[ω² / (ω² + σ̄²_A)], σ̄²_A = Σ_g w_g sigma2_gj averaged over the group-pattern block A.

## Constants (`ModelConfig`)

τ² and β maximize the summed consensus evidence, s_ω the summed group-pattern evidence (`hyper.py`:
per-DHS evidence tables, then one fit over all DHSs). On 10 DHSs: τ² = 0.4, β = 1.0, s_ω = 0.2, the
same at every block-length setting tried. Block lengths are resolution choices: the evidence prefers
blocks of a few bp in both layers (geometric, means 8 and 5 bp), which fragments calls; defaults are
consensus gamma(mean 12, shape 3, min 2) and group pattern gamma(mean 20, shape 3, min 2).

## Sequence correction (`seqcorr.py`)

After the 6-mer correction a base-level residual remains: about 80% of it reproduces between random halves of the
samples and about 91% of that is common to all groups; it has no helical periodicity and does not concentrate in
protected stretches. It is predicted by local sequence: the residual r_j (11-bp high-pass of
log2((Σ obs + 0.5) / (Σ exp + 0.5)), sums over all samples) is regressed (ridge) on one-hot bases and dinucleotides at
offsets −6..6, and the expected becomes exp · 2^f(seq). On 10 DHSs the leave-one-chromosome-out R² is 0.34. Calls
barely change (the residual is shared, so it cancels in the contrasts); the consensus tracks become cleaner.

## Fitting the global parameters

Each pipeline run writes per-DHS statistics (`stats.<dhs>.npz`): sample sums of obs and of the uncorrected 6-mer exp,
the sequence of [start − 6, end + 7), and the evidence tables of the constants' grids. One collecting step fits
`seqcorr.json` and `model_config.json`; the next run uses them. The sequence correction is always fitted against the
uncorrected expected; the constants are fitted on the expected the run used.

## Tests

`tests/` check each layer's recursion and sampled boundary frequencies against enumeration of all
partitions, one block's evidence against a dense Gaussian density, NPZ round trips with
realignment to a sub-interval, and recovery of a planted sequence effect by `seqcorr`.
