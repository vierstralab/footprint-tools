# differential_bayesian

Per DHS: which bases are footprinted in any group, and which groups are more or less
depleted than the others, with posterior probabilities. All quantities are log2.

## Pipeline

| Step | Module | Output |
|---|---|---|
| Group and sample summaries from counts | `summaries.py` | `Summaries`: x, v (G × n); theta_mean, theta_sd (S × n) |
| Nested segmentation, exact posterior draws | `nested.py` | depth draws, boundary probabilities, ICC |
| Tracks and calls | `results.py` | `NestedFit` |
| Loaders, NPZ storage | `api.py` | `SummaryLoader`, `NestedLoader`, `save_data_results`, `load_data_results` |
| Global constants, once per release | `hyper.py` | `NestedConfig` JSON |

## Summaries

Counts follow a negative binomial around the 6-mer expectation scaled by 2^θ. Sample depletions
θ ~ N(μ_g, σ²), with an inverse-gamma prior on σ² fitted per group to per-base sample variances.
θ and σ² are integrated out on grids; x and v are the posterior mean and variance of μ_g under a
flat prior. Per-sample θ is summarized under the prior N(x_g, σ² scale + v_g), for display only.
Bases with zero expected cleavage carry no information.

## Nested model

Inside inner block B nested in outer block A:

    x_gj = b'_B + δ_gA + ε_j + e_gj,   ε_j ~ N(0, τ²),   e_gj ~ N(0, κ v_gj)
    b'_B ~ N(0, β),   δ_A ~ N(0, ω² K) with Σ_g w_g δ_g = 0,   ω ~ half-Cauchy(s_ω) on a grid

with w_g ∝ 1 / mean_j(κ v_gj), P = I − 1wᵀ, K = PPᵀ restricted to G − 1 groups. The weighted mean
m_j = Σ w_g x_gj depends only on b′ (inner blocks); the contrasts P x_j depend only on δ (outer
blocks). Block evidences are closed-form from prefix sums; one eigendecomposition per outer block
covers the ω grid. The sum over all nested segmentations is two forward recursions; the inner
prior is normalized within each outer block (q_out(L) / u(L)). Draws are exact backward samples.

Block lengths: discretized gamma (shape 3, 2–150 bp), means 20 bp (inner) and 70 bp (outer).
Both are resolution choices: the evidence prefers blocks of a few bp at both levels, and calls
change little between 20 and 70 bp for the outer level.

## Outputs (`NestedFit`, per base)

- `depth_mean`, `depth_lo`, `depth_hi`: posterior mean and 5–95% interval of b′ + δ_g.
- `p_any`: P(min_g depth < 0). Consensus footprints: runs with `p_any` ≥ `prob` (0.95).
- `delta_mean`, `p_more`, `p_less`: Δ_g = depth_g − median of other groups;
  P(Δ_g < −delta), P(Δ_g > delta), delta = 0.5. Calls: runs with probability ≥ `prob`.
- `p_outer`, `p_inner`: P(outer boundary), P(inner-only boundary) between j and j + 1.
- `icc`: E[ω² / (ω² + σ̄²)], σ̄² = 0.14.

## Constants (`NestedConfig`)

τ² = 0.4, s_ω = 0.15, β = 0.5 by summed log evidence over DHSs; κ = 2.2 by split-half
calibration (two halves of each group's samples: SD of the standardized disagreement in Δ = 1).
κ corrects block-level overconfidence from residuals correlated over neighbouring bases.

## Tests

`tests/` check the recursions against enumeration of all nested segmentations, one block's
evidence against a dense Gaussian density, sampler boundary frequencies against exact marginals,
and NPZ round trips with realignment to a sub-interval.
