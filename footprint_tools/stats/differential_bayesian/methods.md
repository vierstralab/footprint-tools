# differential_bayesian

Per DHS: which bases are footprinted in any group, which groups are more or less depleted than
the others, and how much of the variation lies between groups (ICC). All quantities are log2.

## Pipeline

| Step | Module | Output |
|---|---|---|
| Per-sample contrast k_s and floor f_s, once over many DHSs | `samples.py` | `sample_maps.tsv` |
| Group and sample summaries from counts | `summaries.py` | `Summaries`: x, v, sigma2 (G × n); theta_mean, theta_sd (S × n) |
| Two independent segmentations, exact posterior draws | `model.py` | depth draws, boundary probabilities, ICC |
| Tracks and calls | `results.py` | `FootprintFit` |
| Loaders, NPZ storage | `api.py` | `SummaryLoader`, `FootprintLoader`, `save_data_results`, `load_data_results` |
| Global constants and group nuggets, once over all DHSs | `hyper.py` | `ModelConfig` JSON |

## Sample maps (`samples.py`)

A sample sees the common depletion θ as ρ_s(θ) = f_s + (1 − f_s) 2^(k_s θ): k_s stretches it, f_s is the
share of the expected cleavage that binding does not remove. Per DHS, each base gets a reference depth
log2((ΣO + ½)/(ΣE + ½)) over a 3-bp window and the other half of the samples (halves split within groups);
each sample's observed and expected cuts are summed in 0.25-log2 bins of that depth (−6 to 3), over bases where
it expects ≥ 4 cuts. Over DHSs rich in deep footprints (above the median share of expected cuts in bins < −2),
log2(O/E) per sample and bin is fitted to log2 ρ_s(d_b) with a latent depth per bin, by alternating weighted
least squares; median k_s = 1, f_s ∈ [0, 0.6], k_s ∈ [0.2, 3].

## Summaries

Counts follow a negative binomial with mean m_s(e) and size r_s(e) from each sample's dispersion model,
at e = ρ_s(θ) · expected. Below the first breakpoint the mean is extended through the origin and 1/r is held
at its first-breakpoint value, so the floor is f_s rather than the fitted intercept. Sample depletions
θ ~ N(μ_g, σ²), with an inverse-gamma prior on σ² fitted per group to per-base log-ratio variances; θ and σ²
are integrated out on grids. Per group and base: x, v = posterior mean and variance of μ_g under a flat prior;
sigma2 = posterior mean of σ². Per-sample θ is summarized under N(x_g, σ² scale + v_g), for display only.
Bases with zero expected cleavage carry no information.

## Model

    x_gj = b'_B + δ_gA + ε_j + e_gj,   ε_j ~ N(0, τ²),   e_gj ~ N(0, v_gj + t_g)
    b'_B ~ N(0, β),   δ_A ~ N(0, ω² K) with Σ_g w_g δ_g = 0,   ω ~ half-Cauchy(s_ω) on a grid

B is a consensus block, A a group-pattern block; the two partitions are independent. t_g is the group
nugget: variance of group g's summaries beyond v that the shared τ² does not cover. With
w_g ∝ 1 / mean_j(v_gj) and P = I − 1wᵀ, the weighted mean m_j = Σ w_g x_gj depends only on b′ and
the contrasts P x_j only on δ, so the evidence factorizes: Z = Z_consensus · Z_pattern. Block
evidences are closed-form from prefix sums; each layer's sum over partitions is one forward
recursion, and draws are exact backward samples. Block lengths: discretized gammas, one per layer.

## Outputs (`FootprintFit`, per base)

- `depth_mean`, `depth_lo`, `depth_hi`: posterior mean and 5–95% interval of b′ + δ_g.
- `p_any`: P(min_g depth < 0). Consensus footprints: runs with `p_any` ≥ `prob` (0.95).
- `delta_mean`, `p_more`, `p_less`: groups are compared on u = log2(c + (1 − c) 2^depth), c = `call_floor` (0.1),
  the depth as a sample with a 10% floor would see it, so differences deep below what counts can resolve do not
  become calls. Δ_g = u_g − median of the other groups' u; P(Δ_g < −delta), P(Δ_g > delta), delta = log2 1.5.
  Calls: runs with probability ≥ `prob` and at least `min_len` (6) bases.
- `p_pattern`, `p_consensus`: boundary probabilities of each partition between j and j + 1.
- `icc`: E[ω² / (ω² + σ̄²_A)], σ̄²_A = Σ_g w_g sigma2_gj averaged over the group-pattern block A.

## Constants (`ModelConfig`, `hyper.py`)

τ² and β maximize the summed consensus evidence, s_ω the summed group-pattern evidence (per-DHS tables, one fit
over all DHSs); defaults τ² = 0.6, β = 1, s_ω = 0.2. Block lengths are resolution choices: consensus
gamma(mean 24.8, shape 1.9, min 2), group pattern gamma(mean 20, shape 3, min 2). The nuggets t_g are fitted by
moments from base-to-base differences of the contrasts c_gj = d_gj − Σ_h w_h d_hj (d = x_j+1 − x_j), where block
levels and τ² cancel: Var c_g = (1 − 2w_g)(V_g + 2t_g) + Σ_h w_h² (V_h + 2t_h), V = v_j + v_j+1; differences beyond
5 SD (block edges) are dropped, and the weighted equations are summed over DHSs.

## Running

1. `-entry sample_maps`: per-DHS profiles of a random subset of DHSs → `sample_maps.tsv`.
2. Main run with `--sample_maps`: results, and per-DHS statistics (`stats.<dhs>.npz`: evidence tables, nugget
   equations, group names) → `model_config.json`.
3. Rerun with `--sample_maps` and `--model_config`.

## Tests

`tests/` check each layer's recursion and sampled boundary frequencies against enumeration of all
partitions, one block's evidence against a dense Gaussian density, NPZ round trips with
realignment to a sub-interval, the NB likelihood against scipy, and recovery of planted sample maps and nuggets.
