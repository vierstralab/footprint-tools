Auto-generated

# Differential-footprint model documentation

This document describes the statistical model implemented in `differential_bayesian` v0.15.0.

```text
posterior.py        Generic finite-grid posterior utilities.
integration.py      Shared numerical integration helpers.
segmentation.py     Generic explicit-duration segmentation model.
differential.py     Group-level likelihoods over (mu, sigma2).
variance_ratio.py   Common-mean / variance-ratio hierarchy over (mu0, eta).
eta.py              ICC/eta segmentation and reconstructed mu0 posterior.
coefficients.py     Standardized group coefficients and coefficient-derived counts.
theta.py            Optional sample-level theta likelihoods and segmentation.
config.py           Public configuration dataclasses.
api.py              PlotDataLoader-based user-facing workflow.
io.py               Serializable mixin for to_dict/from_dict and to_npz/from_npz.
```

## Notation

Samples are indexed by $i$, groups or cell types by $g$, and genomic positions by $t$.

The main latent variables are:

- $\theta_{i,g,t}$: optional sample-level latent signal.
- $\mu_{g,t}$: group-level latent mean signal.
- $\sigma^2_{g,t}$: within-group variance.
- $\mu_{0,t}$: shared/common mean across groups.
- $\eta_t=\log\lambda_t$: log between-group variance ratio.
- $\mathrm{ICC}_t=\lambda_t/(1+\lambda_t)$: ICC-like transform of the variance ratio.
- $z_{g,t}$: standardized group coefficient, with the exact interpretation depending on the coefficient model.

All main continuous latent variables are represented on finite grids. A posterior mass on a grid point should be read as mass assigned to the bin around that grid point, not as a density value at an infinitesimal point.

## Public object map

The usual workflow attaches the following objects to a `DataBundle`-like `data` object:

```python
# group-level model
data.differential

# group mean segmentation
data.segmentation

# common mean / variance-ratio hierarchy
data.variance_ratio
data.eta_segmentation
data.mu0_segmentation

# group coefficients
data.common_coefficient_likelihood
data.common_coefficient_segmentation
data.zero_coefficient_likelihood
data.zero_coefficient_segmentation

# coefficient-derived group counts
data.kfp_zero
data.kfp_dev

# optional sample-level model
data.theta
data.theta_segmentation
```

The objects intentionally distinguish likelihoods from segmentations. Likelihood objects own their grids, priors, and pointwise posterior calculations. Segmentation objects own spatial posterior summaries and, for freshly fitted results, traceback caches for posterior path sampling.

---

# 1. Finite-grid posterior convention

Most posterior-returning methods produce a `GridPosterior`.

```python
posterior.x                # state grid
posterior.log_mass         # normalized log posterior mass, state on the last axis
posterior.mass
posterior.mean
posterior.map
posterior.quantile(q)
posterior.prob_gt(c)
posterior.prob_lt(c)
posterior.prob_abs_gt(c)
```

If $x_1,\ldots,x_M$ is a grid, `log_mass[..., j]` stores

$$
\log P(X=x_j\mid y)
$$

under the package's discretized model. Continuous priors such as normal priors are discretized by integrating the continuous density over bin intervals. If $b_j$ and $b_{j+1}$ are neighboring bin boundaries, then

$$
P(X=x_j)=\int_{b_j}^{b_{j+1}} p(x)\,dx.
$$

The helper `normal_grid_log_mass(x, mean, sd)` constructs these bin masses using normal CDF differences. The first and last bins are open-ended.

A key convention is that likelihood arrays and priors are stored separately. For example, `loglik_mu` stores the evidence over $\mu$ without the $\mu$ prior, while `log_mu_prior` stores the prior. Pointwise posterior construction adds the prior once per base. Segmentation adds the same prior once per segment.
---

# 2. Serialization

Serializable model objects inherit from `Serializable` in `io.py` and expose:

```python
obj.to_dict()
cls.from_dict(d)
obj.to_npz(path)
cls.from_npz(path)
```

`to_npz()` and `from_npz()` are thin wrappers around `to_dict()` and `from_dict()`. This keeps file serialization and in-memory dictionary serialization consistent.

Typical saved attributes are declared by each class:

```python
class SomeResult(Serializable):
    save_attrs = ("x", "loglik", "log_prior")
    optional_save_attrs = ("theta_x",)
```

`Segmentation` and `LengthPrior` own their own serialization. Saved segmentations are summary-only: they contain posterior state masses, boundary probabilities, and log partitions. They do not save the large traceback caches used by `.sample()` and `.sample_prior()`.

---

# 3. Observation likelihood over sample-level theta

The lowest-level count model is external to this package and is evaluated by `footprint_tools.stats.differential.compute_logpmf_values(...)`. For sample $i$ at position $t$, it provides a custom log likelihood over the theta grid:

$$
\ell_{i,t}(\theta_q)=\log p(y_{i,t}\mid \theta_q).
$$

Here $y_{i,t}$ denotes the observed count information together with expected/background signal and the sample-specific dispersion model. The differential-footprint module does not require this observation likelihood to be Gaussian. It only requires an array of log likelihoods on a finite theta grid.

The theta grid is owned by `DifferentialConfig`, because the same grid is used to build the group-level likelihood and the optional sample-level theta likelihood.

```python
DifferentialConfig(
    mu_min=-6.0,
    mu_max=6.0,
    n_mu=481,
    sig2_min=1e-3,
    sig2_max=10.0,
    n_sig2=181,
    theta_step=0.1,
    theta_tail_sd=5.0,
)
```

If no explicit `theta_grid` is supplied, the theta range is

$$
\left[
\mu_{\min}-a\sqrt{\sigma^2_{\max}},\
\mu_{\max}+a\sqrt{\sigma^2_{\max}}
\right],
$$

where $a=\texttt{theta\_tail\_sd}$, and the spacing is approximately `theta_step`.

---

# 4. Group-level likelihood over mu and sigma2

Implemented in `differential.py`.

## 4.1 Model

For samples in group $g$:

$$
\theta_{i,g,t}\mid \mu_{g,t},\sigma^2_{g,t}
\sim
\mathcal N(\mu_{g,t},\sigma^2_{g,t}).
$$

The count likelihood is supplied externally as $p(y_{i,g,t}\mid\theta_{i,g,t})$. The group-level likelihood integrates out the sample-specific $\theta$ values:

$$
L_{g,t}(\mu,\sigma^2)
=
\prod_{i\in g}
\int
p(y_{i,g,t}\mid\theta)
\,p(\theta\mid\mu,\sigma^2)
\,d\theta.
$$

On the finite theta grid, the Gaussian conditional distribution is converted into bin probabilities:

$$
K_{q,j,k}
=
P(\theta\in B_q\mid \mu_j,\sigma^2_k),
$$

where $B_q$ is the theta-bin around $\theta_q$. The integral is approximated as

$$
\int p(y_i\mid\theta)p(\theta\mid\mu_j,\sigma_k^2)d\theta
\approx
\sum_q p(y_i\mid\theta_q)K_{q,j,k}.
$$

The raw likelihood stored in `loglik_mu_sig2` is

$$
\log L_{g,t}(\mu_j,\sigma^2_k).
$$

The sigma-collapsed likelihood stored in `loglik_mu` is

$$
\log L_{g,t}(\mu_j)
=
\log\sum_k
L_{g,t}(\mu_j,\sigma^2_k)P_g(\sigma^2_k).
$$

## 4.2 Stored object

```python
data.differential.group_names

data.differential.mu_x

data.differential.sig2_x

data.differential.theta_x

data.differential.loglik_mu          # group x position x mu

data.differential.loglik_mu_sig2     # group x position x mu x sigma2

data.differential.log_sig2_prior     # group x sigma2

data.differential.log_mu_prior       # mu
```

## 4.3 Within-group variance prior

For each group, an inverse-chi-square prior is fit empirically to per-position variance estimates of

$$
\log(1+\mathrm{obs})-\log(1+\mathrm{exp}).
$$

The continuous density is evaluated on the `sig2_x` grid and multiplied by quadrature weights. If $w_k$ is the quadrature weight around $\sigma^2_k$, then the stored prior mass is proportional to

$$
P_g(\sigma^2_k)\propto p_{\mathrm{Inv}\chi^2}(\sigma^2_k\mid\nu_g,s_g^2)w_k.
$$

It is normalized on the finite grid. This prior is group-specific and is used when collapsing `loglik_mu_sig2` into `loglik_mu`.

## 4.4 Group mean prior

The prior over $\mu$ is empirical and shared across groups. First, a rough common mean track is obtained by maximizing the summed group likelihood across groups. Then a normal grid prior is fit with mean and standard deviation estimated from that rough common track:

$$
P(\mu_j)\approx
\int_{B_j}
\mathcal N(\mu\mid m_\mu,s_\mu^2)d\mu.
$$

The standard deviation has a floor, controlled by `mu_prior_sd_floor`, to prevent a degenerate prior.

## 4.5 Pointwise group posterior

`data.differential.posterior()` returns

$$
P(\mu_{g,t}=\mu_j\mid y)
\propto
L_{g,t}(\mu_j)P(\mu_j).
$$

This posterior is pointwise and does not impose spatial sharing.

## 4.6 Group mean segmentation

`fit_group_means_segmentation(...)` segments each group mean track using the collapsed likelihood `loglik_mu` and the prior `log_mu_prior`. Within a segment $r=[a,b)$ with state $\mu_j$, the segment score is

$$
\log P(\mu_j)+\sum_{t=a}^{b-1}\log L_{g,t}(\mu_j),
$$

plus the length and transition priors described later.

---

# 5. Common mean and variance-ratio hierarchy

Implemented in `variance_ratio.py`.

## 5.1 Model

The top-level hierarchy relates group means to a shared mean:

$$
\mu_{g,t}\mid \mu_{0,t},\sigma^2_{g,t},\eta_t
\sim
\mathcal N\left(\mu_{0,t},e^{\eta_t}\sigma^2_{g,t}\right).
$$

The variance ratio is

$$
\lambda_t=e^{\eta_t},
$$

and the ICC-like transform is

$$
\mathrm{ICC}_t=\frac{\lambda_t}{1+\lambda_t}.
$$

If the grid includes $\eta=-\infty$, then $\lambda=0$ and $\mathrm{ICC}=0$. This is the exact-consistency state where all group means are equal to $\mu_0$ under the random-effects hierarchy.

## 5.2 Group contribution

For a single group, the contribution to the top-level likelihood is

$$
T_{g,t}(\mu_0,\eta)
=
\sum_{j,k}
L_{g,t}(\mu_j,\sigma^2_k)
P_g(\sigma^2_k)
P(\mu_j\mid\mu_0,e^\eta\sigma^2_k).
$$

The full likelihood is the product across groups:

$$
L_t(\mu_0,\eta)=\prod_g T_{g,t}(\mu_0,\eta).
$$

In log space this is stored as `VarianceRatioLikelihood.loglik` with shape `position x mu0 x eta`.

## 5.3 Priors

The $\mu_0$ prior defaults to the same empirical `log_mu_prior` constructed at the group-mean level.

The $\eta$ prior is a spike-slab grid prior if `consistent_mass` is not `None` and `include_zero=True`:

$$
P(\eta=-\infty)=\pi_0,
$$

where $\pi_0=\texttt{consistent\_mass}$. The remaining mass is assigned to a discretized normal slab:

$$
P(\eta_j)\propto
(1-\pi_0)
\int_{B_j}
\mathcal N(\eta\mid m_\eta,s_\eta^2)d\eta.
$$

The defaults are `consistent_mass=0.30`, `eta_prior_mean=-1.0`, and `eta_prior_sd=2.0`.

If the exact zero state is not included, the prior is just the discretized normal slab.

## 5.4 Exact and Gaussian integration modes

`VarianceRatioConfig.method` controls how the inner integration over $\mu_g$ and $\sigma_g^2$ is approximated.

In exact mode, the calculation uses the full grid likelihood `loglik_mu_sig2` and performs the finite sum over both axes. In Gaussian mode, the group-level posterior is approximated by moment-matched Gaussian terms before evaluating the top-level hierarchy. Exact mode is more faithful to the grid posterior; Gaussian mode is faster and often useful for exploratory work.

The optional `variance_floor` prevents very small variance estimates from creating numerically extreme standardized effects.

## 5.5 Pointwise posterior

`data.variance_ratio.posterior()` returns two marginal posteriors:

$$
P(\mu_{0,t}\mid y)
$$

and

$$
P(\mathrm{ICC}_t\mid y).
$$

The joint posterior before marginalization is

$$
P(\mu_0,\eta\mid y)
\propto
L_t(\mu_0,\eta)P(\mu_0)P(\eta).
$$

---

# 6. Mu0 and ICC segmentation

## 6.1 Mu0 segmentation

`fit_mu0_segmentation(...)` integrates eta independently per base and segments the common mean:

$$
L_t(\mu_0)=\sum_\eta L_t(\mu_0,\eta)P(\eta).
$$

Then a standard segmentation is run over the `mu0_x` states.

## 6.2 ICC / eta segmentation

`fit_eta_segmentation(...)` integrates $\mu_0$ independently per base:

$$
L_t(\eta)=\sum_{\mu_0}L_t(\mu_0,\eta)P(\mu_0).
$$

The segmented state is represented as ICC values, but the prior is the original $\eta$ prior. The state grid stored in the segmentation is

$$
\mathrm{ICC}(\eta)=\frac{e^\eta}{1+e^\eta},
$$

with the convention that $\eta=-\infty$ maps to $\mathrm{ICC}=0$.

After segmenting ICC, the package reconstructs a posterior over $\mu_0$ by conditioning on the segmented ICC posterior. This gives `EtaSegmentation.mu0` and `EtaSegmentation.icc` in one object.

```python
data.eta_segmentation.icc

data.eta_segmentation.mu0

data.eta_segmentation.boundary

data.eta_segmentation.sample(n_draws)
```

---

# 7. Coefficient models

Implemented in `coefficients.py`.

The package currently supports two coefficient definitions. Both use the same `CoefficientLikelihood` container but differ in the reference.

## 7.1 Common-reference coefficient

The common-reference coefficient is

$$
z^{(\mu_0)}_{g,t}=\frac{\mu_{g,t}-\mu_{0,t}}{\sigma_{g,t}}.
$$

This coefficient asks whether group $g$ differs from the shared population mean, using the group-specific within-group scale.

For a target group $g$, the model represents $z_g$ explicitly:

$$
\mu_{g,t}=\mu_{0,t}+z_{g,t}\sigma_{g,t}.
$$

The other groups remain integrated through the top-level hierarchy:

$$
\mu_{h,t}\mid\mu_{0,t},\sigma^2_{h,t},\eta_t
\sim
\mathcal N(\mu_{0,t},e^{\eta_t}\sigma^2_{h,t}),
\qquad h\ne g.
$$

The posterior evidence for $z_g$ is proportional to

$$
\sum_{\mu_0,\eta}
P(\mu_0)P(\eta)P(z_g\mid\eta)
T_{g,t}(\mu_0,z_g)
\prod_{h\ne g}T_{h,t}(\mu_0,\eta).
$$

The target group is excluded from the collapsed product only because it is represented explicitly through $z_g$. It is then added back before integrating over $\mu_0$ and $\eta$.

The prior is induced by eta:

$$
P(z_g)=\sum_\eta P(\eta)P(z_g\mid\eta),
$$

where

$$
z_g\mid\eta\sim\mathcal N(0,e^\eta).
$$

If $\eta=-\infty$, then the conditional prior assigns all mass to $z_g=0$.

The stored `loglik` is prior-corrected:

$$
\log L_t(z)=\log P(z\mid y)-\log P(z)+C_t.
$$

This allows pointwise posterior inference to apply $P(z)$ once per base and segmentation to apply $P(z)$ once per segment.

## 7.2 Zero-reference coefficient

The zero-reference coefficient is

$$
z^{(0)}_{g,t}=\frac{\mu_{g,t}}{\sigma_{g,t}}.
$$

This coefficient asks how many within-group standard deviations the group mean is from zero. It does not use other groups because the reference value is fixed.

For a fixed grid value $z$, the likelihood integrates over the group's variance axis:

$$
L_{g,t}(z)=
\sum_k
L_{g,t}(\mu=z\sqrt{\sigma^2_k},\sigma^2_k)
P_g(\sigma^2_k).
$$

The default prior is a discretized normal prior:

$$
z^{(0)}\sim\mathcal N(0,s_z^2),
$$

with `z_prior_sd=2.0` by default.

## 7.3 Coefficient object and API

```python
data.common_coefficient_likelihood.reference   # "mu0"
data.zero_coefficient_likelihood.reference     # "zero"

data.common_coefficient_likelihood.posterior()
data.zero_coefficient_likelihood.posterior()

data.common_coefficient_segmentation.posterior
data.zero_coefficient_segmentation.posterior
```

The common-reference coefficient is useful for cell-type-specific deviation from the shared mean. The zero-reference coefficient is useful as a standardized depletion/elevation score relative to a fixed no-effect signal.

---

# 8. Coefficient-derived group counts

Only coefficient-derived group counts are retained in v0.15.0.

## 8.1 Zero-reference footprint count

The zero-reference footprint count is

$$
K_{\mathrm{fp,zero},t}
=
\sum_g\mathbf 1\left\{z^{(0)}_{g,t}<-c\right\}.
$$

The default threshold is $c=1$. This counts groups whose zero-reference standardized signal is sufficiently negative.

For each group, the posterior event probability is

$$
p_{g,t}=P(z^{(0)}_{g,t}<-c\mid y).
$$

The count posterior is the Poisson-binomial distribution:

$$
P(K=k\mid y)=
[x^k]
\prod_g\left[(1-p_{g,t})+p_{g,t}x\right].
$$

The loader is:

```python
ZeroFootprintCountLoader()._load(
    data,
    threshold=1.0,
    source="segmented",   # or "pointwise"
)
```

The result is stored as `data.kfp_zero`.

## 8.2 Common-reference deviant count

The deviant count is

$$
K_{\mathrm{fp,dev},t}
=
\sum_g\mathbf 1\left\{|z^{(\mu_0)}_{g,t}|>c\right\}.
$$

For each group,

$$
p_{g,t}=P(|z^{(\mu_0)}_{g,t}|>c\mid y),
$$

and the same Poisson-binomial recursion gives the posterior over $0,\ldots,G$.

The loader is:

```python
DeviantFootprintCountLoader()._load(
    data,
    threshold=1.0,
    source="segmented",   # or "pointwise"
)
```

The result is stored as `data.kfp_dev`.

---

# 9. Sample-level theta likelihood

Implemented in `theta.py`.

`ThetaLikelihood` stores per-sample evidence over theta:

```python
data.theta.sample_names

data.theta.sample_groups

data.theta.group_names

data.theta.theta_x

data.theta.loglik              # sample x position x theta

data.theta.log_theta_prior     # sample x theta

data.theta.mode
```

There are two modes.

## 9.1 Sample-only mode

In `sample_only` mode,

$$
\log L_{i,t}(\theta)=\log p(y_{i,t}\mid\theta).
$$

The sample-specific prior over theta is derived from the corresponding group-level prior predictive distribution.

## 9.2 Group-informed mode

In `group_informed` mode, the target sample is combined with leave-one-sample-out group information. Conceptually:

$$
\log L_{i,t}^{\mathrm{GI}}(\theta)
=
\log p(y_{i,t}\mid\theta)
+
\log p(\theta\mid y_{g\setminus i,t})
-
\log p_0(\theta\mid g).
$$

The subtraction of $\log p_0(\theta\mid g)$ is a prior correction. It ensures that the stored `loglik` can be combined with `log_theta_prior` once for pointwise posterior inference or once per segment for theta segmentation.

The pointwise posterior is

$$
P(\theta_{i,t}\mid y)
\propto
L_{i,t}(\theta)P_i(\theta).
$$

## 9.3 Theta segmentation

Theta segmentation uses `data.theta.loglik` and `data.theta.log_theta_prior` with the same generic explicit-duration segmenter. The segmented state is the sample-level theta value.

```python
data.theta_segmentation.posterior

data.theta_segmentation.boundary

data.theta_segmentation.sample(n_draws)
```

---

# 10. Generic explicit-duration segmentation model

Implemented in `segmentation.py`.

The generic segmenter assumes a one-dimensional state grid $x_1,\ldots,x_M$. For a track $r$ and position $t$, the input likelihood is

$$
\ell_{r,t}(x_m)=\log p(y_{r,t}\mid x_m).
$$

A segmentation consists of segments

$$
([a_1,b_1),s_1),\ldots,([a_R,b_R),s_R),
$$

where $s_j$ is a state index and $b_j-a_j$ is the segment length.

For a complete segment $[a,b)$ with state $s$, the emission score is

$$
\log P(x_s)+\sum_{t=a}^{b-1}\ell_t(x_s),
$$

where $P(x_s)$ is the state prior. The state prior is applied once per segment.

## 10.1 Length prior

`LengthPrior` stores explicit log masses for lengths $1,\ldots,L$ and an optional geometric tail. If the explicit mass for length $L$ is $q_L$ and the tail ratio is $r$, then for lengths beyond $L$:

$$
q_{L+h}=q_L r^h,
\qquad h\ge 1.
$$

The constructor normalizes the finite support plus tail so that

$$
\sum_{\ell\ge1}q_\ell=1.
$$

The length prior may be shared across states or state-specific.

```python
length_prior = LengthPrior.from_log_mass(log_mass, infer_tail=True)
length_prior.to_dict()
LengthPrior.from_dict(d)
```

To forbid short segments, set their log mass to $-\infty$ before constructing the `LengthPrior`.

## 10.2 Stationary finite-window correction

The segmenter treats the observed DHS/window as a finite window cut out of a stationary segmentation process. Therefore, the first observed segment may have started before position zero, and the last observed segment may continue past the right edge.

If $q_k(\ell)$ is the complete length distribution in state $k$, its mean is

$$
\bar L_k=\sum_{\ell\ge1}\ell q_k(\ell).
$$

The residual length distribution at a random position is proportional to the survival function:

$$
P(R_k\ge r)\propto P(L_k\ge r).
$$

This stationary correction avoids forcing a segment boundary at the first base of every DHS.

## 10.3 State transitions

The transition prior can be independent of the previous state or locally attractive.

If `transition_sd=None`, the next state is sampled from the state prior. If `forbid_same_state=True`, transitions back to the same state are removed and the rows are renormalized.

If `transition_sd` is finite, transitions prefer nearby states according to a Gaussian kernel on the state grid:

$$
T_{ab}\propto P(x_b)\exp\left[-\frac{(x_b-x_a)^2}{2\tau_T^2}\right].
$$

The implementation uses exact dynamic programming. Independent transitions use an $O(M)$ update. Narrow Gaussian transitions use a sparse/banded representation rather than a fully dense $M^2$ matrix when possible.

## 10.4 Posterior outputs

`segment(...)` returns a `Segmentation`:

```python
seg.names
seg.posterior       # GridPosterior over states, track x position x state
seg.boundary        # track x (position - 1)
seg.log_partition   # track
seg.sample(n_draws)
seg.sample_prior(n_draws)
```

`seg.boundary[i, t]` is

$$
P(\text{boundary between }t\text{ and }t+1\mid y).
$$

The expected number of boundaries in a track is therefore

$$
\sum_t P(B_t=1\mid y).
$$

Posterior path sampling is exact under the fitted discretized segmentation model. It uses stochastic traceback through cached backward messages. A segmentation loaded from NPZ is summary-only and cannot sample new paths because the traceback messages are intentionally not saved.

---

# 11. Priors by object

The following table summarizes where priors enter the current package.

| Object | Prior | Default/source | Applied |
|---|---|---|---|
| `Differential` | $P_g(\sigma^2)$ | Empirical inverse-chi-square per group | During collapse of `loglik_mu_sig2` to `loglik_mu` |
| `Differential` | $P(\mu)$ | Empirical normal grid prior from rough common mean track | Pointwise group posterior and group mean segmentation |
| `VarianceRatioLikelihood` | $P(\mu_0)$ | Defaults to `Differential.log_mu_prior` | Pointwise top-level posterior and mu0/eta segmentation |
| `VarianceRatioLikelihood` | $P(\eta)$ | Spike-slab with optional exact $\eta=-\infty$ state | Pointwise top-level posterior and eta/mu0 segmentation |
| `CommonCoefficientModel` | $P(z\mid\eta)$ and induced $P(z)$ | Induced by eta hierarchy | Common-reference coefficient posterior and segmentation |
| `ZeroCoefficientModel` | $P(z)$ | Discretized normal, default sd 2 | Zero-reference coefficient posterior and segmentation |
| `ThetaLikelihood` | $P_i(\theta)$ | Group prior predictive theta distribution | Theta posterior and theta segmentation |
| `Segmentation` | $P(L)$ | User-supplied `LengthPrior` | Segment duration prior |
| `Segmentation` | $T(s'\mid s)$ | Independent or Gaussian state transition | Between consecutive segments |

---

# 12. Main interpretation guide

## Group means

`data.differential.posterior()` and `data.segmentation.posterior` describe the absolute group signal $\mu_g$. These tracks can be sensitive to reproducible model misspecification because, with enough samples, even small systematic deviations can become very certain.

## ICC / eta

`data.eta_segmentation.icc` describes the degree of between-group heterogeneity relative to within-group variance. It is a variance-ratio track, not a footprint-depth track. It is often smoother than group mean tracks because the fine common mean structure is integrated out when segmenting eta.

## Common-reference coefficients

`data.common_coefficient_likelihood` and `data.common_coefficient_segmentation` describe group deviations from the shared mean:

$$
z^{(\mu_0)}_g=\frac{\mu_g-\mu_0}{\sigma_g}.
$$

These are useful for attributing which groups are high or low relative to the population of groups.

## Zero-reference coefficients

`data.zero_coefficient_likelihood` and `data.zero_coefficient_segmentation` describe standardized signal relative to zero:

$$
z^{(0)}_g=\frac{\mu_g}{\sigma_g}.
$$

These are useful for thresholding depletion or elevation relative to a fixed no-effect point.

## Count summaries

`data.kfp_zero` counts groups with sufficiently negative zero-reference coefficients. `data.kfp_dev` counts groups with sufficiently large common-reference deviations in either direction.

Both are posterior distributions over counts, not just point estimates.

---

# 13. Bayesian event reporting

For coefficient inference, define biologically meaningful events rather than testing exact equality to zero. Examples:

$$
E^+_{g,t}(c)=\{z_{g,t}>c\},
$$

$$
E^-_{g,t}(c)=\{z_{g,t}<-c\},
$$

or

$$
E^{\mathrm{abs}}_{g,t}(c)=\{|z_{g,t}|>c\}.
$$

The posterior event probability is directly available from `GridPosterior` methods. The local posterior error for a selected event is

$$
\ell_i=1-P(E_i\mid y).
$$

For a selected set $S$, the Bayesian false-discovery rate is

$$
\mathrm{BFDR}(S)=\frac{1}{|S|}\sum_{i\in S}\ell_i.
$$

This is the Bayesian analogue of controlling an expected false discovery fraction, but it is conditional on the model. Posterior probabilities are not p-values and are not expected to be uniformly distributed under a null. Calibration should instead be checked through prior predictive checks, posterior predictive checks, simulation-based calibration, and empirical negative controls.

---

# 14. Minimal workflow example

```python
from differential_footprint.api import (
    DifferentialLoader,
    GroupMeanSegmentationLoader,
    VarianceRatioLoader,
    EtaSegmentationLoader,
    Mu0SegmentationLoader,
    CommonCoefficientLikelihoodLoader,
    CommonCoefficientSegmentationLoader,
    ZeroCoefficientLikelihoodLoader,
    ZeroCoefficientSegmentationLoader,
    ZeroFootprintCountLoader,
    DeviantFootprintCountLoader,
)

# 1. Group-level likelihood
data = DifferentialLoader()._load(data)

# 2. Optional group mean segmentation
data = GroupMeanSegmentationLoader()._load(data, length_prior=length_prior)

# 3. Common mean / variance-ratio hierarchy
data = VarianceRatioLoader()._load(data)
data = EtaSegmentationLoader()._load(data, length_prior=length_prior)
data = Mu0SegmentationLoader()._load(data, length_prior=length_prior)

# 4. Coefficients
data = CommonCoefficientLikelihoodLoader()._load(data)
data = CommonCoefficientSegmentationLoader()._load(data, length_prior=length_prior)

data = ZeroCoefficientLikelihoodLoader()._load(data)
data = ZeroCoefficientSegmentationLoader()._load(data, length_prior=length_prior)

# 5. Coefficient-derived counts
data = ZeroFootprintCountLoader()._load(data, threshold=1.0, source="segmented")
data = DeviantFootprintCountLoader()._load(data, threshold=1.0, source="segmented")
```

---

# 15. Saved and unsaved state

The following are saved by NPZ serialization:

```text
Likelihood objects:
  grids
  log likelihood arrays
  log priors
  names
  method/reference metadata

Posterior objects:
  state grid
  normalized log posterior mass

Segmentation objects:
  names
  state grid
  posterior log mass
  boundary probabilities
  log partition

LengthPrior:
  normalized explicit log masses
  geometric tail ratios
```

The following are not saved:

```text
Segmentation traceback caches
Forward/reverse dynamic-programming work arrays
Raw count data
Temporary theta integration buffers
```

Therefore, freshly fitted segmentation objects can draw posterior paths with `.sample()`. Segmentations loaded from NPZ remain valid summary objects but cannot generate new paths without refitting the segmentation.
