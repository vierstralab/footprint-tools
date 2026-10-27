# Supplementary Methods: Bayesian hierarchical inference and spatial segmentation of differential footprints

## Overview

We developed a hierarchical Bayesian model to quantify regulatory footprint signal at three related levels: individual samples, cell-type or sample groups, and the population of groups. The model separates the uncertainty associated with sequencing counts from biological variability among samples within a group and heterogeneity among groups. Spatial consistency along the genomic coordinate is incorporated using an explicit-duration Bayesian segmentation model rather than by smoothing point estimates.

Samples are indexed by $i$, groups by $g=1,\ldots,G$, and genomic positions within the analyzed interval by $t=1,\ldots,N$. The principal latent variables are

$$
\theta_{igt},
$$

the latent signal of sample $i$ belonging to group $g$;

$$
\mu_{gt},
$$

the mean signal of group $g$;

$$
\sigma^2_{gt},
$$

the variance among samples within group $g$;

$$
\mu_{0t},
$$

the common mean across groups;

$$
\lambda_t=e^{\eta_t},
$$

the ratio between among-group and within-group variance; and

$$
\operatorname{ICC}_t=\frac{\lambda_t}{1+\lambda_t},
$$

a bounded transformation of this variance ratio.

Two standardized group effects are additionally defined:

$$
z^{(\mu_0)}_{gt}=\frac{\mu_{gt}-\mu_{0t}}{\sigma_{gt}},
$$

which measures deviation of group $g$ from the common mean, and

$$
z^{(0)}_{gt}=\frac{\mu_{gt}}{\sigma_{gt}},
$$

which measures the signal relative to zero in units of the within-group standard deviation.

All continuous latent quantities are represented numerically on finite one-dimensional grids. Thus, all reported posterior probabilities are posterior probability masses associated with grid bins rather than pointwise evaluations of continuous posterior densities.

---

## 1. Numerical grids and discretization of continuous distributions

The default group-mean grid contains 481 equally spaced values,

$$
\mu_j=-6+0.025j,\qquad j=0,\ldots,480,
$$

and therefore spans

$$
-6\leq\mu\leq6.
$$

The within-group variance grid consists of 181 geometrically spaced values between

$$
10^{-3}\leq\sigma^2\leq10.
$$

Specifically,

$$
\sigma_k^2=10^{-3}\left(\frac{10}{10^{-3}}\right)^{k/180},\qquad k=0,\ldots,180.
$$

The ratio between adjacent variance-grid values is approximately 1.0525003.

Unless an explicit sample-level grid is supplied, the $\theta$ grid is generated from

$$
\theta_{\min}=\mu_{\min}-a\sqrt{\sigma^2_{\max}},
$$

$$
\theta_{\max}=\mu_{\max}+a\sqrt{\sigma^2_{\max}},
$$

with $a=5$. Under the default settings this gives

$$
\theta_{\min}\approx-21.8113883,
$$

$$
\theta_{\max}\approx21.8113883.
$$

The requested nominal spacing is 0.1. The implementation first computes

$$
Q=\left\lceil\frac{\theta_{\max}-\theta_{\min}}{0.1}\right\rceil+1=438
$$

and subsequently uses `linspace`; hence the actual default spacing is approximately

$$
\Delta_\theta=\frac{\theta_{\max}-\theta_{\min}}{437}\approx0.0998233,
$$

rather than exactly 0.1.

The finite part of the $\eta$ grid is

$$
\eta_k=-4+0.25k,\qquad k=0,\ldots,32.
$$

By default an additional state $\eta=-\infty$ is prepended. It corresponds exactly to

$$
\lambda=e^\eta=0
$$

and therefore

$$
\operatorname{ICC}=0.
$$

Thus the default variance-ratio model contains 34 states: one exact consistency state and 33 finite $\eta$ states.

The standardized coefficient grid consists of 33 points

$$
z=-4,-3.75,\ldots,3.75,4.
$$

The grid point nearest zero is explicitly replaced by exactly 0 to avoid floating-point representations of zero.

### Discretization of normal distributions

A continuous normal distribution is converted into probability mass on a grid $x_1,\ldots,x_M$ by defining bin boundaries

$$
b_0=-\infty,
$$

$$
b_j=\frac{x_j+x_{j+1}}{2},\qquad j=1,\ldots,M-1,
$$

and

$$
b_M=+\infty.
$$

For a normal random variable $X\sim\mathcal N(m,s^2)$, the mass associated with grid point $x_j$ is therefore

$$
p_j=\Phi\left(\frac{b_j-m}{s}\right)-\Phi\left(\frac{b_{j-1}-m}{s}\right).
$$

Consequently, probability outside the numerical range of the grid is assigned to the two terminal grid bins rather than discarded. Numerical zero probabilities are clipped at the smallest positive floating-point number before logarithms are taken, and the resulting masses are normalized to sum to one.

---

## 2. Sample-level observation likelihood

For each sample $i$ and genomic position $t$, the count model provides a log likelihood evaluated on the $\theta$ grid,

$$
\ell_{itq}=\log p(y_{it}\mid\theta_q).
$$

Here $y_{it}$ denotes the observed footprint count information together with the corresponding expected/background signal and the sample-specific dispersion model.

The implementation obtains these values using `compute_logpmf_values` from the external `footprint_tools.stats.differential` module. The resulting array has dimensions

$$
Q\times N\times I,
$$

corresponding to $\theta$ state, genomic position, and sample, respectively. No Gaussian approximation to this observation likelihood is required by the hierarchical model.

**Reproducibility note.** The exact probability mass function evaluated by `compute_logpmf_values` is external to the present module. Thus, exact reproduction starting from raw observed and expected counts additionally requires the implementation or mathematical definition of that function. From the array $\ell_{itq}$ onward, the inference procedure is completely specified below.

---

## 3. Group-level hierarchical model

### 3.1 Model for variation among samples

Within each group,

$$
\theta_{igt}\mid\mu_{gt},\sigma^2_{gt}\sim\mathcal N(\mu_{gt},\sigma^2_{gt}).
$$

For a fixed group mean $\mu_j$ and variance $\sigma_k^2$, define

$$
K_{qjk}=P\left(\theta\in B_q\mid\mu_j,\sigma_k^2\right),
$$

where $B_q$ is the bin associated with $\theta_q$. Explicitly,

$$
K_{qjk}=\Phi\left(\frac{b_{q+1}-\mu_j}{\sqrt{\sigma_k^2}}\right)-\Phi\left(\frac{b_q-\mu_j}{\sqrt{\sigma_k^2}}\right).
$$

The implementation renormalizes $K_{\cdot jk}$ after computing the CDF differences to ensure

$$
\sum_qK_{qjk}=1
$$

to numerical precision.

For sample $i$, the likelihood conditional on $(\mu_j,\sigma_k^2)$ is approximated by the finite-grid sum

$$
A_{itjk}=\sum_q\exp(\ell_{itq})K_{qjk}.
$$

Samples are conditionally independent given the group parameters, giving

$$
L_{gt}(\mu_j,\sigma_k^2)=\prod_{i\in g}A_{itjk}.
$$

Equivalently,

$$
\log L_{gt}(\mu_j,\sigma_k^2)=\sum_{i\in g}\log A_{itjk}.
$$

These values are stored without a prior on $\mu$ or $\sigma^2$.

All sums are evaluated using max-shifted exponentiation or log-sum-exp operations to avoid numerical underflow.

### 3.2 Empirical prior for within-group variance

A separate prior on $\sigma^2$ is estimated for every group.

For sample $i$ and position $t$, a preliminary log-ratio statistic is calculated as

$$
r_{it}=\log(1+\operatorname{obs}_{it})-\log(1+\operatorname{exp}_{it}).
$$

For each group $g$, the sample variance of this statistic is calculated independently at every position:

$$
v_{gt}=\operatorname{Var}_{i\in g}(r_{it}),
$$

using the unbiased sample-variance convention with denominator $n_g-1$. Non-finite variance estimates and estimates satisfying

$$
v_{gt}\leq0
$$

are discarded.

An inverse-$\chi^2$ distribution is then fit to the remaining $v_{gt}$ values by maximum likelihood. If its parameters are denoted $\nu_g$ and $s_g^2$, the optimization is performed over

$$
\alpha=\log\nu_g,
$$

$$
\beta=\log s_g^2.
$$

The objective is

$$
-\sum_t\log p_{\operatorname{Inv}\chi^2}(v_{gt}\mid\nu_g,s_g^2).
$$

Optimization uses L-BFGS-B initialized at

$$
\nu_g=5,
$$

$$
s_g^2=\operatorname{median}_t(v_{gt}),
$$

with bounds

$$
-8\leq\log\nu_g\leq12
$$

and

$$
-20\leq\log s_g^2\leq12.
$$

The exact analytic parameterization of the inverse-$\chi^2$ density is inherited from the external `footprint_tools.stats.distributions.invchi2` implementation.

The fitted continuous density is converted into discrete masses on the $\sigma^2$ grid using quadrature weights. For a grid $x_k=\sigma_k^2$,

$$
w_1=\frac{x_2-x_1}{2},
$$

$$
w_K=\frac{x_K-x_{K-1}}{2},
$$

and for internal grid points

$$
w_k=\frac{x_{k+1}-x_{k-1}}{2}.
$$

The discrete variance prior is

$$
\pi^{(\sigma)}_{gk}=\frac{p_{\operatorname{Inv}\chi^2}(\sigma_k^2\mid\nu_g,s_g^2)w_k}{\sum_{k'}p_{\operatorname{Inv}\chi^2}(\sigma_{k'}^2\mid\nu_g,s_g^2)w_{k'}}.
$$

Thus the inverse-$\chi^2$ distribution is effectively truncated to the numerical range

$$
10^{-3}\leq\sigma^2\leq10
$$

and renormalized over that grid.

### 3.3 Marginal likelihood for the group mean

Within-group variance is integrated out to obtain

$$
L_{gt}(\mu_j)=\sum_kL_{gt}(\mu_j,\sigma_k^2)\pi^{(\sigma)}_{gk}.
$$

The implementation stores both

$$
\log L_{gt}(\mu_j,\sigma_k^2)
$$

and

$$
\log L_{gt}(\mu_j).
$$

Importantly, neither quantity contains the prior on $\mu$.

### 3.4 Empirical prior for group means

A preliminary common mean is first estimated independently at every genomic position by maximizing the summed likelihood across groups:

$$
\widetilde{\mu}_{0t}=\arg\max_{\mu_j}\sum_g\log L_{gt}(\mu_j).
$$

No $\mu$ prior is used in this step.

An empirical normal prior is then constructed from the distribution of $\widetilde{\mu}_{0t}$ over positions. Its continuous mean and standard deviation are

$$
m_\mu=\frac{1}{N}\sum_t\widetilde{\mu}_{0t}
$$

and

$$
s_\mu=\max\left[\operatorname{SD}_{t}(\widetilde{\mu}_{0t}),s_{\min}\right],
$$

where NumPy's population-standard-deviation convention is used. Unless an explicit floor is supplied,

$$
s_{\min}=\operatorname{median}_j(\mu_{j+1}-\mu_j)=0.025.
$$

The continuous prior

$$
\mu\sim\mathcal N(m_\mu,s_\mu^2)
$$

is converted to grid-bin masses as described above. This single prior is shared by all groups and is subsequently also used as the default prior for the common mean $\mu_0$.

The pointwise posterior for each group is therefore

$$
P(\mu_{gt}=\mu_j\mid y)=\frac{L_{gt}(\mu_j)\pi^{(\mu)}_j}{\sum_{j'}L_{gt}(\mu_{j'})\pi^{(\mu)}_{j'}}.
$$

No spatial information is used at this stage.

---

## 4. Common-mean and between-group variance hierarchy

To distinguish differences among groups from variation among samples within a group, group means are modeled as

$$
\mu_{gt}\mid\mu_{0t},\sigma^2_{gt},\lambda_t\sim\mathcal N\left(\mu_{0t},\lambda_t\sigma^2_{gt}\right).
$$

We parameterize

$$
\eta_t=\log\lambda_t
$$

and report

$$
\operatorname{ICC}_t=\frac{\lambda_t}{1+\lambda_t}.
$$

The term “ICC” is used here as a convenient bounded variance-ratio transformation; the modeled between-group variance for group $g$ is explicitly $\lambda_t\sigma^2_{gt}$.

### 4.1 Exact finite-grid likelihood

For a finite value of $\eta_k$,

$$
\lambda_k=e^{\eta_k}>0.
$$

For each $\sigma_l^2$, the conditional normal distribution

$$
\mu_g\mid\mu_0,\sigma_l^2,\lambda_k\sim\mathcal N(\mu_0,\lambda_k\sigma_l^2)
$$

is converted to probability masses over the $\mu$ bins. Define

$$
H_{ljuk}=P\left(\mu_g\in B_j^{(\mu)}\mid\mu_0=\mu_u,\lambda_k\sigma_l^2\right).
$$

The contribution of group $g$ to the population-level likelihood is then

$$
T_{gt}(\mu_u,\eta_k)=\sum_l\sum_jL_{gt}(\mu_j,\sigma_l^2)\pi^{(\sigma)}_{gl}H_{ljuk}.
$$

For the exact consistency state

$$
\eta=-\infty,\qquad\lambda=0,
$$

the normal conditional degenerates at $\mu_g=\mu_0$. Rather than numerically approximating a zero-variance normal distribution, the implementation uses

$$
T_{gt}(\mu_u,-\infty)=L_{gt}(\mu_u),
$$

with linear interpolation of the log likelihood if the $\mu_0$ and $\mu$ grids differ. In the default implementation they are the same grid.

The likelihood combining all groups is

$$
L_t(\mu_u,\eta_k)=\prod_{g=1}^{G}T_{gt}(\mu_u,\eta_k),
$$

or equivalently

$$
\log L_t(\mu_u,\eta_k)=\sum_g\log T_{gt}(\mu_u,\eta_k).
$$

The implementation performs these sums directly on the full $(\mu,\sigma^2)$ likelihood grid when `method="exact"`.

### 4.2 Gaussian approximation

For faster computation, the default variance-ratio model uses a Gaussian approximation to each group likelihood.

First, the marginalized group likelihood over $\mu$ is treated as a continuous likelihood sampled on the $\mu$ grid. Let $w_j^{(\mu)}$ denote the trapezoidal quadrature weights. Define

$$
Z_{gt}=\sum_jL_{gt}(\mu_j)w_j^{(\mu)},
$$

$$
\widehat{\mu}_{gt}=\frac{\sum_j\mu_jL_{gt}(\mu_j)w_j^{(\mu)}}{Z_{gt}},
$$

and

$$
V^{(\mu)}_{gt}=\frac{\sum_j(\mu_j-\widehat{\mu}_{gt})^2L_{gt}(\mu_j)w_j^{(\mu)}}{Z_{gt}}.
$$

This calculation uses the likelihood over $\mu$, not the empirical prior on $\mu$.

To prevent a nearly degenerate approximation, the variance is bounded below. If no explicit floor is supplied,

$$
V_{\min}=\frac{\Delta_\mu^2}{12}.
$$

With the default $\Delta_\mu=0.025$,

$$
V_{\min}\approx5.2083\times10^{-5}.
$$

Thus

$$
V^{(\mu)}_{gt}\leftarrow\max(V^{(\mu)}_{gt},V_{\min}).
$$

A posterior-like distribution over $\sigma^2$ is calculated independently as

$$
P(\sigma_l^2\mid y_{gt})\propto\pi^{(\sigma)}_{gl}\sum_jL_{gt}(\mu_j,\sigma_l^2)w_j^{(\mu)},
$$

and

$$
\widehat{\sigma^2}_{gt}=\sum_lP(\sigma_l^2\mid y_{gt})\sigma_l^2.
$$

For $\lambda_k>0$,

$$
T_{gt}(\mu_0,\eta_k)\approx Z_{gt}\,\mathcal N\left(\widehat{\mu}_{gt};\mu_0,V^{(\mu)}_{gt}+\lambda_k\widehat{\sigma^2}_{gt}\right).
$$

The exact-consistency state $\lambda=0$ remains special: the original marginalized likelihood $L_{gt}(\mu_0)$ is used directly rather than the Gaussian approximation.

### 4.3 Prior on the variance ratio

The common mean uses the empirical prior described above,

$$
P(\mu_0=\mu_j)=\pi^{(\mu)}_j.
$$

By default, the $\eta$ prior contains an exact consistency component with probability

$$
P(\eta=-\infty)=0.30.
$$

The remaining probability 0.70 is distributed over the finite $\eta$ grid according to bin masses of

$$
\eta\sim\mathcal N(-1,2^2).
$$

Thus,

$$
P(\eta_k)=0.7\,P_{\mathrm{slab}}(\eta_k),
$$

for finite states, while the exact state receives mass 0.3. The normal slab is first discretized and normalized over the finite grid, using open-ended terminal bins as described above.

### 4.4 Pointwise posterior

The pointwise joint posterior is

$$
P(\mu_{0t}=\mu_u,\eta_t=\eta_k\mid y)\propto L_t(\mu_u,\eta_k)\pi^{(\mu)}_u\pi^{(\eta)}_k.
$$

It is normalized jointly over both state axes.

The common-mean posterior is obtained by summing over $\eta$,

$$
P(\mu_{0t}=\mu_u\mid y)=\sum_kP(\mu_{0t}=\mu_u,\eta_t=\eta_k\mid y),
$$

while the variance-ratio posterior is obtained by summing over $\mu_0$. Its state labels are subsequently transformed from $\eta_k$ into

$$
\operatorname{ICC}_k=\frac{\lambda_k}{1+\lambda_k}.
$$

---

## 5. Explicit-duration spatial segmentation

Pointwise inference treats genomic positions independently. Spatially coherent posterior estimates are obtained using an explicit-duration hidden semi-Markov model in which a latent state remains constant over segments of variable length.

Let the state grid for a particular analysis be

$$
x_1,\ldots,x_M.
$$

For example, $x$ may represent $\mu$, $\operatorname{ICC}$, $z$, or $\theta$.

For each position $t$, the segmentation algorithm receives an emission log likelihood

$$
\ell_t(s)=\log p(y_t\mid x_s).
$$

A segmentation consists of half-open segments. For a segment occupying

$$
[a,b)
$$

with state $s$, define the cumulative emission

$$
C_s(a,b)=\sum_{t=a}^{b-1}\ell_t(s).
$$

Prefix sums are used so that $C_s(a,b)$ can be evaluated in constant time.

### 5.1 Segment-length distribution

For each state $s$, the model specifies explicit unnormalized log weights

$$
a_s(1),\ldots,a_s(L_0).
$$

A single set of weights may be supplied and broadcast to every state, or separate length distributions may be supplied for individual states.

An optional geometric tail extends the distribution beyond $L_0$. Let

$$
r_s=e^{\rho_s},\qquad0\leq r_s<1.
$$

If

$$
w_s(l)=\exp[a_s(l)-\max_ja_s(j)],
$$

the normalization constant is

$$
Z_s=\sum_{l=1}^{L_0}w_s(l)+w_s(L_0)\frac{r_s}{1-r_s}.
$$

The resulting complete duration distribution is

$$
q_s(l)=\frac{w_s(l)}{Z_s},\qquad1\leq l\leq L_0,
$$

and

$$
q_s(L_0+h)=q_s(L_0)r_s^h,\qquad h\geq1.
$$

If no geometric tail is requested, $r_s=0$.

When the tail ratio is inferred from the supplied log masses, the implementation calculates

$$
\rho_s=\frac{1}{n_{\mathrm{tail}}-1}\sum_{j=L_0-n_{\mathrm{tail}}+1}^{L_0-1}\left[a_s(j+1)-a_s(j)\right].
$$

The inferred value is required to satisfy

$$
\rho_s<0.
$$

The current default is $n_{\mathrm{tail}}=2$, although a larger value can be supplied. Segment lengths that should be impossible are represented by assigning their input log mass to $-\infty$.

Define the state-specific survival function

$$
S_s(v)=P(L_s\geq v)=\sum_{l=v}^{\infty}q_s(l),
$$

the mean duration

$$
m_s=E[L_s]=\sum_{l=1}^{\infty}lq_s(l),
$$

and

$$
E_s(v)=\sum_{l=v}^{\infty}(l-v+1)q_s(l).
$$

These quantities are computed analytically for the geometric tail rather than by truncating it numerically.

### 5.2 State transitions

Let $\pi_s$ denote the supplied state prior.

When `transition_sd=None` and transitions to the same state are permitted,

$$
T_{ab}=\pi_b.
$$

Thus the state of each newly initiated segment is independent of the previous segment state.

When transitions to the same state are forbidden,

$$
T_{ab}=\begin{cases}\dfrac{\pi_b}{1-\pi_a},&a\neq b,\\[6pt]0,&a=b.\end{cases}
$$

The stationary distribution of this transition matrix is then

$$
\nu_s\propto\pi_s(1-\pi_s).
$$

The implementation calculates this expression directly for the independent no-self-transition model.

If `transition_sd=\tau_T` is finite, a locally attractive transition model is instead used:

$$
T_{ab}\propto\pi_b\exp\left[-\frac{(x_b-x_a)^2}{2\tau_T^2}\right].
$$

The diagonal is set to zero before normalization when same-state transitions are forbidden. Each row is then normalized to sum to one.

For this general transition matrix, its stationary distribution $\nu$ is computed by power iteration starting from the uniform distribution,

$$
\nu^{(r+1)}=\nu^{(r)}T,
$$

until

$$
\max_s|\nu^{(r+1)}_s-\nu^{(r)}_s|<10^{-14},
$$

or 10,000 iterations have been performed.

### 5.3 Stationary finite-window correction

A crucial feature of the segmentation model is that the analyzed genomic interval is treated as a finite window cut from an underlying stationary segmentation process. Consequently, a segment visible at the left edge is not assumed to begin at the first analyzed nucleotide, and a segment visible at the right edge is not assumed to terminate at the last analyzed nucleotide.

Let

$$
\bar m=\sum_s\nu_sm_s
$$

be the mean segment length under the stationary distribution of segment states.

Consider first a segmentation containing a single observed segment with state $s$ that spans an interval of length $N$. Its prior weight is

$$
\frac{\nu_sE_s(N)}{\bar m}.
$$

Now consider a segmentation with $R\geq2$ observed segments,

$$
[0,b_1),[b_1,b_2),\ldots,[b_{R-1},N),
$$

having states

$$
s_1,\ldots,s_R.
$$

Its prior weight is

$$
\frac{\nu_{s_1}S_{s_1}(b_1)}{\bar m}\left[\prod_{r=2}^{R-1}T_{s_{r-1}s_r}q_{s_r}(b_r-b_{r-1})\right]T_{s_{R-1}s_R}S_{s_R}(N-b_{R-1}).
$$

The first survival term arises because the first observed segment may have begun before the left edge. The final survival term arises because the final observed segment may continue beyond the right edge.

After including sequence evidence, the unnormalized posterior weight of the path is obtained by multiplying the above expression by

$$
\exp\left[\sum_{r=1}^{R}C_{s_r}(b_{r-1},b_r)\right].
$$

This finite-window construction is why the method does not artificially impose boundaries at the beginning or end of every analyzed DHS or footprint interval.

### 5.4 Exact marginalization over segmentations

The posterior partition function is

$$
Z=\sum_{\mathcal S}W(\mathcal S),
$$

where the sum is over every admissible segmentation path $\mathcal S$ and $W(\mathcal S)$ is its unnormalized weight.

The implementation evaluates this sum exactly under the discretized model using forward and reverse dynamic programming in log space. Prefix sums of the emission likelihoods permit segment emission terms to be evaluated without repeatedly summing over their constituent bases. Contributions from segment lengths beyond the explicitly represented duration range are accumulated analytically using the geometric-tail ratio.

The marginal state posterior is

$$
P(X_t=x_s\mid y)=\frac{\sum_{\mathcal S}W(\mathcal S)\mathbf 1\{X_t(\mathcal S)=x_s\}}{Z}.
$$

Similarly, for every adjacent pair of positions,

$$
P(B_t=1\mid y)=\frac{\sum_{\mathcal S}W(\mathcal S)\mathbf 1\{\text{a segment boundary occurs between }t\text{ and }t+1\}}{Z}.
$$

The algorithm returns both quantities together with

$$
\log Z.
$$

Boundary probabilities are calculated directly from matching forward and reverse messages and are clipped at one only to protect against floating-point roundoff.

Posterior segmentation paths can additionally be sampled by stochastic traceback. At each sampled boundary, the next state and segment endpoint are drawn from the exact conditional distributions implied by the reverse messages. Thus path sampling introduces Monte Carlo variability only into the sampled path representation; the reported marginal state and boundary probabilities themselves are calculated analytically by dynamic programming.

---

## 6. Segmentation of group means, the common mean, and ICC

### 6.1 Group mean segmentation

Each group is segmented independently using the marginalized group likelihood

$$
L_{gt}(\mu_j).
$$

The state grid is the $\mu$ grid and the state prior is $\pi^{(\mu)}$.

The default group-mean segmentation uses independent transitions with same-state transitions permitted:

$$
\texttt{transition\_sd=None},
$$

$$
\texttt{forbid\_same\_state=False}.
$$

Thus, apart from the duration model and stationary finite-window correction, new segment states are drawn from the empirical $\mu$ prior.

### 6.2 Common-mean segmentation

For segmentation of $\mu_0$, the variance-ratio state is marginalized independently at each position:

$$
L_t^{(\mu_0)}(\mu_u)=\sum_kL_t(\mu_u,\eta_k)\pi^{(\eta)}_k.
$$

The resulting likelihood is supplied to the same explicit-duration segmenter with state grid $\mu_u$ and state prior $\pi^{(\mu)}_u$.

### 6.3 ICC segmentation

For segmentation of between-group heterogeneity, $\mu_0$ is instead integrated independently at every position:

$$
L_t^{(\eta)}(\eta_k)=\sum_uL_t(\mu_u,\eta_k)\pi^{(\mu)}_u.
$$

The segmentation state labels are

$$
x_k=\operatorname{ICC}(\eta_k),
$$

while the state probabilities supplied to the segmenter are the corresponding $\eta$-prior masses.

Under the default settings,

$$
\texttt{transition\_sd=None},
$$

$$
\texttt{forbid\_same\_state=True}.
$$

Hence successive segments cannot use the same ICC state. If a finite transition standard deviation is supplied, distances between states are calculated on the ICC grid, rather than on the original $\eta$ scale, because ICC values are passed to the generic transition model.

After segmentation, the common-mean posterior is reconstructed without rerunning a joint segmentation. At every position,

$$
P_t(\mu_u\mid\eta_k,y)=\frac{L_t(\mu_u,\eta_k)\pi^{(\mu)}_u}{\sum_{u'}L_t(\mu_{u'},\eta_k)\pi^{(\mu)}_{u'}}.
$$

If

$$
P_{\mathrm{seg}}(\eta_k\mid y)
$$

is the spatially segmented posterior, the reconstructed common-mean posterior is

$$
P_{\mathrm{seg}}(\mu_u\mid y)=\sum_kP_{\mathrm{seg}}(\eta_k\mid y)P_t(\mu_u\mid\eta_k,y).
$$

Thus the segmentation acts spatially on the heterogeneity parameter, while the fine-scale conditional uncertainty in $\mu_0$ is retained.

---

## 7. Standardized group coefficients

Two standardized effect variables are calculated.

### 7.1 Group deviation from the common mean

For each group,

$$
z^{(\mu_0)}_{gt}=\frac{\mu_{gt}-\mu_{0t}}{\sigma_{gt}}.
$$

Equivalently, conditional on a proposed value $z$,

$$
\mu_{gt}=\mu_{0t}+z\sigma_{gt}.
$$

Under the variance-ratio hierarchy,

$$
z_{gt}\mid\eta_t\sim\mathcal N(0,e^{\eta_t}),
$$

because

$$
\mu_{gt}-\mu_{0t}\mid\sigma^2_{gt},\eta_t\sim\mathcal N(0,e^{\eta_t}\sigma^2_{gt}).
$$

Thus the conditional standard deviation of $z$ is

$$
e^{\eta/2}.
$$

For every finite $\eta_k$, this normal distribution is converted to masses on the $z$ grid. At the exact consistency state,

$$
\eta=-\infty,
$$

all conditional probability is placed exactly at

$$
z=0.
$$

The marginal coefficient prior is therefore

$$
\pi^{(z)}_r=\sum_k\pi^{(\eta)}_kP(z_r\mid\eta_k).
$$

With the default $\eta$ prior this produces a coefficient prior containing an exact point mass contribution at zero as well as a continuous mixture over nonzero $z$ states.

For a target group $g$, define its likelihood under the explicit coefficient constraint by

$$
A_{gt}(\mu_u,z_r)=\sum_l\pi^{(\sigma)}_{gl}L_{gt}\left(\mu_u+z_r\sqrt{\sigma_l^2},\sigma_l^2\right).
$$

In exact mode, the required likelihood value is obtained by linear interpolation of the log likelihood along the $\mu$ grid separately for every $\sigma_l^2$. If

$$
\mu_u+z_r\sqrt{\sigma_l^2}
$$

lies outside the numerical $\mu$ range, its log likelihood is assigned $-\infty$.

All other groups remain represented through the variance-ratio hierarchy. The joint evidence associated with $z_r$ for target group $g$ is

$$
J_{gt}(z_r)=\sum_u\pi^{(\mu)}_uA_{gt}(\mu_u,z_r)\sum_k\pi^{(\eta)}_kP(z_r\mid\eta_k)\prod_{h\neq g}T_{ht}(\mu_u,\eta_k).
$$

The target group is excluded from the product over $h\neq g$ because it is already represented explicitly by $A_{gt}(\mu_u,z_r)$.

The stored coefficient likelihood is prior-corrected:

$$
L^{(z)}_{gt}(z_r)=\frac{J_{gt}(z_r)}{\pi^{(z)}_r},
$$

up to a multiplicative constant independent of $z_r$. Consequently,

$$
P(z_{gt}=z_r\mid y)\propto L^{(z)}_{gt}(z_r)\pi^{(z)}_r=J_{gt}(z_r).
$$

This prior correction is important because it permits the prior to be applied once for pointwise inference and once per segment, rather than once per base, when the coefficient track is spatially segmented.

A subtle implementation detail is that the common-reference coefficient inherits the integration method and `variance_floor` from the fitted variance-ratio model. Therefore, under the default configuration, common-reference coefficients use the Gaussian upper-level approximation because `VarianceRatioConfig.method="gaussian"`.

### 7.2 Standardized signal relative to zero

The second coefficient is

$$
z^{(0)}_{gt}=\frac{\mu_{gt}}{\sigma_{gt}}.
$$

For fixed $z_r$,

$$
\mu_{gt}=z_r\sqrt{\sigma^2_{gt}},
$$

and the exact likelihood is

$$
L^{(0)}_{gt}(z_r)=\sum_l\pi^{(\sigma)}_{gl}L_{gt}\left(z_r\sqrt{\sigma_l^2},\sigma_l^2\right).
$$

Again, log likelihoods are linearly interpolated on the $\mu$ grid and contributions outside the $\mu$ range are zero.

Unlike the common-reference coefficient, the zero-reference coefficient does not depend on the other groups or on $\mu_0$.

Its default prior is

$$
z^{(0)}\sim\mathcal N(0,2^2),
$$

discretized into bin masses on the $z$ grid.

The zero-reference coefficient uses the integration method specified by `CoefficientConfig`; this method is `"exact"` by default.

### 7.3 Coefficient segmentation

Both coefficient likelihoods can be supplied to the generic explicit-duration segmenter. The default coefficient segmentation uses

$$
\texttt{transition\_sd=None}
$$

and

$$
\texttt{forbid\_same\_state=True}.
$$

The appropriate coefficient prior is applied once for each inferred segment.

---

## 8. Posterior counts of footprinted or deviant groups

Posterior distributions over the number of groups satisfying a coefficient criterion are derived from the group-specific coefficient posteriors.

For the zero-reference coefficient, define

$$
p_{gt}^{\mathrm{zero}}=P\left(z^{(0)}_{gt}<-c\mid y\right).
$$

With the default threshold $c=1$,

$$
K_{\mathrm{fp,zero},t}=\sum_g\mathbf 1\left\{z^{(0)}_{gt}<-1\right\}.
$$

For common-reference deviations,

$$
p_{gt}^{\mathrm{dev}}=P\left(|z^{(\mu_0)}_{gt}|>c\mid y\right)
$$

and

$$
K_{\mathrm{fp,dev},t}=\sum_g\mathbf 1\left\{|z^{(\mu_0)}_{gt}|>c\right\}.
$$

Given group event probabilities $p_1,\ldots,p_G$, the implementation computes the Poisson-binomial distribution recursively. Initialize

$$
f_0(0)=1.
$$

After adding group $g$,

$$
f_g(k)=f_{g-1}(k)(1-p_g)+f_{g-1}(k-1)p_g,
$$

with terms outside the valid count range taken as zero.

All calculations are carried out in log space using `logaddexp`.

The resulting posterior is defined over

$$
K=0,1,\ldots,G.
$$

The calculation uses the marginal event probabilities as independent Bernoulli probabilities. Thus the implemented Poisson-binomial count distribution does not retain posterior dependence among coefficient events in different groups. This distinction is particularly relevant for the common-reference coefficients, which originate from a hierarchy involving shared $\mu_0$ and $\eta$.

Either pointwise or spatially segmented coefficient posteriors may be supplied to this calculation. In the standard workflow the segmented coefficient posterior is used.

---

## 9. Optional sample-level posterior over theta

Individual-sample inference is performed after fitting the group-level model and reuses the same $\mu$, $\sigma^2$, and $\theta$ grids.

### 9.1 Group prior-predictive distribution of theta

For group $g$, define the prior-predictive distribution

$$
\pi^{(\theta)}_{gq}=\sum_{j,k}\pi^{(\mu)}_j\pi^{(\sigma)}_{gk}K_{qjk}.
$$

This distribution is normalized over $\theta_q$ and is used as the sample-specific prior for every sample belonging to group $g$.

### 9.2 Sample-only mode

In `sample_only` mode, the stored likelihood is simply

$$
L_{it}(\theta_q)=p(y_{it}\mid\theta_q).
$$

The pointwise posterior is therefore

$$
P(\theta_{it}=\theta_q\mid y_{it})\propto p(y_{it}\mid\theta_q)\pi^{(\theta)}_{gq}.
$$

### 9.3 Group-informed leave-one-out mode

The default mode is `group_informed`.

For every sample $i$, define its likelihood marginalized over $\theta$ conditional on group-level parameters:

$$
A_{itjk}=\sum_qp(y_{it}\mid\theta_q)K_{qjk}.
$$

When constructing group information for sample $i$, that sample itself is explicitly excluded. The leave-one-out distribution over group parameters is proportional to

$$
W^{(-i)}_{tjk}\propto\pi^{(\mu)}_j\pi^{(\sigma)}_{gk}\prod_{\substack{h\in g\\h\neq i}}A_{htjk}.
$$

This is projected back to the $\theta$ grid to obtain the normalized leave-one-out posterior-predictive distribution

$$
p^{(-i)}_{gt}(\theta_q)=\sum_{j,k}K_{qjk}W^{(-i)}_{tjk}.
$$

The likelihood stored for sample $i$ is not simply the product of its own likelihood and the leave-one-out predictive distribution. Instead, a prior correction is applied:

$$
L^{\mathrm{GI}}_{it}(\theta_q)=p(y_{it}\mid\theta_q)\frac{p^{(-i)}_{gt}(\theta_q)}{\pi^{(\theta)}_{gq}}.
$$

The pointwise posterior subsequently multiplies this quantity by the group prior predictive exactly once:

$$
P(\theta_{it}=\theta_q\mid y)\propto L^{\mathrm{GI}}_{it}(\theta_q)\pi^{(\theta)}_{gq},
$$

and therefore

$$
P(\theta_{it}=\theta_q\mid y)\propto p(y_{it}\mid\theta_q)p^{(-i)}_{gt}(\theta_q).
$$

This leave-one-out construction ensures that the target sample's observation likelihood enters the posterior exactly once.

The same prior-corrected likelihood can be supplied to the generic segmenter. In that case, $\pi^{(\theta)}_g$ is applied once per segment rather than once at every base.

For memory efficiency, sample-level inference is processed in genomic chunks of eight positions by default, and the resulting sample-by-position-by-$\theta$ likelihood array is stored as `float32`. These choices change memory use but not the intended statistical model.

---

## 10. Posterior summaries and probability thresholds

All finite-grid posterior objects are normalized on their final array axis. If $p_j$ denotes the normalized mass at state $x_j$, the posterior mean is

$$
E[X\mid y]=\sum_jp_jx_j,
$$

and the MAP estimate is

$$
x_{\mathrm{MAP}}=x_{\arg\max_jp_j}.
$$

The posterior $q$-quantile is defined as the first grid value whose cumulative probability reaches or exceeds $q$:

$$
Q(q)=x_{\min\left\{j:\sum_{r\leq j}p_r\geq q\right\}}.
$$

Threshold probabilities are obtained by summing grid masses satisfying the requested strict or non-strict inequality. For example,

$$
P(X>c\mid y)=\sum_{j:x_j>c}p_j,
$$

and

$$
P(|X|>c\mid y)=\sum_{j:|x_j|>c}p_j.
$$

Thus, when a differential base is declared using a criterion such as

$$
P(|z_{gt}|>c\mid y)>0.99,
$$

both inequalities are strict in the implementation: states satisfying exactly $|z|=c$ do not contribute to the event probability.

---

## 11. Default parameter settings

Unless explicitly overridden, the calculations described above use the following numerical parameters:

| Parameter | Default |
|---|---:|
| $\mu_{\min}$ | $-6$ |
| $\mu_{\max}$ | $6$ |
| number of $\mu$ states | 481 |
| $\Delta_\mu$ | 0.025 |
| $\sigma^2_{\min}$ | $10^{-3}$ |
| $\sigma^2_{\max}$ | 10 |
| number of $\sigma^2$ states | 181 |
| nominal $\theta$ step | 0.1 |
| $\theta$ tail margin | $5\sqrt{\sigma^2_{\max}}$ |
| default number of $\theta$ states | 438 |
| $\eta_{\min}$ | $-4$ |
| $\eta_{\max}$ | 4 |
| $\Delta_\eta$ | 0.25 |
| exact $\lambda=0$ state | included |
| prior mass at $\lambda=0$ | 0.30 |
| finite-$\eta$ slab mean | $-1$ |
| finite-$\eta$ slab SD | 2 |
| variance-ratio integration | Gaussian |
| $z_{\min}$ | $-4$ |
| $z_{\max}$ | 4 |
| $\Delta_z$ | 0.25 |
| zero-reference $z$ prior SD | 2 |
| zero-reference coefficient integration | exact |
| group-mean same-state transitions | allowed |
| common-mean same-state transitions | allowed |
| ICC same-state transitions | forbidden |
| coefficient same-state transitions | forbidden |
| sample-$\theta$ same-state transitions | allowed |
| default transition attraction | none |
| sample-$\theta$ mode | group-informed |
| sample-$\theta$ storage | `float32` |

---

## 12. Computational order used for each genomic interval

For each analyzed interval, the full calculation can be reproduced in the following order.

1. Align observed counts, expected counts, sample group labels, and sample-specific dispersion models by sample identifier. At least two groups are required, and every included group must contain at least two samples.

2. Evaluate the external observation likelihood

$$
\log p(y_{it}\mid\theta_q)
$$

for every sample, position, and $\theta$ state.

3. Construct the conditional Gaussian $\theta$-bin kernel

$$
K_{qjk}=P(\theta_q\text{ bin}\mid\mu_j,\sigma_k^2).
$$

4. Within every group, estimate the empirical inverse-$\chi^2$ prior on $\sigma^2$ from the position-wise sample variances of

$$
\log(1+\mathrm{obs})-\log(1+\mathrm{exp}).
$$

5. Integrate each sample's $\theta$ likelihood through $K_{qjk}$, multiply across samples within the group, and obtain

$$
L_{gt}(\mu_j,\sigma_k^2).
$$

6. Integrate $\sigma^2$ using its group-specific prior to obtain

$$
L_{gt}(\mu_j).
$$

7. Obtain the preliminary common-mean track from

$$
\arg\max_{\mu_j}\sum_g\log L_{gt}(\mu_j)
$$

and use that track to construct the empirical normal prior on $\mu$.

8. If desired, calculate pointwise group-mean posteriors or spatially segment each group's $\mu$ likelihood.

9. Construct the population hierarchy

$$
L_t(\mu_0,\eta)=\prod_gT_{gt}(\mu_0,\eta)
$$

using either the exact finite-grid calculation or the Gaussian approximation.

10. Combine this likelihood with the $\mu_0$ and $\eta$ priors to obtain pointwise posterior distributions for the common mean and ICC.

11. For spatial inference of between-group heterogeneity, integrate $\mu_0$ pointwise, segment the resulting ICC likelihood using the explicit-duration model, and reconstruct the $\mu_0$ posterior conditional on the segmented ICC distribution.

12. For each group, construct either

$$
z^{(\mu_0)}=\frac{\mu_g-\mu_0}{\sigma_g}
$$

or

$$
z^{(0)}=\frac{\mu_g}{\sigma_g},
$$

calculate its pointwise likelihood and posterior, and optionally apply the same explicit-duration segmentation model.

13. Convert coefficient posteriors into biologically interpretable event probabilities such as

$$
P(z<-c\mid y)
$$

or

$$
P(|z|>c\mid y).
$$

14. When required, combine the group-specific event probabilities with the Poisson-binomial recursion to obtain a posterior distribution over the number of groups satisfying the event.

15. Optionally calculate sample-level $\theta$ posteriors using either the sample-only model or the default leave-one-sample-out group-informed model, and optionally segment these sample-level tracks.

This sequence maintains a deliberate distinction between likelihood information and prior information. Priors are not embedded repeatedly into intermediate likelihood arrays: each prior is applied exactly at the level at which its corresponding latent state is introduced, and prior-corrected likelihoods are used where necessary so that spatial segmentation applies a state prior once per segment rather than once per nucleotide.
