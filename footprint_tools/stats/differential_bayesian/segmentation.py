"""Stage 1: each sample segmented on its own, without group labels.

y_sj ~ NB(m_s(lambda_sj 2^theta_sj), r_s(lambda_sj 2^theta_sj)) with the sample's dispersion model; theta is piecewise
constant, block lengths ~ discretized gamma, block levels ~ N(0, level_sd^2) on a grid. Exact forward-backward sums over
all segmentations give each base's posterior mean and variance of theta.
"""

from collections.abc import Sequence

import numpy as np
from scipy.special import gammaln, logsumexp

from .config import DEFAULT_DF_CONFIG, DifferentialFootprintConfig


def segment_samples(obs: np.ndarray, exp: np.ndarray, disp_models: Sequence,
                    config: DifferentialFootprintConfig = DEFAULT_DF_CONFIG) -> tuple[np.ndarray, np.ndarray]:
    """obs, exp: samples x positions; disp_models: one per sample. Returns posterior mean and variance (samples x positions)."""
    grid = config.grid()
    level = -0.5 * (grid / config.level_sd) ** 2
    lengths = gamma_lengths(min(config.length_max, obs.shape[1]), config.length_mean, config.length_shape, config.length_min)
    m, v = np.empty(obs.shape), np.empty(obs.shape)
    for lo in range(0, len(obs), config.chunk):
        sl = slice(lo, lo + config.chunk)
        m[sl], v[sl] = forward_backward(nb_loglik(disp_models[sl], obs[sl], exp[sl], grid), grid, level - logsumexp(level), lengths)
    return m, v


def nb_loglik(models: Sequence, obs: np.ndarray, exp: np.ndarray, grid: np.ndarray) -> np.ndarray:
    """log NB(y | m_s(e), r_s(e)) at e = exp * 2^grid, samples x positions x grid. Below its first breakpoint the mean
    is extended through the origin (slope of the first segment) and 1/r is held at its first-breakpoint value."""
    mp, rp = (np.stack([getattr(m, a) for m in models]) for a in ("mu_params", "r_params"))
    e = exp[:, :, None] * np.exp2(grid)
    first = lambda p, i: p[:, i, None, None]
    m = np.maximum(np.where(e < first(mp, 0), first(mp, 6) * e, _pwl(e, mp, 3)), 1e-4)
    r = 1 / np.maximum(_pwl(np.maximum(e, first(rp, 0)), rp, 5), 0.01)
    y = obs[:, :, None]
    ll = gammaln(y + r) - gammaln(r) - gammaln(y + 1) + r * np.log(r / (r + m)) + y * np.log(m / (r + m))
    return np.where(exp[:, :, None] > 0, ll, 0.0)


def _pwl(x: np.ndarray, p: np.ndarray, k: int) -> np.ndarray:
    """footprint-tools piecewise-linear curve with k segments per sample: breaks p[:, :k-1], intercepts p[:, k:2k],
    slopes p[:, 2k:3k]; x is samples x ...."""
    shape = (len(p),) + (1,) * (x.ndim - 1)
    seg = sum((x >= p[:, i].reshape(shape)).astype(int) for i in range(k - 1)).reshape(len(p), -1)
    take = lambda c: np.take_along_axis(c, seg, 1).reshape(x.shape)
    return take(p[:, k:2 * k]) + take(p[:, 2 * k:3 * k]) * x


def gamma_lengths(l_max: int, mean: float, shape: float, l_min: int) -> np.ndarray:
    """Log mass of a discretized gamma on lengths 0..l_max (index = length)."""
    L = np.arange(l_max + 1.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        lq = np.where(L >= l_min, (shape - 1) * np.log(L) - L * shape / mean, -np.inf)
    return lq - logsumexp(lq)


def forward_backward(ll: np.ndarray, grid: np.ndarray, level: np.ndarray, lengths: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Exact posterior moments of a piecewise-constant theta. ll: samples x positions x grid (log-likelihoods);
    level: log prior mass of a block level on the grid; lengths: log prior of block lengths (index = length)."""
    S, n, T = ll.shape
    Ls = np.arange(1, len(lengths))
    C = np.concatenate([np.zeros((S, 1, T)), np.cumsum(ll, 1)], 1)
    E = np.full((S, len(lengths), n + 1), -np.inf)                       # block evidence by length and start
    M1, M2 = np.zeros(E.shape), np.zeros(E.shape)                       # block level moments
    for L in Ls:
        blk = C[:, L:] - C[:, :-L] + level
        z = logsumexp(blk, 2); p = np.exp(blk - z[..., None])
        E[:, L, :n - L + 1], M1[:, L, :n - L + 1], M2[:, L, :n - L + 1] = z, p @ grid, p @ grid ** 2
    Zf, Zb = np.full((S, n + 1), -np.inf), np.full((S, n + 1), -np.inf)
    Zf[:, 0], Zb[:, n] = 0, 0
    for e in range(1, n + 1):
        L = Ls[Ls <= e]; Zf[:, e] = logsumexp(Zf[:, e - L] + lengths[L] + E[:, L, e - L], 1)
    for a in range(n - 1, -1, -1):
        L = Ls[Ls <= n - a]; Zb[:, a] = logsumexp(lengths[L] + E[:, L, a] + Zb[:, a + L], 1)
    D1, D2 = np.zeros((S, n + 1)), np.zeros((S, n + 1))                 # difference arrays over covered bases
    for L in Ls:
        a = np.arange(n - L + 1)
        w = np.exp(Zf[:, a] + lengths[L] + E[:, L, a] + Zb[:, a + L] - Zf[:, n, None])
        for D, M in ((D1, M1), (D2, M2)):
            D[:, a] += w * M[:, L, a]; D[:, a + L] -= w * M[:, L, a]
    m, m2 = np.cumsum(D1, 1)[:, :n], np.cumsum(D2, 1)[:, :n]
    return m, np.maximum(m2 - m ** 2, 1e-4)
