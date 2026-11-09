"""Stage 2, per DHS: p-values for differential occupancy from the stage-1 posteriors (m, v; samples x positions).

For two samples s, t the posterior depths differ by D_st = theta_t - theta_s ~ N(m_t - m_s, v_s + v_t): the difference
m_t - m_s and the distance E|D_st|. Effect of group g against group h, from all pairs s in g, t in h:
    d_gh = mean(m_t - m_s) / mean E|D_st| = (mean_h m - mean_g m) / mean E|D_st|  in [-1, 1],
d_gh > 0: g more depleted; d_hg = -d_gh. A magnitude-weighted Cliff's delta (with exact depths, the distance-weighted
mean sign of the pairwise differences); equal to Dagum's relative economic affluence. Effect of a set S of groups
against the rest: d_S = mean of d_gh over g in S, h not in S.
Statistic per position: T = max over every set of up to k_max groups and both directions of (|d_S| - e) / SD_S,
SD_S from the grouped delete-one jackknife. Reference: T (e = 0) under label permutations, one permutation for all
positions; beyond the permutations, a generalised Pareto tail fitted to the top values.
"""

import itertools

import numpy as np
from numpy.lib.stride_tricks import sliding_window_view
from scipy.special import ndtr


def set_weights(G: int, k_max: int) -> tuple[list[tuple[int, ...]], np.ndarray]:
    """All sets S of up to k_max groups, with d_S = a_S . R for the row sums R_g = sum_h d_gh: summing R over S covers
    the k(G - k) pairs across S and the pairs inside S, which cancel (d_gh + d_hg = 0); dividing by the number of
    cross pairs gives the mean. Sets of up to k_max groups and both signs cover every split when k_max >= G // 2,
    since d of the complement is -d_S."""
    sets = [s for k in range(1, k_max + 1) for s in itertools.combinations(range(G), k)]
    in_set = np.zeros((len(sets), G))
    for i, s in enumerate(sets):
        in_set[i, list(s)] = 1
    k = in_set.sum(1)
    return sets, in_set / (k * (G - k))[:, None]


def expected_abs_diff(m: np.ndarray, v: np.ndarray) -> np.ndarray:
    """E|theta_t - theta_s| for every pair of samples (positions x N x N), theta_s ~ N(m_s, v_s) independent."""
    mu = m.T[:, None, :] - m.T[:, :, None]
    sd = np.sqrt(v.T[:, :, None] + v.T[:, None, :])
    z = mu / sd
    return (mu * (2 * ndtr(z) - 1) + 2 * sd * np.exp(-z ** 2 / 2) / np.sqrt(2 * np.pi)).astype(np.float32)


def group_effects(row_sums: np.ndarray, m: np.ndarray, g: np.ndarray, G: int) -> tuple[np.ndarray, np.ndarray]:
    """d_gh (G x G x positions) and, for each sample s, the change in d_{g(s) h} when s is left out (N x G x positions).
    row_sums: sum of E|D_st| over t in h (N x G x positions); m: posterior means (N x positions)."""
    N, n = m.shape
    onehot = np.eye(G)[g]
    size = onehot.sum(0)[:, None]                                             # n_g
    total = onehot.T @ m                                                      # sum of m over each group
    abs_diff = (onehot.T @ row_sums.reshape(N, -1)).reshape(G, G, n)          # sum of E|D_st|, s in g, t in h
    mean_diff = size[:, None] * total[None] - size[None] * total[:, None]     # n_g n_h (mean_h - mean_g)
    off_diag = ~np.eye(G, dtype=bool)[..., None]
    d = np.where(off_diag, mean_diff / np.where(off_diag, abs_diff, 1), 0)
    # s left out of its group g: n_g - 1 samples, total_g - m_s, and its own pairs removed from abs_diff
    mean_diff_loo = (size[g] - 1)[:, None] * total[None] - size[None] * (total[g] - m)[:, None]
    with np.errstate(divide='ignore', invalid='ignore'):
        dd = mean_diff_loo / (abs_diff[g] - row_sums) - d[g]
    dd[np.arange(N), g] = 0
    return d, dd


def jackknife_cov(dd: np.ndarray, g: np.ndarray, G: int) -> np.ndarray:
    """Grouped delete-one jackknife covariance of the row sums R (G x G x positions)."""
    N = len(g)
    dR = -dd                                      # R_h changes through d_{h g(s)} = -d_{g(s) h}
    dR[np.arange(N), g] = dd.sum(1)               # R_{g(s)} through all of its own row
    for h in range(G):
        i = g == h
        n_h = i.sum()
        dR[i] = np.sqrt((n_h - 1) / n_h) * (dR[i] - dR[i].mean(0))
    return np.einsum('sgj,shj->ghj', dR, dR)


def set_stats(delta: np.ndarray, m: np.ndarray, labelings: list[np.ndarray], G: int, a: np.ndarray):
    """For each labeling (group index per sample): d_S and its jackknife SD for every set (sets x positions).
    delta: E|D_st| (positions x N x N); m: posterior means (N x positions)."""
    n, N = delta.shape[:2]
    Z = np.concatenate([np.eye(G, dtype=np.float32)[g] for g in labelings], 1)                     # N x (L G)
    row_sums = (delta.reshape(-1, N) @ Z).reshape(n, N, len(labelings), G).transpose(2, 1, 3, 0)   # L x N x G x n
    for g, rs in zip(labelings, row_sums):
        d, dd = group_effects(rs, m, g, G)
        C = jackknife_cov(dd, g, G)
        d_S = a @ d.sum(1)
        var = ((a @ C.reshape(G, -1)).reshape(len(a), G, n) * a[:, :, None]).sum(1)                # a_S' C a_S
        yield d_S, np.sqrt(np.maximum(var, 1e-12))


def log10_p(M: np.ndarray, T: np.ndarray, tail_frac: float, min_hits: int) -> tuple[np.ndarray, np.ndarray]:
    """log10 P(M >= T) per column of M (replicates x columns): (1 + hits) / (B + 1) when at least min_hits replicates
    reach T, else a generalised Pareto tail fitted to the top tail_frac by probability-weighted moments (Knijnenburg et al. 2009).
    xi < 0 is taken as 0. Also returns the fitted xi."""
    B = len(M)
    n_tail = max(int(tail_frac * B), 20)
    M_sorted = np.sort(M, 0)
    u = M_sorted[-n_tail - 1].astype(float)                                   # tail threshold
    y = M_sorted[-n_tail:] - u                                                # exceedances, ascending
    i = np.arange(n_tail)[:, None]
    a0, a1 = y.mean(0), (y * (n_tail - 1 - i) / (n_tail - 1)).mean(0)         # probability-weighted moments
    xi_fit = 2 - a0 / np.maximum(a0 - 2 * a1, 1e-12)
    xi = np.clip(xi_fit, 0, None)
    scale = np.maximum(np.where(xi > 0, 2 * a0 * a1 / np.maximum(a0 - 2 * a1, 1e-12), a0), 1e-9)
    x = np.maximum(T - u, 0)
    with np.errstate(divide='ignore', invalid='ignore'):
        log_sf = np.where(xi > 0, -np.log1p(xi * x / scale) / np.where(xi > 0, xi, 1), -x / scale)
    hits = (M >= T).sum(0)
    empirical = np.log10((1 + hits) / (B + 1))
    tail = np.log10(n_tail / B) + log_sf / np.log(10)
    return np.where(hits >= min_hits, empirical, tail), xi_fit


def differential_pvalues(m: np.ndarray, v: np.ndarray, g: np.ndarray, G: int, config) -> dict[str, np.ndarray]:
    """g: group index per sample. Per position, for the call threshold (config.d_min) and each direction: the best set,
    its d on the called side and log10 p. For every reporting threshold: log10 p of the best hypothesis per position and
    per DHS (best window of config.min_len positions)."""
    _, a = set_weights(G, config.k_max)
    delta = expected_abs_diff(m, v)
    stats = lambda labelings: set_stats(delta, m, labelings, G, a)
    (d, SD), = stats([g])

    rng, M = np.random.default_rng(config.seed), []                           # largest T (e = 0) per permutation
    for b in range(0, config.n_perm, config.perm_batch):
        perms = [rng.permutation(g) for _ in range(min(config.perm_batch, config.n_perm - b))]
        M += [(np.abs(x) / sd).max(0) for x, sd in stats(perms)]
    M = np.stack(M)
    best_window = lambda x: sliding_window_view(x, config.min_len, axis=-1).min(-1).max(-1)
    M_dhs = best_window(M)[:, None]
    log_p = lambda null, T: log10_p(null, T, config.tail_frac, config.min_hits)[0]

    cols, out = np.arange(m.shape[1]), {}
    for kind, sign in (('more', 1), ('less', -1)):
        t = (sign * d - config.d_min) / SD
        best = t.argmax(0)
        out[f'lp_{kind}'] = log_p(M, t[best, cols])
        out[f'set_{kind}'] = best
        out[f'd_{kind}'] = sign * d[best, cols]

    lp, lp_dhs = [], []
    for e in config.report_d:
        T = ((np.abs(d) - e) / SD).max(0)
        lp.append(log_p(M, T))
        lp_dhs.append(log_p(M_dhs, np.array([best_window(T)]))[0])
    out.update(lp=np.stack(lp), lp_dhs=np.array(lp_dhs), xi=log10_p(M, M[0], config.tail_frac, config.min_hits)[1],
               c95=np.quantile(M, .95, 0))
    return out
