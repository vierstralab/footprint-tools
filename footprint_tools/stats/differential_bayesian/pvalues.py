"""Stage 2, per DHS: p-values for differential occupancy from the stage-1 posteriors (m, v; samples x positions).

Effect of a set S of groups against the rest (W_S > 1/2: S more depleted):
    W_gh = sum_{s in g, t in h} E[(m_t - m_s)_+] / sum E|m_t - m_s|, the two depths independent normals;
    W_S = mean of W_gh over g in S, h not in S = (sum_{g in S} R_g - k(k-1)/2) / (k (G - k)), R_g = sum_h W_gh.
Statistic per base: T = max over every set of up to k_max groups and both directions of (|W_S - 1/2| - e) / SD_S,
SD_S from the grouped delete-one jackknife. Reference: T (e = 0) under label permutations, one permutation for all
bases; beyond the permutations, a generalised Pareto tail fitted to the top values.
"""

import itertools

import numpy as np
from numpy.lib.stride_tricks import sliding_window_view
from scipy.special import ndtr


def set_weights(G: int, k_max: int) -> tuple[list[tuple[int, ...]], np.ndarray, np.ndarray]:
    """All sets of up to k_max groups and the weights a_S, offsets c_S with W_S = a_S . R - c_S."""
    sets = [s for k in range(1, k_max + 1) for s in itertools.combinations(range(G), k)]
    ind = np.zeros((len(sets), G))
    for i, s in enumerate(sets):
        ind[i, list(s)] = 1
    k = ind.sum(1)
    return sets, ind / (k * (G - k))[:, None], (k - 1) / (2 * (G - k))


def pair_terms(m: np.ndarray, v: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """E[(m_t - m_s)_+] and E|m_t - m_s| for every pair of samples, positions x N x N."""
    mu = m.T[:, None, :] - m.T[:, :, None]
    sd = np.sqrt(v.T[:, :, None] + v.T[:, None, :]); z = mu / sd
    pos = mu * ndtr(z) + sd * np.exp(-z ** 2 / 2) / np.sqrt(2 * np.pi)
    return pos.astype(np.float32), (2 * pos - mu).astype(np.float32)


def set_stats(pos: np.ndarray, ab: np.ndarray, labelings: list[np.ndarray], G: int, a: np.ndarray, c: np.ndarray):
    """For each labeling (group index per sample): W_S and its jackknife SD for every set (sets x positions)."""
    n, N = pos.shape[:2]
    Z = np.concatenate([np.eye(G, dtype=np.float32)[g] for g in labelings], 1)                  # N x (L G)
    PZ, AZ = [(x.reshape(-1, N) @ Z).reshape(n, N, len(labelings), G).transpose(2, 1, 3, 0) for x in (pos, ab)]
    off_diag = ~np.eye(G, dtype=bool)[..., None]
    for g, pz, az in zip(labelings, PZ, AZ):                                                     # pz: N x G x n
        onehot = np.eye(G)[g]
        num, den = [(onehot.T @ x.reshape(N, -1)).reshape(G, G, n) for x in (pz, az)]
        W = np.where(off_diag, num / np.where(off_diag, den, 1), 0)
        dW = (num[g] - pz) / (den[g] - az) - W[g]                     # change in W_{g(s) h} when sample s is left out
        dW[np.arange(N), g] = 0
        D = -dW; D[np.arange(N), g] = dW.sum(1)                       # change in every R_h
        for h in range(G):
            i = g == h; D[i] = np.sqrt((i.sum() - 1) / i.sum()) * (D[i] - D[i].mean(0))
        C = np.einsum('ngj,nhj->ghj', D, D)                           # jackknife covariance of R, G x G x n
        var = ((a @ C.reshape(G, -1)).reshape(len(a), G, n) * a[:, :, None]).sum(1)
        yield a @ W.sum(1) - c[:, None], np.sqrt(np.maximum(var, 1e-12))


def log10_p(M: np.ndarray, T: np.ndarray, tail_frac: float, min_hits: int) -> tuple[np.ndarray, np.ndarray]:
    """log10 P(M >= T) per column of M (replicates x columns): (1 + hits) / (B + 1) when at least min_hits replicates
    reach T, else a generalised Pareto tail fitted to the top tail_frac by probability-weighted moments, xi < 0
    taken as 0. Also returns the fitted xi."""
    B = len(M); n_exc = max(int(tail_frac * B), 20)
    s = np.sort(M, 0); u = s[-n_exc - 1].astype(float); y = s[-n_exc:] - u
    i = np.arange(n_exc)[:, None]
    a0, a1 = y.mean(0), (y * (n_exc - 1 - i) / (n_exc - 1)).mean(0)
    xi_fit = 2 - a0 / np.maximum(a0 - 2 * a1, 1e-12)
    xi = np.clip(xi_fit, 0, None)
    scale = np.maximum(np.where(xi > 0, 2 * a0 * a1 / np.maximum(a0 - 2 * a1, 1e-12), a0), 1e-9)
    x = np.maximum(T - u, 0)
    with np.errstate(divide='ignore', invalid='ignore'):
        tail = np.where(xi > 0, -np.log1p(xi * x / scale) / np.where(xi > 0, xi, 1), -x / scale) / np.log(10)
    hits = (M >= T).sum(0)
    return np.where(hits >= min_hits, np.log10((1 + hits) / (B + 1)), np.log10(n_exc / B) + tail), xi_fit


def differential_pvalues(m: np.ndarray, v: np.ndarray, g: np.ndarray, G: int, config) -> dict[str, np.ndarray]:
    """g: group index per sample. Per position, for the call threshold (config.w_min) and each direction: the best set,
    its W on the called side and log10 p. For every reporting threshold: log10 p of the best hypothesis per position and
    per DHS (best window of config.min_len bases)."""
    _, a, c = set_weights(G, config.k_max)
    pos, ab = pair_terms(m, v)
    (W, SD), = set_stats(pos, ab, [g], G, a, c)
    rng, M = np.random.default_rng(config.seed), []
    for b in range(0, config.n_perm, config.perm_batch):
        perms = [rng.permutation(g) for _ in range(min(config.perm_batch, config.n_perm - b))]
        M += [(np.abs(w - .5) / sd).max(0) for w, sd in set_stats(pos, ab, perms, G, a, c)]
    M = np.stack(M); Mw = sliding_window_view(M, config.min_len, axis=1).min(2).max(1, keepdims=True)
    p = lambda T: log10_p(M, T, config.tail_frac, config.min_hits)[0]
    cols, out = np.arange(m.shape[1]), {}
    for kind, sign in (('more', 1), ('less', -1)):
        t = (sign * (W - .5) - (config.w_min - .5)) / SD; i = t.argmax(0)
        out[f'lp_{kind}'], out[f'set_{kind}'], out[f'w_{kind}'] = p(t[i, cols]), i, .5 + sign * (W[i, cols] - .5)
    lp, lp_dhs = [], []
    for w in config.report_w:
        T = ((np.abs(W - .5) - (w - .5)) / SD).max(0)
        lp.append(p(T))
        lp_dhs.append(log10_p(Mw, np.array([sliding_window_view(T, config.min_len).min(1).max()]), config.tail_frac, config.min_hits)[0][0])
    out.update(lp=np.stack(lp), lp_dhs=np.array(lp_dhs), xi=log10_p(M, M[0], config.tail_frac, config.min_hits)[1],
               c95=np.quantile(M, .95, 0))
    return out
