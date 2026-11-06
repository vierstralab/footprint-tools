"""Stage 2, per DHS: group tracks, discriminability and footprint probability from the samples' posteriors (m, v).

pi_gh = mean over s in g, t in h of Phi((m_t - m_s) / sqrt(v_s + v_t)): P(a sample of g is more depleted than one of h).
A_g = mean over h != g of pi_gh (more depleted) or of 1 - pi_gh (less depleted); 95% lower bound by resampling
samples within groups. Delta_g = mean depth of g - median over h != g of h's mean depth.
p_any = P(min_g mean depth < -fold) = 1 - prod_g Phi((mean_g + fold) / se_g).
"""

import numpy as np
from scipy.special import ndtr


def group_tracks(m: np.ndarray, v: np.ndarray, idx: list[np.ndarray]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Group mean depth, its standard error, and Delta (groups x positions)."""
    mean = np.stack([m[i].mean(0) for i in idx])
    se = np.stack([np.sqrt((m[i].var(0) + v[i].mean(0)) / len(i)) for i in idx])
    delta = np.stack([mean[g] - np.median(np.delete(mean, g, 0), 0) for g in range(len(idx))])
    return mean, se, delta


def discriminability(m: np.ndarray, v: np.ndarray, idx: list[np.ndarray], n_boot: int, seed: int) -> dict[str, np.ndarray]:
    """A for "more" and "less" depleted, and their 95% bootstrap lower bounds (groups x positions)."""
    G, n = len(idx), m.shape[1]; rng = np.random.default_rng(seed)
    w = [rng.multinomial(len(i), np.full(len(i), 1 / len(i)), n_boot) / len(i) for i in idx]   # resampling weights
    a, boot = np.zeros((G, n)), np.zeros((G, n_boot, n))       # sums of pi_gh over the other groups, and per replicate
    for g in range(G):
        for h in range(g + 1, G):
            ig, ih = idx[g], idx[h]
            P = ndtr((m[ih][None] - m[ig][:, None]) / np.sqrt(v[ig][:, None] + v[ih][None]))      # n_g x n_h x n
            p = P.mean((0, 1))
            b = np.einsum('rtj,rt->rj', (w[g] @ P.reshape(len(ig), -1)).reshape(n_boot, len(ih), n), w[h])
            a[g] += p; a[h] += 1 - p; boot[g] += b; boot[h] += 1 - b
    a, boot = a / (G - 1), boot / (G - 1)
    return dict(a_more=a, a_less=1 - a, lo_more=np.quantile(boot, .05, axis=1), lo_less=np.quantile(1 - boot, .05, axis=1))


def p_footprint(mean: np.ndarray, se: np.ndarray, fold: float) -> np.ndarray:
    """P(at least one group's mean depth < -fold), groups treated as independent normals."""
    return 1 - np.prod(ndtr((mean + fold) / se), 0)


def runs(mask: np.ndarray, min_len: int = 1) -> list[tuple[int, int]]:
    """Half-open (start, end) runs of True of at least min_len in a 1-D boolean array."""
    edges = np.flatnonzero(np.diff(np.r_[0, np.asarray(mask, int), 0]))
    return [(s, e) for s, e in zip(edges[::2].tolist(), edges[1::2].tolist()) if e - s >= min_len]
