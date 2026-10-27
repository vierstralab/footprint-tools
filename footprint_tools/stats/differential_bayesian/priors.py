"""Length and scale priors of the nested model; all returned as log masses."""

import numpy as np
from scipy.special import logsumexp


def gamma_length_log_mass(l_max: int, mean: float, shape: float, l_min: int) -> np.ndarray:
    """Discretized gamma on lengths 0..l_max (index = length), zero below l_min."""
    L = np.arange(l_max + 1, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        lq = (shape - 1) * np.log(L) - L * shape / mean
    lq[L < l_min] = -np.inf
    return lq - logsumexp(lq)


def half_cauchy_log_mass(omega2: np.ndarray, scale: float) -> np.ndarray:
    """Half-Cauchy(scale) mass on omega of the grid cells around sqrt(omega2)."""
    om = np.sqrt(omega2)
    edges = np.r_[0.0, np.sqrt(om[:-1] * om[1:]), np.inf]
    mass = np.diff(2 / np.pi * np.arctan(edges / scale))
    return np.log(np.maximum(mass, 1e-300) / mass.sum())


def tiling_log_weight(lq: np.ndarray) -> np.ndarray:
    """u(L) = summed prior weight of all tilings of length L by blocks with weights q."""
    lu = np.full(lq.size, -np.inf)
    lu[0] = 0.0
    for d in range(1, lq.size):
        lu[d] = logsumexp(lu[:d] + lq[d - np.arange(d)])
    return lu


def outer_length_log_weight(lq_inner: np.ndarray, lq_outer: np.ndarray) -> np.ndarray:
    """Outer-block weight q_out(L) / u(L): normalizes the inner prior within each outer block."""
    lu = tiling_log_weight(lq_inner)
    out = np.full(lq_outer.size, -np.inf)
    ok = np.isfinite(lq_outer) & np.isfinite(lu)
    out[ok] = lq_outer[ok] - lu[ok]
    return out
