"""Two independent segmentations of group summaries x, v (groups x positions, log2).

Each base splits into a consensus value (depends on the consensus level b') and G - 1 contrasts
(depend on the group deviations delta). Consensus blocks hold b' constant; group-pattern blocks hold
delta constant. Block evidences are closed-form, sums over all partitions are exact, draws are exact.
"""

import numpy as np
from scipy.special import logsumexp

from .config import DEFAULT_MODEL, ModelConfig

OMEGA2 = np.geomspace(1e-4, 4.0, 25)


def gamma_lengths(l_max, mean, shape, l_min):
    """Log mass of a discretized gamma on lengths 0..l_max (index = length)."""
    L = np.arange(l_max + 1.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        lq = np.where(L >= l_min, (shape - 1) * np.log(L) - L * shape / mean, -np.inf)
    return lq - logsumexp(lq)


def half_cauchy(omega2, scale):
    """Log half-Cauchy mass of the grid cells around sqrt(omega2)."""
    om = np.sqrt(omega2)
    mass = np.diff(np.arctan(np.r_[0, np.sqrt(om[:-1] * om[1:]), np.inf] / scale))
    return np.log(mass / mass.sum())


def _cum(a):
    return np.concatenate([np.zeros((1,) + a.shape[1:]), np.cumsum(a, axis=0)])


class FootprintModel:
    def __init__(self, x, v, sigma2, config: ModelConfig = DEFAULT_MODEL, omega2=OMEGA2):
        x, v = np.asarray(x, float), np.asarray(v, float)
        self.G, self.n = G, n = x.shape
        if G < 2:
            raise ValueError("at least two groups are required")
        self.config, self.omega2 = config, omega2
        self.Lm = Lm = min(config.l_max, n)
        w = 1 / v.mean(axis=1)
        self.w = w = w / w.sum()
        P = np.eye(G) - np.outer(np.ones(G), w)
        self.S2 = _cum(w @ np.asarray(sigma2, float))   # within-group variance, prefix sums for block averages

        # consensus layer: m_j ~ N(b', V_j), b' ~ N(0, b_var)
        m, V = w @ x, config.tau2 + (w ** 2) @ v
        self.Pa, self.Ph, Pq = _cum(1 / V), _cum(m / V), _cum(m ** 2 / V)
        self.E_b = np.full((Lm + 1, n + 1), -np.inf)
        for L in range(1, Lm + 1):
            a = np.arange(n - L + 1)
            A, h, q = (S[a + L] - S[a] for S in (self.Pa, self.Ph, Pq))
            prec = A + 1 / config.b_var
            self.E_b[L, a] = -0.5 * q + 0.5 * h ** 2 / prec - 0.5 * np.log(config.b_var * prec)

        # group-pattern layer: c_j ~ N(delta, C_j), delta ~ N(0, omega^2 K); evidence kept per omega
        c = (P @ x)[:-1]
        C = np.einsum("gh,hj,lh->jgl", P, v, P)[:, :-1, :-1]
        Ci = np.linalg.inv(C)
        self.PA, self.PH = _cum(Ci), _cum(np.einsum("jgh,hj->jg", Ci, c))
        self.PQ = _cum(np.einsum("gj,jgh,hj->j", c, Ci, c))
        self.Linv = np.linalg.inv(np.linalg.cholesky(np.linalg.inv((P @ P.T)[:-1, :-1])))
        self.Q = np.vstack([np.eye(G - 1), -w[:-1] / w[-1]])
        self.E_d_omega = np.full((Lm + 1, n + 1, omega2.size), -np.inf)
        for L in range(1, Lm + 1):
            a = np.arange(n - L + 1)
            self.E_d_omega[L, a] = self._pattern(a, L)[3]

        self.const = (-0.5 * np.log(2 * np.pi * V).sum() - 0.5 * np.linalg.slogdet(2 * np.pi * C)[1].sum()
                      + n * np.log(abs(np.linalg.det(np.vstack([w, P[:-1]])))))

    def _pattern(self, a, L):
        """Eigen-quantities and per-omega log evidence of group-pattern blocks [a, a + L)."""
        A, h, q = (S[a + L] - S[a] for S in (self.PA, self.PH, self.PQ))
        lam, vec = np.linalg.eigh(self.Linv @ A @ self.Linv.T)
        lam = np.maximum(lam, 0.0)
        U = np.einsum("ij,sjk->sik", self.Linv.T, vec)
        g = np.einsum("sik,si->sk", U, h)
        w2 = self.omega2[:, None]
        lm = (-0.5 * q[:, None] + 0.5 * np.sum(g[:, None] ** 2 / (lam[:, None] + 1 / w2), axis=2)
              - 0.5 * np.sum(np.log1p(lam[:, None] * w2), axis=2))
        return lam, U, g, lm

    def layers(self, config: ModelConfig | None = None):
        """Evidence tables and length priors of both layers for the given constants (lengths, omega_scale)."""
        cfg = config or self.config
        lq_b = gamma_lengths(self.Lm, cfg.cons_mean, cfg.cons_shape, cfg.cons_min)
        lq_d = gamma_lengths(self.Lm, cfg.pattern_mean, cfg.pattern_shape, cfg.pattern_min)
        lw = half_cauchy(self.omega2, cfg.omega_scale)
        return (self.E_b, lq_b), (logsumexp(self.E_d_omega + lw, axis=2), lq_d), lw

    def log_evidence(self, config: ModelConfig | None = None):
        """(consensus, group-pattern) log evidence; their sum plus self.const is log p(x)."""
        (E_b, lq_b), (E_d, lq_d), _ = self.layers(config)
        return _forward(E_b, lq_b)[-1], _forward(E_d, lq_d)[-1]

    def sample(self, n_draws=1000, seed=0):
        """Exact draws. Returns depth (draws x n x G), P(group-pattern boundary) and P(consensus boundary)
        between bases j and j + 1, and the posterior mean ICC per base (block-averaged within-group variance)."""
        rng = np.random.default_rng(seed)
        n, G = self.n, self.G
        (E_b, lq_b), (E_d, lq_d), lw = self.layers()
        b, marks_b = np.zeros((n_draws, n)), np.zeros((n_draws, n), bool)
        for t, a, L in _backward(E_b, lq_b, _forward(E_b, lq_b), n_draws, rng):
            prec = self.Pa[a + L] - self.Pa[a] + 1 / self.config.b_var
            b[t, a] = (self.Ph[a + L] - self.Ph[a]) / prec + rng.standard_normal(a.size) / np.sqrt(prec)
            marks_b[t, a] = True
        delta, icc, marks_d = np.zeros((n_draws, n, G)), np.zeros((n_draws, n)), np.zeros((n_draws, n), bool)
        for t, a, L in _backward(E_d, lq_d, _forward(E_d, lq_d), n_draws, rng):
            delta[t, a], icc[t, a] = self._draw_pattern(a, L, lw, rng)
            marks_d[t, a] = True
        fb, fd = _fill_index(marks_b), _fill_index(marks_d)
        depth = np.take_along_axis(b, fb, 1)[..., None] + np.take_along_axis(delta, fd[..., None], 1)
        return (depth.astype(np.float32), marks_d[:, 1:].mean(0), marks_b[:, 1:].mean(0),
                np.take_along_axis(icc, fd, 1).mean(0))

    def _draw_pattern(self, a, L, lw, rng):
        blocks, inv = np.unique(np.stack([a, L], 1), axis=0, return_inverse=True)
        lam, U, g, lm = (z[inv.ravel()] for z in self._pattern(blocks[:, 0], blocks[:, 1]))
        k = _pick(lm + lw, rng)
        dd = 1 / (lam + 1 / self.omega2[k, None])
        z = g * dd + rng.standard_normal(dd.shape) * np.sqrt(dd)
        om2 = self.omega2[k]
        return np.einsum("gi,sik,sk->sg", self.Q, U, z), om2 / (om2 + (self.S2[a + L] - self.S2[a]) / L)


def _forward(E, lq):
    """Tool B: log Z(e) summed over all partitions of [0, e)."""
    n, Ls = E.shape[1] - 1, np.arange(1, E.shape[0])
    Z = np.full(n + 1, -np.inf)
    Z[0] = 0.0
    for e in range(1, n + 1):
        L = Ls[Ls <= e]
        Z[e] = logsumexp(Z[e - L] + lq[L] + E[L, e - L])
    return Z


def _backward(E, lq, Z, n_draws, rng):
    """Exact partition draws, vectorized over draws; yields (draw indices, block starts, block lengths)."""
    Ls, c = np.arange(1, E.shape[0]), np.full(n_draws, E.shape[1] - 1)
    while np.any(c > 0):
        t = np.flatnonzero(c > 0)
        a_ = np.maximum(c[t, None] - Ls, 0)
        W = np.where(Ls <= c[t, None], Z[a_] + lq[Ls] + E[Ls, a_], -np.inf)
        L = Ls[_pick(W, rng)]
        yield t, c[t] - L, L
        c[t] -= L


def _fill_index(marks):
    return np.maximum.accumulate(np.where(marks, np.arange(marks.shape[1]), 0), axis=1)


def _pick(log_w, rng):
    """One categorical draw per row (Gumbel-max)."""
    return np.argmax(log_w + rng.gumbel(size=log_w.shape), axis=1)
