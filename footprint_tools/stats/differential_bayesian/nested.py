"""Nested product-partition model for group summaries x, v (groups x positions, log2).

Outer blocks hold the group contrasts constant; inner blocks, nested inside them, hold
the consensus level constant. Block evidences are closed-form, the sum over all nested
segmentations is exact, and posterior draws are exact.
"""

import numpy as np
from scipy.special import logsumexp

from .config import DEFAULT_NESTED, NestedConfig

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


def outer_weights(lq, lq_out):
    """q_out(L) / u(L), u(L) = prior weight of all inner tilings of L: the inner prior is normalized per outer block."""
    lu = np.zeros(lq.size)
    for d in range(1, lq.size):
        lu[d] = logsumexp(lu[:d] + lq[d:0:-1])
    with np.errstate(invalid="ignore"):
        return np.where(np.isfinite(lq_out), lq_out - lu, -np.inf)


def _cum(a):
    return np.concatenate([np.zeros((1,) + a.shape[1:]), np.cumsum(a, axis=0)])


class NestedModel:
    def __init__(self, x, v, config: NestedConfig = DEFAULT_NESTED, omega2=OMEGA2, log_omega_prior=None):
        x, ve = np.asarray(x, float), config.kappa * np.asarray(v, float)
        self.G, self.n = G, n = x.shape
        if G < 2:
            raise ValueError("at least two groups are required")
        self.Lm = Lm = min(config.l_max, n)
        self.omega2 = omega2
        self.log_omega_prior = half_cauchy(omega2, config.omega_scale) if log_omega_prior is None else log_omega_prior
        self.lq = gamma_lengths(Lm, config.inner_mean, config.len_shape, config.len_min)
        self.lo = outer_weights(self.lq, gamma_lengths(Lm, config.outer_mean, config.len_shape, config.len_min))
        self.icc_grid = omega2 / (omega2 + config.sigma2_bar)
        self.b_var = config.b_var

        w = 1 / ve.mean(axis=1)
        self.w = w = w / w.sum()
        P = np.eye(G) - np.outer(np.ones(G), w)

        # consensus layer: m_j ~ N(b', V_j)
        m, V = w @ x, config.tau2 + (w ** 2) @ ve
        self.Pa, self.Ph, Pq = _cum(1 / V), _cum(m / V), _cum(m ** 2 / V)
        self.mb = np.full((Lm + 1, n + Lm + 1), -np.inf)
        for L in range(1, Lm + 1):
            a = np.arange(n - L + 1)
            A, h, q = (S[a + L] - S[a] for S in (self.Pa, self.Ph, Pq))
            prec = A + 1 / self.b_var
            self.mb[L, a] = -0.5 * q + 0.5 * h ** 2 / prec - 0.5 * np.log(self.b_var * prec)

        # contrast layer: c_j ~ N(delta, C_j), delta ~ N(0, omega^2 K)
        c = (P @ x)[:-1]
        C = np.einsum("gh,hj,lh->jgl", P, ve, P)[:, :-1, :-1]
        Ci = np.linalg.inv(C)
        self.PA, self.PH = _cum(Ci), _cum(np.einsum("jgh,hj->jg", Ci, c))
        self.PQ = _cum(np.einsum("gj,jgh,hj->j", c, Ci, c))
        self.Linv = np.linalg.inv(np.linalg.cholesky(np.linalg.inv((P @ P.T)[:-1, :-1])))
        self.Q = np.vstack([np.eye(G - 1), -w[:-1] / w[-1]])
        self.mdp = np.full((Lm + 1, n + 1), -np.inf)
        for L in range(1, Lm + 1):
            a = np.arange(n - L + 1)
            self.mdp[L, a] = logsumexp(self._outer(a, L)[3], axis=1)

        self.const = (-0.5 * np.log(2 * np.pi * V).sum() - 0.5 * np.linalg.slogdet(2 * np.pi * C)[1].sum()
                      + n * np.log(abs(np.linalg.det(np.vstack([w, P[:-1]])))))
        self._forward()

    def _outer(self, a, L):
        """Eigen-quantities and per-omega log evidence (incl. prior) of outer blocks [a, a + L)."""
        A, h, q = (S[a + L] - S[a] for S in (self.PA, self.PH, self.PQ))
        lam, vec = np.linalg.eigh(self.Linv @ A @ self.Linv.T)
        lam = np.maximum(lam, 0.0)
        U = np.einsum("ij,sjk->sik", self.Linv.T, vec)
        g = np.einsum("sik,si->sk", U, h)
        w2 = self.omega2[:, None]
        lm = (-0.5 * q[:, None] + 0.5 * np.sum(g[:, None] ** 2 / (lam[:, None] + 1 / w2), axis=2)
              - 0.5 * np.sum(np.log1p(lam[:, None] * w2), axis=2) + self.log_omega_prior)
        return lam, U, g, lm

    def _forward(self):
        """f(a, d): all inner tilings of [a, a + d); F(c): all nested segmentations of [0, c)."""
        n, Lm = self.n, self.Lm
        f = np.full((n + 1, Lm + 1), -np.inf)
        f[:, 0] = 0.0
        a = np.arange(n + 1)[:, None]
        for d in range(1, Lm + 1):
            dp = np.arange(d)
            f[:, d] = logsumexp(f[:, :d] + self.lq[d - dp] + self.mb[d - dp, a + dp], axis=1)
        F = np.full(n + 1, -np.inf)
        F[0] = 0.0
        for c in range(1, n + 1):
            L = np.arange(1, min(Lm, c) + 1)
            F[c] = logsumexp(F[c - L] + self.lo[L] + self.mdp[L, c - L] + f[c - L, L])
        self.f, self.F = f, F
        self.log_evidence = float(F[n] + self.const)

    def sample(self, n_draws=1000, seed=0):
        """Exact backward sampling, vectorized over draws.

        Returns depth (draws x n x G), P(outer boundary) and P(inner-only boundary) between
        bases j and j + 1, and the posterior mean ICC per base.
        """
        rng = np.random.default_rng(seed)
        n, G, f = self.n, self.G, self.f
        b, delta = np.zeros((n_draws, n)), np.zeros((n_draws, n, G))
        icc = np.zeros((n_draws, n))
        inner, outer = np.zeros((n_draws, n), bool), np.zeros((n_draws, n), bool)
        Ls = np.arange(1, self.Lm + 1)
        c = np.full(n_draws, n)
        while np.any(c > 0):
            t = np.flatnonzero(c > 0)
            a_ = np.maximum(c[t, None] - Ls, 0)
            W = np.where(Ls <= c[t, None], self.F[a_] + self.lo[Ls] + self.mdp[Ls, a_] + f[a_, Ls], -np.inf)
            a = c[t] - Ls[_pick(W, rng)]
            outer[t, a] = True
            delta[t, a], icc[t, a] = self._draw_contrasts(a, c[t] - a, rng)
            e = c[t].copy()
            while np.any(e > a):
                k = e > a
                tk, ak, ek = t[k], a[k], e[k]
                s = np.maximum(ek[:, None] - Ls, ak[:, None])
                W = np.where(ek[:, None] - Ls >= ak[:, None], f[ak[:, None], s - ak[:, None]] + self.lq[Ls] + self.mb[Ls, s], -np.inf)
                s = s[np.arange(tk.size), _pick(W, rng)]
                prec = self.Pa[ek] - self.Pa[s] + 1 / self.b_var
                b[tk, s] = (self.Ph[ek] - self.Ph[s]) / prec + rng.standard_normal(s.size) / np.sqrt(prec)
                inner[tk, s] = True
                e[k] = s
            c[t] = a
        fill_in = np.maximum.accumulate(np.where(inner, np.arange(n), 0), axis=1)
        fill_out = np.maximum.accumulate(np.where(outer, np.arange(n), 0), axis=1)
        depth = np.take_along_axis(b, fill_in, 1)[..., None] + np.take_along_axis(delta, fill_out[..., None], 1)
        return (depth.astype(np.float32), outer[:, 1:].mean(0), (inner & ~outer)[:, 1:].mean(0),
                np.take_along_axis(icc, fill_out, 1).mean(0))

    def _draw_contrasts(self, a, L, rng):
        blocks, inv = np.unique(np.stack([a, L], 1), axis=0, return_inverse=True)
        lam, U, g, lm = (z[inv.ravel()] for z in self._outer(blocks[:, 0], blocks[:, 1]))
        r = np.exp(lm - logsumexp(lm, axis=1, keepdims=True))
        dd = 1 / (lam + 1 / self.omega2[_pick(lm, rng), None])
        z = g * dd + rng.standard_normal(dd.shape) * np.sqrt(dd)
        return np.einsum("gi,sik,sk->sg", self.Q, U, z), r @ self.icc_grid


def _pick(log_w, rng):
    """One categorical draw per row (Gumbel-max)."""
    return np.argmax(log_w + rng.gumbel(size=log_w.shape), axis=1)
