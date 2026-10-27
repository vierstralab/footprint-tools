"""Nested product-partition model for group footprint summaries.

Outer blocks hold the group contrasts constant; inner blocks, nested inside them,
hold the consensus level constant. Block evidences are closed-form, the sum over
all nested segmentations is exact, and posterior draws are exact.
"""

from dataclasses import dataclass

import numpy as np
from scipy.special import logsumexp

from .config import DEFAULT_NESTED, NestedConfig
from .priors import gamma_length_log_mass, half_cauchy_log_mass, outer_length_log_weight


def _cum(a: np.ndarray) -> np.ndarray:
    return np.concatenate([np.zeros((1,) + a.shape[1:]), np.cumsum(a, axis=0)])


@dataclass(frozen=True)
class Draws:
    depth: np.ndarray          # draws x n x G
    p_outer: np.ndarray        # n - 1, boundary before base j + 1
    p_inner: np.ndarray        # n - 1, inner-only boundary
    icc: np.ndarray            # n


class NestedModel:
    def __init__(
        self,
        x: np.ndarray,
        v: np.ndarray,
        config: NestedConfig = DEFAULT_NESTED,
        omega2: np.ndarray | None = None,
        log_omega_prior: np.ndarray | None = None,
    ):
        x, v = np.asarray(x, float), np.asarray(v, float)
        if x.shape != v.shape or x.ndim != 2:
            raise ValueError("x and v must both be groups x positions")
        if not np.all(np.isfinite(x)) or np.any(~(v > 0)):
            raise ValueError("x must be finite and v positive")
        self.config = config
        self.G, self.n = G, n = x.shape
        self.l_max = Lm = min(config.l_max, n)
        self.omega2 = config.omega2() if omega2 is None else np.asarray(omega2, float)
        self.log_omega_prior = (
            half_cauchy_log_mass(self.omega2, config.omega_scale)
            if log_omega_prior is None else np.asarray(log_omega_prior, float)
        )
        self.lq = gamma_length_log_mass(Lm, config.inner_mean, config.len_shape, config.len_min)
        lq_out = gamma_length_log_mass(Lm, config.outer_mean, config.len_shape, config.len_min)
        self.lo = outer_length_log_weight(self.lq, lq_out)

        ve = config.kappa * v + config.tau2_contrast
        w = 1 / ve.mean(axis=1)
        self.w = w = w / w.sum()
        P = np.eye(G) - np.outer(np.ones(G), w)
        self._consensus(w @ x, config.tau2 + (w ** 2) @ ve, config.b_var)
        if G > 1:
            self._contrasts(P, (P @ x)[:-1], ve)
        else:
            self.mdp = np.zeros((Lm + 1, n + 1))
            self.const_c = 0.0
        self.icc_grid = self.omega2 / (self.omega2 + config.sigma2_bar)
        T = np.vstack([w, P[:-1]])
        self.const = self.const_m + self.const_c + n * np.log(abs(np.linalg.det(T)))

    def _consensus(self, m, V, b_var):
        n, Lm = self.n, self.l_max
        self.Pa, self.Ph = _cum(1 / V), _cum(m / V)
        Pq = _cum(m ** 2 / V)
        self.b_var = b_var
        self.mb = np.full((Lm + 1, n + Lm + 1), -np.inf)
        for L in range(1, Lm + 1):
            a = np.arange(n - L + 1)
            A, h, q = self.Pa[a + L] - self.Pa[a], self.Ph[a + L] - self.Ph[a], Pq[a + L] - Pq[a]
            prec = A + 1 / b_var
            self.mb[L, : n - L + 1] = -0.5 * q + 0.5 * h ** 2 / prec - 0.5 * np.log(b_var * prec)
        self.const_m = -0.5 * np.log(2 * np.pi * V).sum()

    def _contrasts(self, P, c, ve):
        n, Lm, G = self.n, self.l_max, self.G
        C = np.einsum("gh,hj,lh->jgl", P, ve, P)[:, :-1, :-1]
        Ci = np.linalg.inv(C)
        self.PA = _cum(Ci)
        self.PH = _cum(np.einsum("jgh,hj->jg", Ci, c))
        self.PQ = _cum(np.einsum("gj,jgh,hj->j", c, Ci, c))
        K = (P @ P.T)[:-1, :-1]
        self.Linv = np.linalg.inv(np.linalg.cholesky(np.linalg.inv(K)))
        self.Q = np.vstack([np.eye(G - 1), -self.w[:-1] / self.w[-1]])
        self.mdp = np.full((Lm + 1, n + 1), -np.inf)
        for L in range(1, Lm + 1):
            a = np.arange(n - L + 1)
            self.mdp[L, : n - L + 1] = logsumexp(self._outer_terms(a, np.full(a.size, L))[3], axis=1)
        self.const_c = -0.5 * np.linalg.slogdet(2 * np.pi * C)[1].sum()

    def _outer_terms(self, a, L):
        """Eigen-quantities and per-omega log evidence (incl. prior) of outer blocks [a, a+L)."""
        A = self.PA[a + L] - self.PA[a]
        h = self.PH[a + L] - self.PH[a]
        q = self.PQ[a + L] - self.PQ[a]
        lam, vec = np.linalg.eigh(self.Linv @ A @ self.Linv.T)
        lam = np.maximum(lam, 0.0)
        U = np.einsum("ij,sjk->sik", self.Linv.T, vec)
        g = np.einsum("sik,si->sk", U, h)
        w2 = self.omega2[None, :, None]
        lm = (-0.5 * q[:, None]
              + 0.5 * np.sum(g[:, None] ** 2 / (lam[:, None] + 1 / w2), axis=2)
              - 0.5 * np.sum(np.log1p(lam[:, None] * w2), axis=2)
              + self.log_omega_prior[None])
        return lam, U, g, lm

    def partition(self):
        """Inner tiling weights f(a, d) and outer forward sums F(c), in log space (cached)."""
        if getattr(self, "_part", None) is None:
            self._part = self._forward()
        return self._part

    def _forward(self):
        n, Lm = self.n, self.l_max
        f = np.full((n + 1, Lm + 1), -np.inf)
        f[:, 0] = 0.0
        idx = np.arange(n + 1)[:, None]
        for d in range(1, Lm + 1):
            dp = np.arange(d)[None]
            f[:, d] = logsumexp(f[:, :d] + self.lq[d - dp] + self.mb[d - dp, idx + dp], axis=1)
        F = np.full(n + 1, -np.inf)
        F[0] = 0.0
        for c in range(1, n + 1):
            L = np.arange(1, min(Lm, c) + 1)
            F[c] = logsumexp(F[c - L] + self.lo[L] + self.mdp[L, c - L] + f[c - L, L])
        return f, F

    def log_evidence(self) -> float:
        return float(self.partition()[1][-1] + self.const)

    def sample(self, n_draws: int = 1000, rng=None) -> Draws:
        rng = np.random.default_rng(rng)
        n, Lm, G = self.n, self.l_max, self.G
        f, F = self.partition()
        T = np.arange(n_draws)
        inner_val = np.zeros((n_draws, n))
        inner_mark = np.zeros((n_draws, n), bool)
        outer_val = np.zeros((n_draws, n, G), np.float32)
        outer_icc = np.full((n_draws, n), np.nan)
        outer_mark = np.zeros((n_draws, n), bool)
        Ls = np.arange(1, Lm + 1)

        c = np.full(n_draws, n)
        while np.any(c > 0):
            t = T[c > 0]
            ct = c[t, None]
            a_ = np.maximum(ct - Ls, 0)
            W = F[a_] + self.lo[Ls] + self.mdp[Ls, a_] + f[a_, Ls]
            W[Ls > ct] = -np.inf
            L = Ls[np.argmax(W + rng.gumbel(size=W.shape), axis=1)]
            a = c[t] - L
            outer_mark[t, a] = True
            if G > 1:
                outer_val[t, a], outer_icc[t, a] = self._draw_contrasts(a, L, rng)
            self._draw_inner(t, a, c[t].copy(), f, inner_val, inner_mark, rng)
            c[t] = a

        fill_in = np.maximum.accumulate(np.where(inner_mark, np.arange(n), 0), axis=1)
        fill_out = np.maximum.accumulate(np.where(outer_mark, np.arange(n), 0), axis=1)
        b = np.take_along_axis(inner_val, fill_in, axis=1)
        delta = np.take_along_axis(outer_val, fill_out[..., None], axis=1)
        depth = (b[..., None] + delta).astype(np.float32)
        icc = np.take_along_axis(outer_icc, fill_out, axis=1).mean(axis=0)
        return Draws(
            depth=depth,
            p_outer=outer_mark[:, 1:].mean(axis=0),
            p_inner=(inner_mark & ~outer_mark)[:, 1:].mean(axis=0),
            icc=icc,
        )

    def _draw_contrasts(self, a, L, rng):
        blocks, inverse = np.unique(np.stack([a, L], 1), axis=0, return_inverse=True)
        inverse = inverse.ravel()
        lam, U, g, lm = self._outer_terms(blocks[:, 0], blocks[:, 1])
        r = np.exp(lm - logsumexp(lm, axis=1, keepdims=True))
        icc = (r @ self.icc_grid)[inverse]
        cdf = np.cumsum(r[inverse], axis=1)
        k = np.minimum((cdf < rng.random((a.size, 1)) * cdf[:, -1:]).sum(axis=1), self.omega2.size - 1)
        dd = 1 / (lam[inverse] + 1 / self.omega2[k, None])
        z = g[inverse] * dd + rng.standard_normal(dd.shape) * np.sqrt(dd)
        delta = np.einsum("gi,sik,sk->sg", self.Q, U[inverse], z)
        return delta, icc

    def _draw_inner(self, t, a, e, f, inner_val, inner_mark, rng):
        Lm = self.l_max
        ds = np.arange(1, Lm + 1)
        while np.any(e > a):
            act = e > a
            tt, aa, ee = t[act], a[act], e[act]
            s = ee[:, None] - ds
            ok = s >= aa[:, None]
            sc = np.where(ok, s, aa[:, None])
            W = f[aa[:, None], sc - aa[:, None]] + self.lq[ds] + self.mb[ds, sc]
            W[~ok] = -np.inf
            s = sc[np.arange(tt.size), np.argmax(W + rng.gumbel(size=W.shape), axis=1)]
            prec = self.Pa[ee] - self.Pa[s] + 1 / self.b_var
            mean = (self.Ph[ee] - self.Ph[s]) / prec
            inner_val[tt, s] = mean + rng.standard_normal(s.size) / np.sqrt(prec)
            inner_mark[tt, s] = True
            e[act] = s
