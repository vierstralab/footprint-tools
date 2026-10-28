import itertools

import numpy as np
from scipy.special import logsumexp
from scipy.stats import multivariate_normal

from ..config import NestedConfig
from ..nested import NestedModel
from ..results import fit_nested, runs

SMALL = NestedConfig(l_max=10, len_min=1, inner_mean=2.0, outer_mean=3.0, n_draws=200)


def _data(G, n, seed=0, flat_profile=False):
    rng = np.random.default_rng(seed)
    x = rng.normal(-0.5, 0.6, (G, n))
    x[0, n // 2:] -= 1.0
    profile = rng.uniform(0.5, 2.0, n) if flat_profile else rng.uniform(0.5, 2.0, (G, n))
    return x, rng.uniform(0.01, 0.05, (G, 1)) * profile


def _compositions(n):
    for k in range(n):
        for cuts in itertools.combinations(range(1, n), k):
            yield list(zip((0,) + cuts, cuts + (n,)))


def _enumerate(M):
    """Log weight, outer and inner-only boundary indicators of every nested segmentation."""
    rows = []
    for outer in _compositions(M.n):
        for inners in itertools.product(*(list(_compositions(e - a)) for a, e in outer)):
            w, ob, ib = 0.0, np.zeros(M.n - 1), np.zeros(M.n - 1)
            for (a, e), inner in zip(outer, inners):
                w += M.lo[e - a] + M.mdp[e - a, a]
                ob[a - 1] = a > 0
                for s, t in inner:
                    w += M.lq[t - s] + M.mb[t - s, a + s]
                    ib[a + s - 1] = s > 0
            rows.append((w, ob, ib))
    w, ob, ib = map(np.array, zip(*rows))
    return w, ob, ib


def test_partition_matches_enumeration():
    M = NestedModel(*_data(3, 6), SMALL)
    assert np.isclose(M.F[-1], logsumexp(_enumerate(M)[0]), atol=1e-8)


def test_block_evidence_matches_dense_gaussian():
    G, n, omega2 = 4, 5, 0.3
    x, v = _data(G, n, flat_profile=True)
    cfg = NestedConfig(l_max=n, len_min=1)
    M = NestedModel(x, v, cfg, omega2=np.array([omega2]), log_omega_prior=np.zeros(1))
    P = np.eye(G) - np.outer(np.ones(G), M.w)
    cov = (cfg.b_var + omega2 * np.kron(P @ P.T, np.ones((n, n))) + cfg.tau2 * np.kron(np.ones((G, G)), np.eye(n))
           + np.diag(cfg.kappa * v.ravel()))
    dense = multivariate_normal(np.zeros(G * n), cov).logpdf(x.ravel())
    assert np.isclose(M.mb[n, 0] + M.mdp[n, 0] + M.const, dense, atol=1e-8)


def test_sampler_boundary_marginals():
    M = NestedModel(*_data(3, 6, seed=1), SMALL)
    w, ob, ib = _enumerate(M)
    p = np.exp(w - logsumexp(w))
    _, p_outer, p_inner, _ = M.sample(20000, 0)
    assert np.allclose(p_outer, p @ ob, atol=0.02)
    assert np.allclose(p_inner, p @ (ib - ob).clip(0), atol=0.02)


def test_shared_signal_gives_footprints_but_no_calls():
    rng = np.random.default_rng(3)
    x = np.repeat([-1.5, 0.2, -1.0], 20) + rng.normal(0, 0.05, (5, 60))
    fit = fit_nested(x, np.full_like(x, 0.01), list("abcde"), NestedConfig(n_draws=300))
    assert fit.calls == [] and fit.footprints[0][0] == 0


def test_runs():
    assert runs(np.array([0, 1, 1, 0, 1], bool)) == [(1, 3), (4, 5)]
