import itertools

import numpy as np
from scipy.special import logsumexp
from scipy.stats import multivariate_normal

from ..config import ModelConfig
from ..model import FootprintModel
from ..results import fit_footprints, runs

SMALL = ModelConfig(l_max=10, cons_min=1, cons_mean=2.0, pattern_min=1, pattern_mean=3.0, n_draws=200)


def _data(G, n, seed=0, flat_profile=False):
    rng = np.random.default_rng(seed)
    x = rng.normal(-0.5, 0.6, (G, n))
    x[0, n // 2:] -= 1.0
    profile = rng.uniform(0.5, 2.0, n) if flat_profile else rng.uniform(0.5, 2.0, (G, n))
    return x, rng.uniform(0.01, 0.05, (G, 1)) * profile


def _enumerate(E, lq, n):
    """Log weight and boundary indicators of every partition of [0, n)."""
    rows = []
    for k in range(n):
        for cuts in itertools.combinations(range(1, n), k):
            blocks = list(zip((0,) + cuts, cuts + (n,)))
            b = np.zeros(n - 1); b[[c - 1 for c in cuts]] = 1
            rows.append((sum(lq[e - s] + E[e - s, s] for s, e in blocks), b))
    w, b = map(np.array, zip(*rows))
    return w, b


def test_forward_and_sampler_match_enumeration():
    x, v = _data(3, 6, seed=1)
    M = FootprintModel(x, v, np.ones_like(x), SMALL)
    (E_b, lq_b), (E_d, lq_d), _ = M.layers()
    _, p_pattern, p_consensus, _ = M.sample(20000, 0)
    for (E, lq), z, p_hat in (((E_b, lq_b), M.log_evidence()[0], p_consensus), ((E_d, lq_d), M.log_evidence()[1], p_pattern)):
        w, b = _enumerate(E, lq, M.n)
        assert np.isclose(z, logsumexp(w), atol=1e-8)
        assert np.allclose(p_hat, np.exp(w - logsumexp(w)) @ b, atol=0.02)


def test_block_evidence_matches_dense_gaussian():
    G, n, omega2 = 4, 5, 0.3
    x, v = _data(G, n, flat_profile=True)
    cfg = ModelConfig(l_max=n)
    M = FootprintModel(x, v, np.ones_like(x), cfg, omega2=np.array([omega2]))
    P = np.eye(G) - np.outer(np.ones(G), M.w)
    cov = (cfg.b_var + omega2 * np.kron(P @ P.T, np.ones((n, n))) + cfg.tau2 * np.kron(np.ones((G, G)), np.eye(n))
           + np.diag(v.ravel()))
    dense = multivariate_normal(np.zeros(G * n), cov).logpdf(x.ravel())
    assert np.isclose(M.E_b[n, 0] + M.E_d_omega[n, 0, 0] + M.const, dense, atol=1e-8)


def test_shared_signal_gives_footprints_but_no_calls():
    rng = np.random.default_rng(3)
    x = np.repeat([-1.5, 0.2, -1.0], 20) + rng.normal(0, 0.05, (5, 60))
    fit = fit_footprints(x, np.full_like(x, 0.01), np.full_like(x, 0.4), list("abcde"), ModelConfig(n_draws=300))
    assert fit.calls == [] and fit.footprints[0][0] == 0
    assert np.all(fit.icc < 0.2)


def test_runs():
    assert runs(np.array([0, 1, 1, 0, 1], bool)) == [(1, 3), (4, 5)]
