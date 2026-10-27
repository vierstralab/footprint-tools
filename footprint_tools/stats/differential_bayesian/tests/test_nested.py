import itertools
import json

import numpy as np
import pytest
from scipy.special import logsumexp
from scipy.stats import multivariate_normal

from ..calls import runs
from ..config import NestedConfig
from ..nested import NestedModel
from ..results import NestedFit, fit_nested

SMALL = NestedConfig(l_max=10, len_min=1, inner_mean=2.0, outer_mean=3.0, n_draws=200)


def _data(G, n, seed=0, flat_profile=False):
    rng = np.random.default_rng(seed)
    x = rng.normal(-0.5, 0.6, (G, n))
    x[0, n // 2:] -= 1.0
    if flat_profile:
        v = np.outer(rng.uniform(0.01, 0.05, G), rng.uniform(0.5, 2.0, n))
    else:
        v = rng.uniform(0.01, 0.05, (G, n))
    return x, v


def _compositions(n):
    for k in range(n):
        for cuts in itertools.combinations(range(1, n), k):
            yield list(zip((0,) + cuts, cuts + (n,)))


def _brute(M):
    """Enumerate every nested segmentation: log weights and boundary indicators."""
    logw, outer_b, inner_b = [], [], []
    for outer in _compositions(M.n):
        for inners in itertools.product(*(list(_compositions(e - a)) for a, e in outer)):
            w = 0.0
            ob, ib = np.zeros(M.n - 1), np.zeros(M.n - 1)
            for (a, e), inner in zip(outer, inners):
                w += M.lo[e - a] + M.mdp[e - a, a]
                if a > 0:
                    ob[a - 1] = 1
                for s, t in inner:
                    w += M.lq[t - s] + M.mb[t - s, a + s]
                    if s > 0:
                        ib[a + s - 1] = 1
            logw.append(w)
            outer_b.append(ob)
            inner_b.append(ib)
    return np.array(logw), np.array(outer_b), np.array(inner_b)


def test_partition_matches_enumeration():
    x, v = _data(3, 6)
    M = NestedModel(x, v, SMALL)
    logw, _, _ = _brute(M)
    assert np.isclose(M.partition()[1][-1], logsumexp(logw), atol=1e-8)


def test_single_block_evidence_matches_dense_gaussian():
    G, n, omega2 = 4, 5, 0.3
    x, v = _data(G, n, flat_profile=True)
    cfg = NestedConfig(l_max=n, len_min=1)
    M = NestedModel(x, v, cfg, omega2=np.array([omega2]), log_omega_prior=np.zeros(1))
    model_logz = M.mb[n, 0] + M.mdp[n, 0] + M.const

    P = np.eye(G) - np.outer(np.ones(G), M.w)
    J = np.ones((n, n))
    cov = (cfg.b_var * np.ones((G * n, G * n))
           + omega2 * np.kron(P @ P.T, J)
           + cfg.tau2 * np.kron(np.ones((G, G)), np.eye(n))
           + np.diag((cfg.kappa * v).ravel()))
    dense = multivariate_normal(np.zeros(G * n), cov).logpdf(x.ravel())
    assert np.isclose(model_logz, dense, atol=1e-8)


def test_sampler_boundary_marginals():
    x, v = _data(3, 6, seed=1)
    M = NestedModel(x, v, SMALL)
    logw, ob, ib = _brute(M)
    p = np.exp(logw - logsumexp(logw))
    d = M.sample(20000, 0)
    assert np.allclose(d.p_outer, p @ ob, atol=0.02)
    assert np.allclose(d.p_inner, p @ np.clip(ib - ob, 0, 1), atol=0.02)


def test_single_group():
    x, v = _data(1, 30)
    fit = fit_nested(x, v, ["only"], SMALL)
    assert fit.calls == []
    assert np.all(np.isnan(fit.icc))
    assert np.all(fit.delta_mean == 0)


def test_no_group_differences_without_contrast_signal():
    rng = np.random.default_rng(3)
    shared = np.repeat([-1.5, 0.2, -1.0], 20)
    x = shared + rng.normal(0, 0.05, (5, 60))
    fit = fit_nested(x, np.full_like(x, 0.01), list("abcde"), NestedConfig(n_draws=300))
    assert fit.calls == []
    assert fit.footprints and fit.footprints[0][0] == 0


def test_roundtrip():
    x, v = _data(3, 40)
    fit = fit_nested(x, v, ["a", "b", "c"], SMALL.with_(store_draws=True))
    back = NestedFit.from_dict({k: np.asarray(val) for k, val in fit.to_dict().items()})
    assert back.group_names == fit.group_names
    assert back.config == fit.config
    assert np.array_equal(back.p_more, fit.p_more)
    assert back.draws.shape == (SMALL.n_draws, 40, 3)
    assert json.loads(back.config_json)["l_max"] == SMALL.l_max


def test_runs():
    assert runs(np.array([0, 1, 1, 0, 1], bool)) == [(1, 3), (4, 5)]
    assert runs(np.zeros(3, bool)) == []


def test_removed_loaders_point_to_replacement():
    pytest.importorskip("genome_tools")
    from .. import api
    with pytest.raises(ImportError, match="NestedSegmentationLoader"):
        api.GroupMeanSegmentationLoader


def test_save_load_and_slice(tmp_path):
    gt = pytest.importorskip("genome_tools")
    from genome_tools.plotting.modular_plot.api import DataBundle
    from .. import api
    from ..base_differential import Differential

    G, n = 3, 40
    x, v = _data(G, n)
    mu_x = np.linspace(-3, 3, 61)
    loglik = -0.5 * (mu_x[None, None] - x[..., None]) ** 2 / v[..., None]
    diff = Differential.from_dict(dict(
        group_names=np.array(["a", "b", "c"]), mu_x=mu_x, sig2_x=np.geomspace(1e-3, 1, 5),
        loglik_mu=loglik, log_sig2_prior=np.zeros((G, 5)), log_mu_prior=np.zeros(61)))
    assert np.allclose(diff.mu_mean, x, atol=0.02)

    data = DataBundle(interval=gt.GenomicInterval("chr1", 1000, 1000 + n))
    data.differential = diff
    data = api.NestedSegmentationLoader()._load(data, config=SMALL)
    path = tmp_path / "dhs.npz"
    api.save_data_results(data, path)

    back = api.load_data_results(path)
    assert np.array_equal(back.nested.p_any, data.nested.p_any)
    assert back.nested.calls == data.nested.calls

    sub = DataBundle(interval=gt.GenomicInterval("chr1", 1010, 1030))
    sub = api.load_data_results(path, sub)
    assert sub.nested.p_any.shape == (20,)
    assert sub.nested.p_outer.shape == (19,)
    assert np.allclose(sub.nested.p_more, data.nested.p_more[10:30])
