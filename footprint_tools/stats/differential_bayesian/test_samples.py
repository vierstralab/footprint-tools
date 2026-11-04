from types import SimpleNamespace

import numpy as np
from scipy.stats import nbinom

from ..hyper import fit_nugget, nugget_stats
from ..samples import BINS, fit_maps, profile_stats, rho
from ..summaries import nb_loglik


def test_fit_maps_recovers_contrast_and_floor():
    rng = np.random.default_rng(0)
    k, f = rng.uniform(0.7, 1.3, 12), rng.uniform(0.0, 0.1, 12)
    k /= np.median(k)
    d = BINS[:-1] + 0.125
    E = np.full((12, d.size), 2e4)
    O = rng.poisson(E * rho(d, k[:, None], f[:, None]))
    k_hat, f_hat = fit_maps(O, E)
    assert np.abs(k_hat - k).max() < 0.05 and np.abs(f_hat - f).max() < 0.01


def test_profile_stats_counts_only_binned_bases():
    rng = np.random.default_rng(1)
    exp = rng.poisson(8, (6, 200)).astype(float)
    obs = rng.poisson(exp)
    O, E = profile_stats(obs, exp, [f"s{i}" for i in range(6)], ["a", "b"] * 3)
    assert O.shape == E.shape == (6, BINS.size - 1)
    assert O.sum() <= obs.sum() and E.sum() <= exp[exp >= 4].sum()


def test_nb_loglik_matches_scipy_above_breakpoints():
    mu = np.r_[1.0, 5.0, 20.0, 0.1, 0.2, 0.5, 0.9, 0.95, 0.98]
    r = np.r_[1.0, 3.0, 10.0, 30.0, 50.0, 0.3, 0.2, 0.1, 0.05, 0.02, 0.0, 0.0, 0.0, 0.0, 0.0]
    model = SimpleNamespace(mu_params=mu, r_params=r)
    exp, obs = np.array([[4.0, 12.0, 60.0]]), np.array([[3.0, 9.0, 70.0]])
    share = np.array([[1.0, 0.5]])
    ll = nb_loglik([model], obs, exp, share)
    for t in range(2):
        e = exp[0] * share[0, t]
        seg_m, seg_r = (e[:, None] >= mu[:2]).sum(1), (e[:, None] >= r[:4]).sum(1)
        m, size = mu[3 + seg_m] + mu[6 + seg_m] * e, 1 / (r[5 + seg_r] + r[10 + seg_r] * e)
        assert np.allclose(ll[t, :, 0], nbinom.logpmf(obs[0], size, size / (size + m)))


def test_nugget_recovers_planted_values():
    rng = np.random.default_rng(2)
    t = np.array([0.0, 0.3, 0.1])
    v = rng.uniform(0.05, 0.2, (3, 20000))
    x = rng.normal(0, 1, 20000) + rng.normal(0, np.sqrt(v + t[:, None]))
    num, A = nugget_stats(x, v)
    h = {"group_names": np.array(["a", "b", "c"]), "nugget.num": num, "nugget.A": A}
    assert np.allclose(list(fit_nugget([h]).values()), t, atol=0.03)
