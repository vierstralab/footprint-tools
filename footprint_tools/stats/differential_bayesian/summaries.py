"""Per-group and per-sample depletion summaries from the negative-binomial cleavage model.

theta_sj ~ N(mu_gj, sig2), sig2 ~ inv-chi2 fitted per group; counts ~ NB(2^theta * expected).
Group: posterior mean and variance of mu (flat prior), the input of the nested model.
Sample: posterior mean and SD of theta under the prior N(group mean, sig2 scale + group variance).
"""

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.special import logsumexp, ndtr
from scipy.stats import invgamma

from footprint_tools.stats.differential import compute_logpmf_values

from .config import DEFAULT_SUMMARY, SummaryConfig


@dataclass(frozen=True)
class Summaries:
    group_names: tuple[str, ...]
    x: np.ndarray               # G x n, E[mu | counts]
    v: np.ndarray               # G x n, Var[mu | counts]
    sample_names: tuple[str, ...]
    sample_groups: tuple[str, ...]
    theta_mean: np.ndarray      # S x n
    theta_sd: np.ndarray        # S x n


def fit_summaries(groups: pd.Series, obs: pd.DataFrame, exp: pd.DataFrame, disp_models: pd.Series,
                  config: SummaryConfig = DEFAULT_SUMMARY) -> Summaries:
    """groups: sample -> group label; obs, exp: samples x positions; disp_models: sample -> NB model."""
    idx, labels = groups.index, groups.to_numpy().astype(str)
    names = tuple(pd.unique(labels))
    obs, exp = (np.ascontiguousarray(a.loc[idx].to_numpy(float)) for a in (obs, exp))
    mu, sig2, theta = config.grids()
    loglik = compute_logpmf_values(disp_models.loc[idx].to_numpy(), obs, exp, theta[0], theta[-1], theta.size)
    loglik[:, exp.T <= 0] = 0.0   # no expected cleavage: the base carries no information on theta
    kernel = _kernel(mu, sig2, theta)
    log2_ratio = np.log2((obs + 1) / (exp + 1))

    x, v = np.empty((len(names), obs.shape[1])), np.empty((len(names), obs.shape[1]))
    theta_mean, theta_sd = np.empty(obs.shape), np.empty(obs.shape)
    for g, name in enumerate(names):
        s = np.flatnonzero(labels == name)
        if s.size < 2:
            raise ValueError(f"group {name} needs at least two samples")
        sig2_scale, log_prior = _sig2_prior(np.var(log2_ratio[s], axis=0, ddof=1), sig2)
        ll_mu = _group_loglik(kernel, loglik[:, :, s], log_prior, config.chunk)
        x[g], v[g] = _moments(mu, ll_mu)
        prior = -0.5 * (theta - x[g][:, None]) ** 2 / (sig2_scale + v[g][:, None])
        m, var = _moments(theta, loglik[:, :, s].transpose(2, 1, 0) + prior)
        theta_mean[s], theta_sd[s] = m, np.sqrt(var)
    return Summaries(names, x, v, tuple(map(str, idx)), tuple(labels), theta_mean, theta_sd)


def _kernel(mu, sig2, theta):
    """P(theta bin | mu, sig2), shape (sig2 * mu) x theta."""
    edges = np.r_[-np.inf, (theta[1:] + theta[:-1]) / 2, np.inf]
    cdf = ndtr((edges - mu[None, :, None]) / np.sqrt(sig2)[:, None, None])
    mass = np.diff(cdf, axis=-1)
    return (mass / mass.sum(-1, keepdims=True)).reshape(-1, theta.size)


def _sig2_prior(variances, sig2):
    """Inverse-gamma (scaled inv-chi2) fitted to per-base sample variances; log mass on the sig2 grid."""
    variances = variances[np.isfinite(variances) & (variances > 0)]
    a, _, scale = invgamma.fit(variances, floc=0)
    log_mass = invgamma.logpdf(sig2, a, scale=scale) + np.log(sig2)   # geometric grid: d sig2 ~ sig2
    return scale / a, log_mass - logsumexp(log_mass)


def _group_loglik(kernel, loglik, log_prior, chunk):
    """log p(group counts | mu) per position: theta and sig2 integrated out. loglik is theta x n x samples."""
    n_sig2, n_mu = log_prior.size, kernel.shape[0] // log_prior.size
    out = np.empty((loglik.shape[1], n_mu))
    for lo in range(0, loglik.shape[1], chunk):
        block = loglik[:, lo:lo + chunk]
        shift = block.max(axis=0)
        with np.errstate(divide="ignore"):
            raw = np.log(kernel @ np.exp(block - shift).reshape(block.shape[0], -1))
        raw = (raw.reshape(n_sig2, n_mu, *block.shape[1:]) + shift).sum(-1)
        out[lo:lo + chunk] = logsumexp(raw + log_prior[:, None, None], axis=0).T
    return out


def _moments(grid, log_weight):
    p = np.exp(log_weight - logsumexp(log_weight, axis=-1, keepdims=True))
    mean = p @ grid
    var = np.maximum(p @ grid ** 2 - mean ** 2, np.diff(grid).mean() ** 2 / 12)
    return mean, var
