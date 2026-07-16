from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
import pandas as pd
import scipy.optimize
from scipy.special import logsumexp, ndtr

from footprint_tools.stats import differential as differential_core
from footprint_tools.stats.distributions import invchi2

from .config import (
    DEFAULT_DIFFERENTIAL,
    DEFAULT_MEAN_SEGMENTATION,
    DifferentialConfig,
    MeanSegmentationConfig,
)
from .posterior import GridPosterior, normal_grid_log_mass, normalize_log_mass
from .segmentation import LengthPrior, Segmentation, segment


@dataclass(frozen=True, slots=True)
class Differential:
    group_names: tuple[str, ...]
    mu_x: np.ndarray
    sig2_x: np.ndarray
    loglik_mu: np.ndarray
    loglik_mu_sig2: np.ndarray
    log_sig2_prior: np.ndarray
    log_mu_prior: np.ndarray
    theta_x: np.ndarray | None = None

    @property
    def group_mu(self) -> np.ndarray:
        return self.mu_x[np.argmax(self.loglik_mu, axis=2)]

    @property
    def common_mu(self) -> np.ndarray:
        return self.mu_x[np.argmax(self.loglik_mu.sum(axis=0), axis=1)]

    def posterior(self, log_mu_prior: np.ndarray | None = None) -> GridPosterior:
        """Pointwise posterior over mu, shape group x position x mu."""
        prior = self.log_mu_prior if log_mu_prior is None else normalize_log_mass(
            log_mu_prior, self.mu_x.size
        )
        return GridPosterior(self.mu_x, self.loglik_mu + prior[None, None, :])

    def to_npz(self, path) -> None:
        np.savez_compressed(
            path,
            group_names=self.group_names,
            mu_x=self.mu_x,
            sig2_x=self.sig2_x,
            loglik_mu=self.loglik_mu,
            loglik_mu_sig2=self.loglik_mu_sig2,
            log_sig2_prior=self.log_sig2_prior,
            log_mu_prior=self.log_mu_prior,
            **({"theta_x": self.theta_x} if self.theta_x is not None else {}),
        )

    @classmethod
    def from_npz(cls, path) -> "Differential":
        with np.load(path, allow_pickle=False) as x:
            return cls(
                tuple(x["group_names"].tolist()),
                x["mu_x"],
                x["sig2_x"],
                x["loglik_mu"],
                x["loglik_mu_sig2"],
                x["log_sig2_prior"],
                x["log_mu_prior"],
                x["theta_x"] if "theta_x" in x.files else None,
            )


class DifferentialModel:
    def __init__(self, config: DifferentialConfig = DEFAULT_DIFFERENTIAL):
        self.config = config

    def fit(
        self,
        groups_data: pd.Series,
        obs: pd.DataFrame,
        exp: pd.DataFrame,
        disp_models: pd.Series,
        selected_groups: Sequence[str] | None = None,
    ) -> Differential:
        names, indices, obs, exp, models = self._align(
            groups_data, obs, exp, disp_models, selected_groups
        )
        mu_x = self.config.mu_x()
        sig2_x = self.config.sig2_x()
        theta_x = self.config.theta_x()
        theta_kernel = self._theta_kernel(mu_x, sig2_x, theta_x)
        sig2_weights = _quadrature_weights(sig2_x)
        nb = differential_core.compute_logpmf_values(
            models,
            obs,
            exp,
            float(theta_x[0]),
            float(theta_x[-1]),
            theta_x.size,
        )

        g, n, m, s = len(names), obs.shape[1], mu_x.size, sig2_x.size
        loglik_mu = np.empty((g, n, m))
        loglik_mu_sig2 = np.empty((g, n, m, s))
        log_sig2_prior = np.empty((g, s))
        log_ratio = np.log1p(obs) - np.log1p(exp)

        for i, sample_index in enumerate(indices):
            nu, prior_sig2 = _fit_invchi2(log_ratio[sample_index])
            prior = _sig2_prior(sig2_x, sig2_weights, nu, prior_sig2)
            collapsed, raw = self._integrate_group(
                theta_kernel, prior, nb, sample_index
            )
            loglik_mu[i] = collapsed
            loglik_mu_sig2[i] = raw
            log_sig2_prior[i] = prior

        common_mu = mu_x[np.argmax(loglik_mu.sum(axis=0), axis=1)]
        log_mu_prior = make_group_mean_log_prior(
            mu_x, common_mu, self.config.mu_prior_sd_floor
        )
        return Differential(
            names,
            mu_x,
            sig2_x,
            loglik_mu,
            loglik_mu_sig2,
            log_sig2_prior,
            log_mu_prior,
            theta_x,
        )

    @staticmethod
    def _align(groups, obs, exp, disp_models, selected):
        names = (
            tuple(pd.unique(groups))
            if selected is None
            else tuple(dict.fromkeys(selected))
        )
        groups = groups[groups.isin(names)]
        if len(names) < 2:
            raise ValueError("at least two groups are required")
        index = groups.index
        obs = np.ascontiguousarray(obs.loc[index].to_numpy(float))
        exp = np.ascontiguousarray(exp.loc[index].to_numpy(float))
        models = disp_models.loc[index].to_numpy()
        labels = groups.to_numpy()
        indices = tuple(np.flatnonzero(labels == name) for name in names)
        if any(i.size < 2 for i in indices):
            raise ValueError("every group needs at least two samples")
        return names, indices, obs, exp, models

    def _integrate_group(self, kernel, prior, nb, sample_index):
        s, m, t = kernel.shape
        kernel = kernel.reshape(s * m, t)
        n = nb.shape[1]
        collapsed = np.empty((n, m))
        raw_out = np.empty((n, m, s))
        for lo in range(0, n, self.config.base_chunk_size):
            hi = min(lo + self.config.base_chunk_size, n)
            x = np.take(nb[:, lo:hi], sample_index, axis=2)
            shift = np.max(x, axis=0)
            scaled = np.exp(x - shift[None]).reshape(t, -1)
            with np.errstate(divide="ignore"):
                raw = np.log(kernel @ scaled).reshape(
                    s, m, hi - lo, sample_index.size
                )
            raw = raw.sum(axis=-1) + shift.sum(axis=-1)[None, None]
            raw_out[lo:hi] = raw.transpose(2, 1, 0)
            collapsed[lo:hi] = logsumexp(
                raw + prior[:, None, None], axis=0
            ).T
        return collapsed, raw_out

    @staticmethod
    def _theta_kernel(mu_x, sig2_x, theta_x):
        bounds = np.r_[-np.inf, (theta_x[:-1] + theta_x[1:]) / 2, np.inf]
        z = (
            bounds[None, None] - mu_x[None, :, None]
        ) / np.sqrt(sig2_x)[:, None, None]
        mass = np.maximum(np.diff(ndtr(z), axis=-1), 0.0)
        return np.ascontiguousarray(mass / mass.sum(axis=-1, keepdims=True))


def make_group_mean_log_prior(
    mu_x: np.ndarray,
    common_mu: np.ndarray,
    sd_floor: float | None = None,
) -> np.ndarray:
    """Empirical normal prior for pointwise and segmented group means."""
    floor = np.median(np.diff(mu_x)) if sd_floor is None else sd_floor
    return normal_grid_log_mass(
        mu_x,
        float(np.mean(common_mu)),
        max(float(np.std(common_mu)), float(floor)),
    )


def fit_group_means_segmentation(
    differential: Differential,
    length_prior: LengthPrior,
    config: MeanSegmentationConfig = DEFAULT_MEAN_SEGMENTATION,
    log_mu_prior: np.ndarray | None = None,
) -> Segmentation:
    prior = (
        differential.log_mu_prior
        if log_mu_prior is None
        else normalize_log_mass(log_mu_prior, differential.mu_x.size)
    )
    return segment(
        differential.loglik_mu,
        differential.mu_x,
        differential.group_names,
        length_prior,
        prior,
        config.transition_sd,
        config.forbid_same_state,
    )


def _quadrature_weights(x):
    w = np.empty_like(x)
    w[0] = (x[1] - x[0]) / 2
    w[-1] = (x[-1] - x[-2]) / 2
    w[1:-1] = (x[2:] - x[:-2]) / 2
    return w


def _fit_invchi2(values):
    variance = np.var(values, axis=0, ddof=1)
    variance = variance[np.isfinite(variance) & (variance > 0)]

    def objective(log_parameters):
        nu, sig2 = np.exp(log_parameters)
        return -float(invchi2.log_likelihood(variance, nu, sig2))

    result = scipy.optimize.minimize(
        objective,
        np.log([5.0, np.median(variance)]),
        method="L-BFGS-B",
        bounds=[(-8.0, 12.0), (-20.0, 12.0)],
    )
    return tuple(np.exp(result.x))


def _sig2_prior(sig2_x, weights, nu, sig2):
    log_mass = np.array(
        [invchi2.logpdf(x, nu, sig2) for x in sig2_x]
    ) + np.log(weights)
    return log_mass - logsumexp(log_mass)
