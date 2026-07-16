"""Sample-level latent-signal likelihoods and segmentation."""

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.special import logsumexp, ndtr

from footprint_tools.stats import differential as differential_core

from .config import (
    DEFAULT_THETA,
    DEFAULT_THETA_SEGMENTATION,
    ThetaConfig,
    ThetaMode,
    ThetaSegmentationConfig,
)
from .differential import Differential
from .posterior import GridPosterior, normalize_log_mass
from .segmentation import LengthPrior, Segmentation, segment


@dataclass(frozen=True, slots=True)
class ThetaLikelihood:
    """Per-base evidence over latent sample signals.

    ``loglik`` has shape ``sample x position x theta``.  In ``sample_only``
    mode it is the sample count likelihood ``log p(y_i | theta)``.  In
    ``group_informed`` mode it additionally contains the leave-one-sample-out
    group evidence, divided by the group prior-predictive theta mass.  The
    stored ``log_theta_prior`` is therefore applied once for pointwise
    inference and once per segment during segmentation.
    """

    sample_names: tuple[str, ...]
    sample_groups: tuple[str, ...]
    group_names: tuple[str, ...]
    theta_x: np.ndarray
    loglik: np.ndarray
    log_theta_prior: np.ndarray
    mode: ThetaMode = "group_informed"

    def __post_init__(self) -> None:
        theta_x = np.asarray(self.theta_x, dtype=float)
        loglik = np.asarray(self.loglik)
        prior = np.asarray(self.log_theta_prior, dtype=float)
        n_sample = len(self.sample_names)
        if len(self.sample_groups) != n_sample:
            raise ValueError("sample_groups must match sample_names")
        if loglik.ndim != 3 or loglik.shape[0] != n_sample:
            raise ValueError("loglik must have shape sample x position x theta")
        if loglik.shape[2] != theta_x.size:
            raise ValueError("theta grid does not match the likelihood state axis")
        if prior.shape != (n_sample, theta_x.size):
            raise ValueError("log_theta_prior must have shape sample x theta")
        prior = prior - logsumexp(prior, axis=1, keepdims=True)
        object.__setattr__(self, "theta_x", theta_x)
        object.__setattr__(self, "loglik", np.ascontiguousarray(loglik))
        object.__setattr__(self, "log_theta_prior", np.ascontiguousarray(prior))

    def posterior(self, log_theta_prior: np.ndarray | None = None) -> GridPosterior:
        """Pointwise posterior, with shape sample x position x theta."""
        prior = (
            self.log_theta_prior
            if log_theta_prior is None
            else _normalize_sample_prior(
                log_theta_prior, len(self.sample_names), self.theta_x.size
            )
        )
        return GridPosterior(
            self.theta_x,
            self.loglik + prior[:, None, :],
        )

    def to_npz(self, path) -> None:
        np.savez_compressed(
            path,
            sample_names=np.asarray(self.sample_names, dtype=str),
            sample_groups=np.asarray(self.sample_groups, dtype=str),
            group_names=np.asarray(self.group_names, dtype=str),
            theta_x=self.theta_x,
            loglik=self.loglik,
            log_theta_prior=self.log_theta_prior,
            mode=np.asarray(self.mode),
        )

    @classmethod
    def from_npz(cls, path) -> "ThetaLikelihood":
        with np.load(path, allow_pickle=False) as x:
            return cls(
                tuple(x["sample_names"].tolist()),
                tuple(x["sample_groups"].tolist()),
                tuple(x["group_names"].tolist()),
                x["theta_x"],
                x["loglik"],
                x["log_theta_prior"],
                str(x["mode"].item()),
            )


class ThetaModel:
    def __init__(self, config: ThetaConfig = DEFAULT_THETA):
        self.config = config

    def fit(
        self,
        groups_data: pd.Series,
        obs: pd.DataFrame,
        exp: pd.DataFrame,
        disp_models: pd.Series,
        differential: Differential,
        log_mu_prior: np.ndarray | None = None,
    ) -> ThetaLikelihood:
        return fit_theta_likelihood(
            groups_data,
            obs,
            exp,
            disp_models,
            differential,
            self.config,
            log_mu_prior,
        )


def fit_theta_likelihood(
    groups_data: pd.Series,
    obs: pd.DataFrame,
    exp: pd.DataFrame,
    disp_models: pd.Series,
    differential: Differential,
    config: ThetaConfig = DEFAULT_THETA,
    log_mu_prior: np.ndarray | None = None,
) -> ThetaLikelihood:
    """Construct sample-level theta evidence from the raw count data.

    The target sample is excluded while constructing its group-informed theta
    distribution, so its count likelihood is used exactly once.  Temporary
    sample-by-``(mu, sig2)`` arrays are processed in position chunks and are not
    retained.
    """
    storage_dtype = np.dtype(config.storage_dtype)

    group_names = differential.group_names
    selected = groups_data[groups_data.isin(group_names)]
    if selected.empty:
        raise ValueError("no samples match the differential groups")
    missing = [name for name in group_names if not np.any(selected.to_numpy() == name)]
    if missing:
        raise ValueError(f"groups have no matching samples: {missing}")

    index = selected.index
    sample_names = tuple(str(value) for value in index)
    sample_groups = tuple(str(value) for value in selected.to_numpy())
    obs_array = np.ascontiguousarray(obs.loc[index].to_numpy(dtype=float))
    exp_array = np.ascontiguousarray(exp.loc[index].to_numpy(dtype=float))
    models = disp_models.loc[index].to_numpy()
    labels = selected.to_numpy()
    group_indices = tuple(
        np.flatnonzero(labels == name) for name in group_names
    )

    if differential.theta_x is None:
        raise ValueError(
            "Differential does not contain theta_x; refit it with the current "
            "DifferentialModel before fitting sample-level theta"
        )
    theta_x = np.asarray(differential.theta_x, dtype=float)

    mu_prior = (
        differential.log_mu_prior
        if log_mu_prior is None
        else normalize_log_mass(log_mu_prior, differential.mu_x.size)
    )
    kernel = _theta_kernel(differential.mu_x, differential.sig2_x, theta_x)
    kernel_flat = np.ascontiguousarray(
        kernel.reshape(differential.sig2_x.size * differential.mu_x.size, theta_x.size)
    )
    group_theta_prior = _group_theta_prior(
        kernel_flat,
        mu_prior,
        differential.log_sig2_prior,
    )

    sample_loglik = differential_core.compute_logpmf_values(
        models,
        obs_array,
        exp_array,
        float(theta_x[0]),
        float(theta_x[-1]),
        theta_x.size,
    )
    sample_loglik = np.asarray(sample_loglik, dtype=float)
    expected_shape = (theta_x.size, obs_array.shape[1], obs_array.shape[0])
    if sample_loglik.shape != expected_shape:
        raise ValueError(
            "compute_logpmf_values returned shape "
            f"{sample_loglik.shape}, expected {expected_shape}"
        )

    n_sample = obs_array.shape[0]
    n_position = obs_array.shape[1]
    output = np.empty((n_sample, n_position, theta_x.size), dtype=storage_dtype)
    sample_prior = np.empty((n_sample, theta_x.size), dtype=float)

    for group, indices in enumerate(group_indices):
        sample_prior[indices] = group_theta_prior[group]
        for lo in range(0, n_position, config.position_chunk_size):
            hi = min(lo + config.position_chunk_size, n_position)
            x = np.take(sample_loglik[:, lo:hi], indices, axis=2)
            if config.mode == "sample_only":
                output[indices, lo:hi] = x.transpose(2, 1, 0).astype(
                    storage_dtype, copy=False
                )
                continue

            log_a = _sample_mu_sig2_loglik(
                kernel_flat, x, differential.sig2_x.size, differential.mu_x.size
            )
            group_total = log_a.sum(axis=0)
            prior_joint = (
                differential.log_sig2_prior[group, :, None]
                + mu_prior[None, :]
            )

            for local, sample in enumerate(indices):
                loo = group_total - log_a[local]
                log_predictive = _theta_predictive(
                    kernel_flat,
                    loo + prior_joint[None, :, :],
                )
                corrected = (
                    x[:, :, local].T
                    + log_predictive
                    - group_theta_prior[group][None, :]
                )
                output[sample, lo:hi] = corrected.astype(
                    storage_dtype, copy=False
                )

    return ThetaLikelihood(
        sample_names,
        sample_groups,
        group_names,
        theta_x,
        output,
        sample_prior,
        config.mode,
    )


def fit_theta_segmentation(
    theta: ThetaLikelihood,
    length_prior: LengthPrior,
    config: ThetaSegmentationConfig = DEFAULT_THETA_SEGMENTATION,
    log_theta_prior: np.ndarray | None = None,
) -> Segmentation:
    """Segment each sample's theta track using the generic segmenter."""
    prior = theta.log_theta_prior if log_theta_prior is None else log_theta_prior
    return segment(
        theta.loglik,
        theta.theta_x,
        theta.sample_names,
        length_prior,
        prior,
        config.transition_sd,
        config.forbid_same_state,
    )


def _theta_kernel(mu_x, sig2_x, theta_x):
    bounds = np.r_[-np.inf, (theta_x[:-1] + theta_x[1:]) / 2, np.inf]
    z = (
        bounds[None, None, :] - np.asarray(mu_x)[None, :, None]
    ) / np.sqrt(np.asarray(sig2_x))[:, None, None]
    mass = np.maximum(np.diff(ndtr(z), axis=-1), 0.0)
    mass /= mass.sum(axis=-1, keepdims=True)
    return np.ascontiguousarray(mass)


def _group_theta_prior(kernel_flat, log_mu_prior, log_sig2_prior):
    n_group = log_sig2_prior.shape[0]
    out = np.empty((n_group, kernel_flat.shape[1]))
    for group in range(n_group):
        log_weight = (
            log_sig2_prior[group, :, None] + log_mu_prior[None, :]
        ).reshape(-1)
        shift = np.max(log_weight)
        mass = kernel_flat.T @ np.exp(log_weight - shift)
        mass = np.maximum(mass, np.finfo(float).tiny)
        out[group] = np.log(mass) + shift
        out[group] -= logsumexp(out[group])
    return out


def _sample_mu_sig2_loglik(kernel_flat, sample_loglik, n_sig2, n_mu):
    """Return sample x position x sig2 x mu integrated log likelihoods."""
    n_theta, n_position, n_sample = sample_loglik.shape
    shift = np.max(sample_loglik, axis=0)
    scaled = np.exp(sample_loglik - shift[None]).reshape(
        n_theta, n_position * n_sample
    )
    integrated = np.maximum(kernel_flat @ scaled, np.finfo(float).tiny)
    value = np.log(integrated)
    value = value.reshape(n_sig2, n_mu, n_position, n_sample)
    value += shift[None, None]
    return value.transpose(3, 2, 0, 1)


def _theta_predictive(kernel_flat, log_weight):
    """Normalize p(theta, y_loo) for each position.

    ``log_weight`` has shape position x sig2 x mu, or position x (sig2*mu)
    after flattening.
    """
    n_position = log_weight.shape[0]
    flat = log_weight.reshape(n_position, -1)
    shift = np.max(flat, axis=1)
    scaled = np.exp(flat - shift[:, None]).T
    mass = (kernel_flat.T @ scaled).T
    mass = np.maximum(mass, np.finfo(float).tiny)
    value = np.log(mass) + shift[:, None]
    value -= logsumexp(value, axis=1, keepdims=True)
    return value


def _normalize_sample_prior(value, n_sample, n_state):
    value = np.asarray(value, dtype=float)
    if value.ndim == 1:
        if value.shape != (n_state,):
            raise ValueError("theta prior state axis does not match")
        value = np.broadcast_to(value, (n_sample, n_state)).copy()
    elif value.shape != (n_sample, n_state):
        raise ValueError(
            f"theta prior must have shape ({n_state},) or ({n_sample}, {n_state})"
        )
    value -= logsumexp(value, axis=1, keepdims=True)
    return value
