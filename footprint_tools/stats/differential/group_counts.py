"""Posterior summaries of footprinted and deviant groups."""

from dataclasses import dataclass

import numpy as np
from scipy.special import logsumexp, ndtr

from .differential import Differential
from .integration import gaussian_summary
from .posterior import GridPosterior, normalize_log_mass
from .variance_ratio import VarianceRatioLikelihood


@dataclass(frozen=True, slots=True)
class SoftFootprintCount:
    """Exact posterior moments of predictive group footprint prevalence."""

    group_names: tuple[str, ...]
    prevalence: np.ndarray
    prevalence_variance: np.ndarray
    threshold: float

    @property
    def mean(self) -> np.ndarray:
        return self.prevalence.sum(axis=0)

    @property
    def expected_count(self) -> np.ndarray:
        return self.mean

    @property
    def strongest(self) -> np.ndarray:
        return self.prevalence.max(axis=0)


def infer_ksoft(
    differential: Differential,
    threshold: float = 0.0,
    log_mu_prior: np.ndarray | None = None,
    position_chunk_size: int = 64,
) -> SoftFootprintCount:
    """Exact moments of ``P(theta_new < threshold | data)`` by group."""
    if position_chunk_size < 1:
        raise ValueError("position_chunk_size must be positive")
    prior = (
        differential.log_mu_prior
        if log_mu_prior is None
        else normalize_log_mass(log_mu_prior, differential.mu_x.size)
    )
    kernel = ndtr(
        (threshold - differential.mu_x[:, None])
        / np.sqrt(differential.sig2_x[None])
    )
    prevalence = np.empty(differential.loglik_mu.shape[:2])
    variance = np.empty_like(prevalence)

    for lo in range(0, prevalence.shape[1], position_chunk_size):
        hi = min(lo + position_chunk_size, prevalence.shape[1])
        joint = (
            differential.loglik_mu_sig2[:, lo:hi]
            + prior[None, None, :, None]
            + differential.log_sig2_prior[:, None, None]
        )
        denominator = logsumexp(joint, axis=(2, 3))
        mean = np.exp(
            logsumexp(joint, axis=(2, 3), b=kernel[None, None])
            - denominator
        )
        second = np.exp(
            logsumexp(joint, axis=(2, 3), b=(kernel * kernel)[None, None])
            - denominator
        )
        prevalence[:, lo:hi] = mean
        variance[:, lo:hi] = np.maximum(second - mean * mean, 0.0)

    return SoftFootprintCount(
        differential.group_names, prevalence, variance, threshold
    )



def binary_count_log_evidence(
    log_evidence_0: np.ndarray,
    log_evidence_1: np.ndarray,
) -> np.ndarray:
    """Combine independent binary contributions into evidence over event count.

    The first axis indexes groups. All remaining axes are arbitrary conditioning
    dimensions. The returned final axis indexes counts from zero through the
    number of groups.
    """
    log_evidence_0 = np.asarray(log_evidence_0, dtype=float)
    log_evidence_1 = np.asarray(log_evidence_1, dtype=float)
    if log_evidence_0.shape != log_evidence_1.shape:
        raise ValueError("binary evidence arrays must have identical shapes")
    if log_evidence_0.ndim < 1 or log_evidence_0.shape[0] < 1:
        raise ValueError("the first axis must contain at least one group")

    n_group = log_evidence_0.shape[0]
    out = np.full(log_evidence_0.shape[1:] + (n_group + 1,), -np.inf)
    out[..., 0] = 0.0

    for group in range(n_group):
        old = out.copy()
        out.fill(-np.inf)
        log0 = log_evidence_0[group][..., None]
        log1 = log_evidence_1[group][..., None]
        out[..., : group + 1] = old[..., : group + 1] + log0
        out[..., 1 : group + 2] = np.logaddexp(
            out[..., 1 : group + 2],
            old[..., : group + 1] + log1,
        )
    return out


def _count_posterior(probability: np.ndarray) -> GridPosterior:
    probability = np.clip(np.asarray(probability, dtype=float), 0.0, 1.0)
    with np.errstate(divide="ignore"):
        log_event = np.log(probability)
        log_no_event = np.log1p(-probability)
    return GridPosterior(
        np.arange(probability.shape[0] + 1, dtype=float),
        binary_count_log_evidence(log_no_event, log_event),
    )


def infer_kfp(
    mu_posterior: GridPosterior,
    threshold: float = 0.0,
) -> GridPosterior:
    """Posterior count of groups with ``mu < threshold`` at each position."""
    if mu_posterior.log_mass.ndim != 3:
        raise ValueError("mu_posterior must have shape group x position x mu")
    return _count_posterior(mu_posterior.prob_lt(threshold))


def infer_kfp_zero(
    zero_z_posterior: GridPosterior,
    threshold: float = 1.0,
) -> GridPosterior:
    """Posterior count of groups with ``mu / sigma < -threshold``."""
    if threshold < 0:
        raise ValueError("threshold must be nonnegative")
    if zero_z_posterior.log_mass.ndim != 3:
        raise ValueError(
            "zero_z_posterior must have shape group x position x z"
        )
    return _count_posterior(zero_z_posterior.prob_lt(-threshold))


def infer_kdev(
    differential: Differential,
    variance_ratio: VarianceRatioLikelihood,
    threshold: float = 1.0,
    variance_floor: float | None = None,
    position_chunk_size: int = 64,
    log_mu0_prior: np.ndarray | None = None,
    log_eta_prior: np.ndarray | None = None,
) -> GridPosterior:
    """Gaussian posterior count of groups with ``|d_g| > threshold``.

    Here ``d_g = (mu_g - mu0) / sigma_g`` uses the common ``mu0`` and variance
    ratio from the hierarchical variance-ratio model. The posterior count is
    obtained conditionally on ``(mu0, eta)`` and then those shared nuisance
    states are integrated out.

    ``variance_floor`` should match the value used to construct a Gaussian
    ``VarianceRatioLikelihood`` when a non-default floor was supplied.
    """
    if threshold < 0:
        raise ValueError("threshold must be nonnegative")
    if position_chunk_size < 1:
        raise ValueError("position_chunk_size must be positive")
    if differential.group_names != variance_ratio.group_names:
        raise ValueError("differential and variance-ratio group names differ")
    if differential.loglik_mu.shape[1] != variance_ratio.loglik.shape[0]:
        raise ValueError("differential and variance-ratio position counts differ")

    mu0_prior = (
        variance_ratio.log_mu0_prior
        if log_mu0_prior is None
        else normalize_log_mass(log_mu0_prior, variance_ratio.mu0_x.size)
    )
    eta_prior = (
        variance_ratio.log_eta_prior
        if log_eta_prior is None
        else normalize_log_mass(log_eta_prior, variance_ratio.eta_x.size)
    )

    if variance_floor is None:
        variance_floor = variance_ratio.variance_floor
    summary = gaussian_summary(differential, variance_floor)
    ratio_x = variance_ratio.ratio_x
    n_group, n_position = summary.mu.shape
    log_count = np.empty((n_position, n_group + 1))

    for lo in range(0, n_position, position_chunk_size):
        hi = min(lo + position_chunk_size, n_position)

        mu = summary.mu[:, lo:hi, None]
        mu_var = summary.mu_var[:, lo:hi, None]
        sig2 = summary.sig2[:, lo:hi, None]
        sigma = np.sqrt(sig2)
        evidence = summary.log_evidence[:, lo:hi, None]
        global_mu = variance_ratio.mu0_x[None, None, :]
        delta = mu - global_mu
        chunk_count = np.full((hi - lo, n_group + 1), -np.inf)

        for eta_index, ratio in enumerate(ratio_x):
            marginal_var = mu_var + ratio * sig2
            log_group_evidence = evidence - 0.5 * (
                np.log(2 * np.pi * marginal_var)
                + delta * delta / marginal_var
            )

            if ratio == 0:
                probability_0 = np.ones_like(log_group_evidence)
            else:
                posterior_mean = ratio * sigma * delta / marginal_var
                posterior_sd = np.sqrt(ratio * mu_var / marginal_var)
                lower = (-threshold - posterior_mean) / posterior_sd
                upper = (threshold - posterior_mean) / posterior_sd
                probability_0 = ndtr(upper) - ndtr(lower)
                probability_0 = np.clip(probability_0, 0.0, 1.0)

            probability_1 = np.clip(1.0 - probability_0, 0.0, 1.0)
            with np.errstate(divide="ignore"):
                log_0 = log_group_evidence + np.log(probability_0)
                log_1 = log_group_evidence + np.log(probability_1)

            conditional_count = binary_count_log_evidence(log_0, log_1)
            integrated_mu0 = logsumexp(
                conditional_count + mu0_prior[None, :, None],
                axis=1,
            )
            chunk_count = np.logaddexp(
                chunk_count,
                integrated_mu0 + eta_prior[eta_index],
            )

        log_count[lo:hi] = chunk_count

    return GridPosterior(
        np.arange(n_group + 1, dtype=float),
        log_count,
    )
