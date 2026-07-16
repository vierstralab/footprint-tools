from dataclasses import dataclass

import numpy as np
from scipy.special import logsumexp

from .config import (
    DEFAULT_COEFFICIENT,
    DEFAULT_COEFFICIENT_SEGMENTATION,
    CoefficientConfig,
    CoefficientSegmentationConfig,
)
from .differential import Differential
from .integration import coefficient_target_terms, variance_ratio_group_terms
from .posterior import GridPosterior, normal_grid_log_mass, normalize_log_mass
from .segmentation import LengthPrior, Segmentation, segment
from .variance_ratio import VarianceRatioLikelihood


@dataclass(frozen=True, slots=True)
class CoefficientLikelihood:
    group_names: tuple[str, ...]
    z_x: np.ndarray
    loglik: np.ndarray
    log_z_prior: np.ndarray
    reference: str

    def posterior(self) -> GridPosterior:
        return GridPosterior(self.z_x, self.loglik + self.log_z_prior[None, None])

    def to_npz(self, path) -> None:
        np.savez_compressed(
            path,
            group_names=self.group_names,
            z_x=self.z_x,
            loglik=self.loglik,
            log_z_prior=self.log_z_prior,
            reference=np.array(self.reference),
        )

    @classmethod
    def from_npz(cls, path) -> "CoefficientLikelihood":
        with np.load(path, allow_pickle=False) as x:
            if "reference" not in x.files:
                reference = "mu0"
            elif x["reference"].dtype.kind in "fiu":
                reference = "zero" if float(x["reference"]) == 0 else "unknown"
            else:
                reference = str(x["reference"])
            if reference not in ("mu0", "zero"):
                raise ValueError("unknown coefficient reference")
            return cls(
                tuple(x["group_names"].tolist()),
                x["z_x"],
                x["loglik"],
                x["log_z_prior"],
                reference,
            )


class CommonCoefficientModel:
    """Group effects relative to the shared mean: ``(mu_g - mu0) / sigma_g``."""

    def __init__(self, config: CoefficientConfig = DEFAULT_COEFFICIENT):
        self.config = config

    def fit(
        self,
        differential: Differential,
        variance_ratio: VarianceRatioLikelihood,
        log_mu0_prior: np.ndarray | None = None,
        log_eta_prior: np.ndarray | None = None,
    ) -> CoefficientLikelihood:
        if differential.group_names != variance_ratio.group_names:
            raise ValueError("differential and variance-ratio group names differ")

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
        z_x = self.config.z_x()
        log_z_given_eta = _z_given_eta_log_mass(z_x, variance_ratio.eta_x)
        log_z_prior = logsumexp(eta_prior[:, None] + log_z_given_eta, axis=0)
        log_z_prior -= logsumexp(log_z_prior)

        group_terms = variance_ratio_group_terms(
            differential,
            variance_ratio.mu0_x,
            variance_ratio.eta_x,
            variance_ratio.method,
            variance_ratio.variance_floor,
        )
        target = coefficient_target_terms(
            differential,
            variance_ratio.mu0_x,
            z_x,
            variance_ratio.method,
            variance_ratio.variance_floor,
        )
        loglik = _common_reference_loglik(
            group_terms,
            target,
            mu0_prior,
            eta_prior,
            log_z_given_eta,
            log_z_prior,
        )
        return CoefficientLikelihood(
            differential.group_names,
            z_x,
            loglik,
            log_z_prior,
            "mu0",
        )


class ZeroCoefficientModel:
    """Group effects relative to zero: ``mu_g / sigma_g``."""

    def __init__(self, config: CoefficientConfig = DEFAULT_COEFFICIENT):
        self.config = config

    def fit(
        self,
        differential: Differential,
        log_z_prior: np.ndarray | None = None,
    ) -> CoefficientLikelihood:
        z_x = self.config.z_x()
        prior = (
            normal_grid_log_mass(z_x, 0.0, self.config.z_prior_sd)
            if log_z_prior is None
            else normalize_log_mass(log_z_prior, z_x.size)
        )
        loglik = coefficient_target_terms(
            differential,
            np.array([0.0]),
            z_x,
            self.config.method,
            self.config.variance_floor,
        )[:, :, 0]
        return CoefficientLikelihood(
            differential.group_names,
            z_x,
            loglik,
            prior,
            "zero",
        )


def _z_given_eta_log_mass(z_x: np.ndarray, eta_x: np.ndarray) -> np.ndarray:
    zero = np.flatnonzero(z_x == 0.0)
    if zero.size != 1:
        raise ValueError("z grid must contain exactly one zero state")
    out = np.full((eta_x.size, z_x.size), -np.inf)
    for i, eta in enumerate(eta_x):
        if np.isneginf(eta):
            out[i, zero[0]] = 0.0
        else:
            out[i] = normal_grid_log_mass(z_x, 0.0, np.exp(eta / 2))
    return out


def _common_reference_loglik(
    group_terms: np.ndarray,
    target: np.ndarray,
    mu0_prior: np.ndarray,
    eta_prior: np.ndarray,
    log_z_given_eta: np.ndarray,
    log_z_prior: np.ndarray,
) -> np.ndarray:
    finite = np.isfinite(group_terms)
    total = np.where(finite, group_terms, 0.0).sum(axis=0)
    invalid = (~finite).sum(axis=0)
    z_given_eta = np.exp(log_z_given_eta)
    out = np.empty(target.shape[:2] + (target.shape[-1],))

    for g in range(group_terms.shape[0]):
        other = total - np.where(finite[g], group_terms[g], 0.0)
        other = np.where(invalid - (~finite[g]) == 0, other, -np.inf)
        reference = _mix_eta(other + eta_prior[None, None], z_given_eta)
        joint = logsumexp(
            reference + target[g] + mu0_prior[None, :, None],
            axis=1,
        )
        out[g] = joint - log_z_prior[None]
    return out


def _mix_eta(log_weight: np.ndarray, conditional_mass: np.ndarray) -> np.ndarray:
    maximum = np.max(log_weight, axis=-1)
    valid = np.isfinite(maximum)
    scaled = np.zeros_like(log_weight)
    scaled[valid] = np.exp(log_weight[valid] - maximum[valid, None])
    mass = scaled.reshape(-1, scaled.shape[-1]) @ conditional_mass
    with np.errstate(divide="ignore"):
        out = np.log(mass).reshape(maximum.shape + (conditional_mass.shape[1],))
    out[valid] += maximum[valid, None]
    out[~valid] = -np.inf
    return out


def fit_coefficient_segmentation(
    likelihood: CoefficientLikelihood,
    length_prior: LengthPrior,
    config: CoefficientSegmentationConfig = DEFAULT_COEFFICIENT_SEGMENTATION,
) -> Segmentation:
    return segment(
        likelihood.loglik,
        likelihood.z_x,
        likelihood.group_names,
        length_prior,
        likelihood.log_z_prior,
        config.transition_sd,
        config.forbid_same_state,
    )
