from dataclasses import dataclass
from typing import ClassVar

import numpy as np
from scipy.special import logsumexp

from .config import (
    DEFAULT_MEAN_SEGMENTATION,
    DEFAULT_VARIANCE_RATIO,
    MeanSegmentationConfig,
    VarianceRatioConfig,
)
from .differential import Differential
from .integration import variance_ratio_group_terms
from .io import Serializable, optional_float, tuple_str
from .posterior import (
    GridPosterior,
    normal_grid_log_mass,
    normalize_log_mass,
    spike_slab_log_mass,
)
from .segmentation import LengthPrior, Segmentation, segment


@dataclass(frozen=True, slots=True)
class VarianceRatioPosterior:
    mu0: GridPosterior
    icc: GridPosterior


@dataclass(frozen=True, slots=True)
class VarianceRatioLikelihood(Serializable):
    save_attrs: ClassVar[tuple[str, ...]] = (
        "group_names", "mu0_x", "eta_x", "loglik", "log_mu0_prior",
        "log_eta_prior", "method", "variance_floor"
    )

    group_names: tuple[str, ...]
    mu0_x: np.ndarray
    eta_x: np.ndarray
    loglik: np.ndarray
    log_mu0_prior: np.ndarray
    log_eta_prior: np.ndarray
    method: str = "gaussian"
    variance_floor: float | None = None

    @property
    def ratio_x(self) -> np.ndarray:
        return np.where(np.isneginf(self.eta_x), 0.0, np.exp(self.eta_x))

    @property
    def icc_x(self) -> np.ndarray:
        ratio = self.ratio_x
        return ratio / (1 + ratio)

    def marginal_loglik_mu0(
        self,
        log_eta_prior: np.ndarray | None = None,
    ) -> np.ndarray:
        """Per-base likelihood over mu0 after integrating eta."""
        eta_prior = (
            self.log_eta_prior
            if log_eta_prior is None
            else normalize_log_mass(log_eta_prior, self.eta_x.size)
        )
        return logsumexp(
            self.loglik + eta_prior[None, None, :],
            axis=2,
        )

    def posterior(
        self,
        log_mu0_prior: np.ndarray | None = None,
        log_eta_prior: np.ndarray | None = None,
    ) -> VarianceRatioPosterior:
        mu0_prior = (
            self.log_mu0_prior
            if log_mu0_prior is None
            else normalize_log_mass(log_mu0_prior, self.mu0_x.size)
        )
        eta_prior = (
            self.log_eta_prior
            if log_eta_prior is None
            else normalize_log_mass(log_eta_prior, self.eta_x.size)
        )
        joint = (
            self.loglik
            + mu0_prior[None, :, None]
            + eta_prior[None, None, :]
        )
        joint -= logsumexp(joint, axis=(1, 2), keepdims=True)
        return VarianceRatioPosterior(
            GridPosterior(self.mu0_x, logsumexp(joint, axis=2)),
            GridPosterior(self.icc_x, logsumexp(joint, axis=1)),
        )

    def to_dict(self) -> dict[str, object]:
        out = Serializable.to_dict(self)
        out["group_names"] = np.asarray(self.group_names, dtype=str)
        out["method"] = np.asarray(self.method)
        out["variance_floor"] = np.asarray(
            np.nan if self.variance_floor is None else self.variance_floor
        )
        return out

    @classmethod
    def from_dict(cls, data: dict[str, object]) -> "VarianceRatioLikelihood":
        return cls(
            tuple_str(data["group_names"]),
            data["mu0_x"],
            data["eta_x"],
            data["loglik"],
            data["log_mu0_prior"],
            data["log_eta_prior"],
            str(np.asarray(data.get("method", "gaussian")).item()),
            optional_float(data.get("variance_floor")),
        )


class VarianceRatioModel:
    def __init__(self, config: VarianceRatioConfig = DEFAULT_VARIANCE_RATIO):
        self.config = config

    def fit(
        self,
        differential: Differential,
        log_mu0_prior: np.ndarray | None = None,
        log_eta_prior: np.ndarray | None = None,
    ) -> VarianceRatioLikelihood:
        mu0_x = differential.mu_x
        eta_x = self.config.eta_x()
        mu0_prior = (
            differential.log_mu_prior
            if log_mu0_prior is None
            else normalize_log_mass(log_mu0_prior, mu0_x.size)
        )
        eta_prior = (
            make_eta_log_prior(eta_x, self.config)
            if log_eta_prior is None
            else normalize_log_mass(log_eta_prior, eta_x.size)
        )
        group_terms = variance_ratio_group_terms(
            differential,
            mu0_x,
            eta_x,
            self.config.method,
            self.config.variance_floor,
        )
        return VarianceRatioLikelihood(
            differential.group_names,
            mu0_x,
            eta_x,
            group_terms.sum(axis=0),
            mu0_prior,
            eta_prior,
            self.config.method,
            self.config.variance_floor,
        )


def make_eta_log_prior(
    eta_x: np.ndarray,
    config: VarianceRatioConfig = DEFAULT_VARIANCE_RATIO,
) -> np.ndarray:
    spike = np.flatnonzero(np.isneginf(eta_x))
    if spike.size:
        if config.consistent_mass is None:
            raise ValueError(
                "eta_x contains an exact zero state but consistent_mass is None"
            )
        return spike_slab_log_mass(
            eta_x,
            int(spike[0]),
            config.consistent_mass,
            config.eta_prior_mean,
            config.eta_prior_sd,
        )
    return normal_grid_log_mass(
        eta_x,
        config.eta_prior_mean,
        config.eta_prior_sd,
    )


def fit_mu0_segmentation(
    likelihood: VarianceRatioLikelihood,
    length_prior: LengthPrior,
    config: MeanSegmentationConfig = DEFAULT_MEAN_SEGMENTATION,
    log_mu0_prior: np.ndarray | None = None,
    log_eta_prior: np.ndarray | None = None,
) -> Segmentation:
    """Segment the global mean after integrating eta independently per base."""
    mu0_prior = (
        likelihood.log_mu0_prior
        if log_mu0_prior is None
        else normalize_log_mass(log_mu0_prior, likelihood.mu0_x.size)
    )
    return segment(
        likelihood.marginal_loglik_mu0(log_eta_prior),
        likelihood.mu0_x,
        ("mu0",),
        length_prior,
        mu0_prior,
        config.transition_sd,
        config.forbid_same_state,
    )
