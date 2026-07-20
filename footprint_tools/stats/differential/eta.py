from dataclasses import dataclass, field
from typing import ClassVar

import numpy as np
from scipy.special import logsumexp

from .config import DEFAULT_ETA_SEGMENTATION, EtaSegmentationConfig
from .io import Serializable, tuple_str
from .posterior import GridPosterior, normalize_log_mass
from .segmentation import LengthPrior, Segmentation, segment
from .variance_ratio import VarianceRatioLikelihood


@dataclass(frozen=True, slots=True)
class EtaSegmentation(Serializable):
    save_attrs: ClassVar[tuple[str, ...]] = (
        "group_names", "eta_x", "mu0_x", "log_mu0", "icc_x",
        "log_icc", "boundary", "log_partition"
    )

    group_names: tuple[str, ...]
    eta_x: np.ndarray
    mu0: GridPosterior
    icc: GridPosterior
    boundary: np.ndarray
    log_partition: float
    _base: Segmentation | None = field(default=None, repr=False, compare=False)

    def sample(self, n_draws=1, rng=None):
        if self._base is None:
            raise RuntimeError(
                "sampling is unavailable after loading a summary-only NPZ"
            )
        return self._base.sample(n_draws, rng)

    def sample_prior(self, n_draws=1, rng=None):
        if self._base is None:
            raise RuntimeError(
                "sampling is unavailable after loading a summary-only NPZ"
            )
        return self._base.sample_prior(n_draws, rng)

    def to_dict(self) -> dict[str, object]:
        return {
            "group_names": np.asarray(self.group_names, dtype=str),
            "eta_x": self.eta_x,
            "mu0_x": self.mu0.x,
            "log_mu0": self.mu0.log_mass,
            "icc_x": self.icc.x,
            "log_icc": self.icc.log_mass,
            "boundary": self.boundary,
            "log_partition": np.asarray(self.log_partition),
        }

    @classmethod
    def from_dict(cls, data: dict[str, object]) -> "EtaSegmentation":
        return cls(
            tuple_str(data["group_names"]),
            data["eta_x"],
            GridPosterior(data["mu0_x"], data["log_mu0"]),
            GridPosterior(data["icc_x"], data["log_icc"]),
            data["boundary"],
            float(np.asarray(data["log_partition"]).item()),
        )


def fit_eta_segmentation(
    likelihood: VarianceRatioLikelihood,
    length_prior: LengthPrior,
    config: EtaSegmentationConfig = DEFAULT_ETA_SEGMENTATION,
    log_mu0_prior: np.ndarray | None = None,
    log_eta_prior: np.ndarray | None = None,
) -> EtaSegmentation:
    mu0_prior = (
        likelihood.log_mu0_prior
        if log_mu0_prior is None
        else normalize_log_mass(log_mu0_prior, likelihood.mu0_x.size)
    )
    eta_prior = (
        likelihood.log_eta_prior
        if log_eta_prior is None
        else normalize_log_mass(log_eta_prior, likelihood.eta_x.size)
    )
    emission = logsumexp(
        likelihood.loglik + mu0_prior[None, :, None], axis=1
    )
    base = segment(
        emission,
        likelihood.icc_x,
        ("icc",),
        length_prior,
        eta_prior,
        config.transition_sd,
        config.forbid_same_state,
    )
    log_icc = base.posterior.log_mass[0]
    log_mu0 = _reconstruct_mu0(
        likelihood.loglik,
        mu0_prior,
        emission,
        log_icc,
    )
    return EtaSegmentation(
        likelihood.group_names,
        likelihood.eta_x,
        GridPosterior(likelihood.mu0_x, log_mu0),
        GridPosterior(likelihood.icc_x, log_icc),
        base.boundary[0],
        float(base.log_partition[0]),
        base,
    )


def _reconstruct_mu0(loglik, log_mu0_prior, emission, log_state):
    conditional = (
        loglik + log_mu0_prior[None, :, None] - emission[:, None, :]
    )
    return logsumexp(log_state[:, None, :] + conditional, axis=2)
