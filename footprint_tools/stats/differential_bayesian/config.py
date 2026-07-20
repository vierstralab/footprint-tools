from dataclasses import dataclass
from typing import Literal

import numpy as np

IntegrationMethod = Literal["gaussian", "exact"]


ThetaMode = Literal["sample_only", "group_informed"]


@dataclass(frozen=True, slots=True)
class ThetaConfig:
    mode: ThetaMode = "group_informed"
    position_chunk_size: int = 8
    storage_dtype: str = "float32"

    def __post_init__(self) -> None:
        if self.mode not in ("sample_only", "group_informed"):
            raise ValueError("mode must be 'sample_only' or 'group_informed'")
        if self.position_chunk_size < 1:
            raise ValueError("position_chunk_size must be positive")
        if np.dtype(self.storage_dtype).kind != "f":
            raise ValueError("storage_dtype must be a floating-point dtype")


@dataclass(frozen=True, slots=True)
class DifferentialConfig:
    mu_min: float = -6.0
    mu_max: float = 6.0
    n_mu: int = 481
    sig2_min: float = 1e-3
    sig2_max: float = 10.0
    n_sig2: int = 181
    theta_step: float = 0.1
    theta_tail_sd: float = 5.0
    theta_grid: tuple[float, ...] | None = None
    base_chunk_size: int = 32
    mu_prior_sd_floor: float | None = None

    def mu_x(self) -> np.ndarray:
        return np.linspace(self.mu_min, self.mu_max, self.n_mu)

    def sig2_x(self) -> np.ndarray:
        return np.geomspace(self.sig2_min, self.sig2_max, self.n_sig2)

    def theta_x(self) -> np.ndarray:
        if self.theta_grid is not None:
            x = np.asarray(self.theta_grid, dtype=float)
            if x.ndim != 1 or x.size < 2 or np.any(np.diff(x) <= 0):
                raise ValueError(
                    "theta_grid must be a strictly increasing one-dimensional grid"
                )
            return x
        margin = self.theta_tail_sd * np.sqrt(self.sig2_max)
        lo, hi = self.mu_min - margin, self.mu_max + margin
        n = int(np.ceil((hi - lo) / self.theta_step)) + 1
        return np.linspace(lo, hi, n)


@dataclass(frozen=True, slots=True)
class VarianceRatioConfig:
    eta_min: float = -4.0
    eta_max: float = 4.0
    eta_step: float = 0.25
    include_zero: bool = True
    method: IntegrationMethod = "gaussian"
    variance_floor: float | None = None
    consistent_mass: float | None = 0.30
    eta_prior_mean: float = -1.0
    eta_prior_sd: float = 2.0

    def eta_x(self) -> np.ndarray:
        finite = np.arange(
            self.eta_min,
            self.eta_max + self.eta_step / 2,
            self.eta_step,
        )
        if self.include_zero and self.consistent_mass is not None:
            return np.r_[-np.inf, finite]
        return finite


@dataclass(frozen=True, slots=True)
class CoefficientConfig:
    z_min: float = -4.0
    z_max: float = 4.0
    z_step: float = 0.25
    z_prior_sd: float = 2.0
    method: IntegrationMethod = "exact"
    variance_floor: float | None = None

    def __post_init__(self) -> None:
        if self.z_prior_sd <= 0:
            raise ValueError("z_prior_sd must be positive")
        if self.method not in ("gaussian", "exact"):
            raise ValueError("method must be 'gaussian' or 'exact'")

    def z_x(self) -> np.ndarray:
        x = np.arange(
            self.z_min,
            self.z_max + self.z_step / 2,
            self.z_step,
        )
        x[np.argmin(np.abs(x))] = 0.0
        return x


@dataclass(frozen=True, slots=True)
class ThetaSegmentationConfig:
    transition_sd: float | None = None
    forbid_same_state: bool = False


@dataclass(frozen=True, slots=True)
class EtaSegmentationConfig:
    transition_sd: float | None = None
    forbid_same_state: bool = True


@dataclass(frozen=True, slots=True)
class CoefficientSegmentationConfig:
    transition_sd: float | None = None
    forbid_same_state: bool = True


@dataclass(frozen=True, slots=True)
class MeanSegmentationConfig:
    transition_sd: float | None = None
    forbid_same_state: bool = False


DEFAULT_DIFFERENTIAL = DifferentialConfig()
DEFAULT_THETA = ThetaConfig()
DEFAULT_VARIANCE_RATIO = VarianceRatioConfig()
DEFAULT_COEFFICIENT = CoefficientConfig()
DEFAULT_THETA_SEGMENTATION = ThetaSegmentationConfig()
DEFAULT_ETA_SEGMENTATION = EtaSegmentationConfig()
DEFAULT_COEFFICIENT_SEGMENTATION = CoefficientSegmentationConfig()
DEFAULT_MEAN_SEGMENTATION = MeanSegmentationConfig()
