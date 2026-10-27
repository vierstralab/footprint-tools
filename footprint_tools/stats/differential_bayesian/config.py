import json
from dataclasses import asdict, dataclass, fields, replace
from typing import Literal

import numpy as np


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
class ThetaSegmentationConfig:
    transition_sd: float | None = None
    forbid_same_state: bool = False


@dataclass(frozen=True, slots=True)
class NestedConfig:
    """Global constants of the nested segmentation model (log2 units)."""

    tau2: float = 0.4
    kappa: float = 2.2
    tau2_contrast: float = 0.0
    omega_scale: float = 0.15
    n_omega: int = 25
    omega2_min: float = 1e-4
    omega2_max: float = 4.0
    b_var: float = 0.5
    inner_mean: float = 20.0
    outer_mean: float = 70.0
    len_shape: float = 3.0
    len_min: int = 2
    l_max: int = 150
    sigma2_bar: float = 0.14
    n_draws: int = 1000
    seed: int = 0
    prob: float = 0.95
    delta: float = 0.5
    store_draws: bool = False

    def __post_init__(self) -> None:
        positive = ("tau2", "kappa", "omega_scale", "b_var", "inner_mean", "outer_mean", "len_shape", "sigma2_bar")
        for name in positive:
            if getattr(self, name) <= 0:
                raise ValueError(f"{name} must be positive")
        if self.tau2_contrast < 0:
            raise ValueError("tau2_contrast must be non-negative")
        if not 1 <= self.len_min <= self.l_max:
            raise ValueError("require 1 <= len_min <= l_max")
        if not 0 < self.prob < 1:
            raise ValueError("prob must lie in (0, 1)")

    def omega2(self) -> np.ndarray:
        return np.geomspace(self.omega2_min, self.omega2_max, self.n_omega)

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, values: dict) -> "NestedConfig":
        known = {f.name for f in fields(cls)}
        return cls(**{k: v for k, v in values.items() if k in known})

    @classmethod
    def from_json(cls, path) -> "NestedConfig":
        with open(path) as handle:
            return cls.from_dict(json.load(handle))

    def to_json(self, path) -> None:
        with open(path, "w") as handle:
            json.dump(self.to_dict(), handle, indent=2)

    def with_(self, **changes) -> "NestedConfig":
        return replace(self, **changes)


DEFAULT_DIFFERENTIAL = DifferentialConfig()
DEFAULT_THETA = ThetaConfig()
DEFAULT_THETA_SEGMENTATION = ThetaSegmentationConfig()
DEFAULT_NESTED = NestedConfig()
