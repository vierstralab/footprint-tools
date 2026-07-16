from dataclasses import dataclass

import numpy as np
from scipy.special import logsumexp, ndtr


@dataclass(frozen=True, slots=True)
class GridPosterior:
    """Normalized posterior on a one-dimensional grid; state is the last axis."""

    x: np.ndarray
    log_mass: np.ndarray

    def __post_init__(self) -> None:
        x = np.asarray(self.x, dtype=float)
        log_mass = np.asarray(self.log_mass, dtype=float)
        if log_mass.shape[-1] != x.size:
            raise ValueError("posterior state axis does not match the grid")
        log_mass = log_mass - logsumexp(log_mass, axis=-1, keepdims=True)
        object.__setattr__(self, "x", x)
        object.__setattr__(self, "log_mass", np.ascontiguousarray(log_mass))

    @classmethod
    def from_log_weight(cls, x: np.ndarray, log_weight: np.ndarray) -> "GridPosterior":
        return cls(x, log_weight)

    @property
    def mass(self) -> np.ndarray:
        return np.exp(self.log_mass)

    @property
    def mean(self) -> np.ndarray:
        return self.mass @ self.x

    @property
    def map(self) -> np.ndarray:
        return self.x[np.argmax(self.log_mass, axis=-1)]

    def quantile(self, probability: float) -> np.ndarray:
        if not 0 <= probability <= 1:
            raise ValueError("probability must lie in [0, 1]")
        cumulative = np.cumsum(self.mass, axis=-1)
        return self.x[np.argmax(cumulative >= probability, axis=-1)]

    def prob_gt(self, threshold: float) -> np.ndarray:
        return self._prob(self.x > threshold)

    def prob_ge(self, threshold: float) -> np.ndarray:
        return self._prob(self.x >= threshold)

    def prob_lt(self, threshold: float) -> np.ndarray:
        return self._prob(self.x < threshold)

    def prob_le(self, threshold: float) -> np.ndarray:
        return self._prob(self.x <= threshold)

    def prob_abs_gt(self, threshold: float) -> np.ndarray:
        return self._prob(np.abs(self.x) > threshold)

    def prob_abs_lt(self, threshold: float) -> np.ndarray:
        return self._prob(np.abs(self.x) < threshold)

    def with_grid(self, x: np.ndarray) -> "GridPosterior":
        return GridPosterior(x, self.log_mass)

    def _prob(self, mask: np.ndarray) -> np.ndarray:
        if not np.any(mask):
            return np.zeros(self.log_mass.shape[:-1])
        return np.exp(logsumexp(self.log_mass[..., mask], axis=-1))


def normalize_log_mass(value: np.ndarray, size: int) -> np.ndarray:
    value = np.asarray(value, dtype=float)
    if value.shape != (size,):
        raise ValueError(f"prior must have shape ({size},)")
    return value - logsumexp(value)


def normal_grid_log_mass(x: np.ndarray, mean: float, sd: float) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    if x.size == 1:
        return np.zeros(1)
    bounds = np.r_[-np.inf, (x[:-1] + x[1:]) / 2, np.inf]
    mass = np.maximum(np.diff(ndtr((bounds - mean) / sd)), np.finfo(float).tiny)
    return np.log(mass) - np.log(mass.sum())


def spike_slab_log_mass(
    x: np.ndarray,
    spike_index: int,
    spike_mass: float,
    slab_mean: float,
    slab_sd: float,
) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    finite = np.isfinite(x)
    slab = np.zeros(x.size)
    finite_x = x[finite]
    finite_mass = np.exp(normal_grid_log_mass(finite_x, slab_mean, slab_sd))
    slab[finite] = finite_mass
    slab[spike_index] = 0.0
    slab /= slab.sum()

    mass = (1 - spike_mass) * slab
    mass[spike_index] = spike_mass
    with np.errstate(divide="ignore"):
        return np.log(mass)
