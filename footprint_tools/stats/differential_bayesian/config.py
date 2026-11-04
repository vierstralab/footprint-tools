import json
from collections.abc import Sequence
from dataclasses import asdict, dataclass, field, replace

import numpy as np


@dataclass(frozen=True)
class DifferentialFootprintSummaryConfig:
    """Grids for the negative-binomial group summaries (log2 units)."""

    mu_max: float = 12.0
    n_mu: int = 401
    sig2_min: float = 1e-3
    sig2_max: float = 10.0
    n_sig2: int = 81
    theta_step: float = 0.1
    chunk: int = 8

    def grids(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        mu = np.linspace(-self.mu_max, self.mu_max, self.n_mu)
        sig2 = np.geomspace(self.sig2_min, self.sig2_max, self.n_sig2)
        half = self.mu_max + 5 * np.sqrt(self.sig2_max)
        theta = np.linspace(-half, half, int(np.ceil(2 * half / self.theta_step)) + 1)
        return mu, sig2, theta


@dataclass(frozen=True)
class DifferentialFootprintModelConfig:
    """Global constants of the two-layer segmentation model and the call rule (log2 units)."""

    tau2: float = 0.6               # shared nugget: base-level jitter common to all groups
    b_var: float = 1.0              # prior variance of the consensus level
    omega_scale: float = 0.2        # half-Cauchy scale of group deviations
    group_nugget: dict = field(default_factory=dict)   # t_g by group name (moments, hyper.py)
    cons_mean: float = 24.8         # consensus block lengths: discretized gamma (footprint index)
    cons_shape: float = 1.9
    cons_min: int = 2
    pattern_mean: float = 20.0      # group-pattern block lengths (resolution choice)
    pattern_shape: float = 3.0
    pattern_min: int = 2
    l_max: int = 150
    n_draws: int = 1000
    seed: int = 0
    prob: float = 0.95              # posterior probability of a call
    delta: float = 0.585            # |Delta| cutoff, log2 1.5
    call_floor: float = 0.1         # c: groups are compared on log2(c + (1 - c) 2^depth)
    min_len: int = 6                # shortest reported call

    def nugget(self, group_names: Sequence[str]) -> np.ndarray:
        return np.array([self.group_nugget.get(str(g), 0.0) for g in group_names])

    def with_(self, **changes) -> "DifferentialFootprintModelConfig":
        return replace(self, **changes)

    def to_json(self, path: str | None = None) -> str | None:
        text = json.dumps(asdict(self))
        if path is None:
            return text
        with open(path, "w") as handle:
            handle.write(text)

    @classmethod
    def from_json(cls, text_or_path: str) -> "DifferentialFootprintModelConfig":
        text = text_or_path if text_or_path.lstrip().startswith("{") else open(text_or_path).read()
        return cls(**json.loads(text))


DEFAULT_DF_SUMMARY = DifferentialFootprintSummaryConfig()
DEFAULT_DF_MODEL = DifferentialFootprintModelConfig()
