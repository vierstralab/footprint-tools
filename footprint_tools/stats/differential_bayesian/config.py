import json
from dataclasses import asdict, dataclass, replace

import numpy as np


@dataclass(frozen=True)
class SummaryConfig:
    """Grids for the negative-binomial group summaries (log2 units)."""

    mu_max: float = 6.0
    n_mu: int = 201
    sig2_min: float = 1e-3
    sig2_max: float = 10.0
    n_sig2: int = 81
    theta_step: float = 0.1
    chunk: int = 8

    def grids(self):
        mu = np.linspace(-self.mu_max, self.mu_max, self.n_mu)
        sig2 = np.geomspace(self.sig2_min, self.sig2_max, self.n_sig2)
        half = self.mu_max + 5 * np.sqrt(self.sig2_max)
        theta = np.linspace(-half, half, int(np.ceil(2 * half / self.theta_step)) + 1)
        return mu, sig2, theta


@dataclass(frozen=True)
class ModelConfig:
    """Global constants of the two-layer segmentation model (log2 units)."""

    tau2: float = 0.4               # nugget: base-level shape shared by all groups
    b_var: float = 1.0              # prior variance of the consensus level
    omega_scale: float = 0.2        # half-Cauchy scale of group deviations
    cons_mean: float = 12.0         # consensus block lengths: discretized gamma (resolution choice)
    cons_shape: float = 3.0
    cons_min: int = 2
    pattern_mean: float = 20.0      # group-pattern block lengths (resolution choice)
    pattern_shape: float = 3.0
    pattern_min: int = 2
    l_max: int = 150
    n_draws: int = 1000
    seed: int = 0
    prob: float = 0.95              # call cutoff
    delta: float = 0.5              # |Delta| cutoff for differential calls

    def with_(self, **changes) -> "ModelConfig":
        return replace(self, **changes)

    def to_json(self, path=None):
        text = json.dumps(asdict(self))
        if path is None:
            return text
        with open(path, "w") as handle:
            handle.write(text)

    @classmethod
    def from_json(cls, text_or_path: str) -> "ModelConfig":
        text = text_or_path if text_or_path.lstrip().startswith("{") else open(text_or_path).read()
        return cls(**json.loads(text))


DEFAULT_SUMMARY = SummaryConfig()
DEFAULT_MODEL = ModelConfig()
