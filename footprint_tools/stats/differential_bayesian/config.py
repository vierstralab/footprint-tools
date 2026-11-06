import json
from dataclasses import asdict, dataclass

import numpy as np


@dataclass(frozen=True)
class DifferentialFootprintConfig:
    """Stage-1 segmentation prior and stage-2 call rule (log2 units)."""

    theta_min: float = -9.0         # depth grid
    theta_max: float = 4.0
    theta_step: float = 0.2
    level_sd: float = 3.0           # block level ~ N(0, level_sd^2), wide on purpose
    length_mean: float = 24.8       # block lengths ~ discretized gamma
    length_shape: float = 1.9
    length_min: int = 2
    length_max: int = 150
    chunk: int = 40                 # samples segmented at once
    a_min: float = 0.75             # call: 95% lower bound of the discriminability A
    fold: float = float(np.log2(1.5))   # call: |Delta| > fold; footprint: some group mean depth < -fold
    p_footprint: float = 0.95       # footprint: P(min_g mean depth < -fold)
    min_len: int = 5                # shortest call or footprint, bp
    n_boot: int = 2000
    seed: int = 0

    def grid(self) -> np.ndarray:
        return np.arange(self.theta_min, self.theta_max + self.theta_step / 2, self.theta_step)

    def to_json(self) -> str:
        return json.dumps(asdict(self))

    @classmethod
    def from_json(cls, text: str) -> "DifferentialFootprintConfig":
        return cls(**json.loads(text))


DEFAULT_DF_CONFIG = DifferentialFootprintConfig()
