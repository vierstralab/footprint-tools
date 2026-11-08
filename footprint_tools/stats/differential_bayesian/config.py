import json
from dataclasses import asdict, dataclass

import numpy as np


@dataclass(frozen=True)
class DifferentialFootprintConfig:
    """Stage-1 segmentation prior, stage-2 test and call rule (depths in log2)."""

    theta_min: float = -9.0         # depth grid
    theta_max: float = 4.0
    theta_step: float = 0.2
    level_sd: float = 3.0           # block level ~ N(0, level_sd^2), wide on purpose
    length_mean: float = 24.8       # block lengths ~ discretized gamma
    length_shape: float = 1.9
    length_min: int = 2
    length_max: int = 150
    chunk: int = 40                 # samples segmented at once
    d_min: float = 0.6              # call: effect d of the set on the called side at least this (d = 2W - 1)
    report_d: tuple = (0.0, 0.2, 0.4, 0.5, 0.6)     # thresholds with p-values kept for the genome-wide report
    k_max: int = 5                  # sets of up to k_max groups
    n_perm: int = 1000              # label permutations
    perm_batch: int = 40            # permutations per matrix product
    tail_frac: float = 0.1          # generalised Pareto tail fitted to the top tail_frac of the permutation maxima
    min_hits: int = 10              # fewer permutations reaching T: use the tail
    fdr: float = 0.05               # Benjamini-Hochberg level, over all positions genome-wide
    min_len: int = 5                # shortest call, group span inside a call, footprint, and DHS window (bp)
    fold: float = float(np.log2(1.5))   # footprint: some group's mean depth < -fold
    p_footprint: float = 0.95       # footprint: P(min_g mean depth < -fold)
    seed: int = 0

    def grid(self) -> np.ndarray:
        return np.arange(self.theta_min, self.theta_max + self.theta_step / 2, self.theta_step)

    def to_json(self) -> str:
        return json.dumps(asdict(self))

    @classmethod
    def from_json(cls, text: str) -> "DifferentialFootprintConfig":
        d = json.loads(text)
        if 'w_min' in d:                                    # saved before W was replaced by d = 2W - 1
            d['d_min'] = round(2 * d.pop('w_min') - 1, 6)
            d['report_d'] = [round(2 * w - 1, 6) for w in d.pop('report_w')]
        return cls(**{**d, 'report_d': tuple(d['report_d'])})


DEFAULT_DF_CONFIG = DifferentialFootprintConfig()
