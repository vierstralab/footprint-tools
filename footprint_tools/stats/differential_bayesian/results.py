"""NestedFit: the stored per-DHS result of the nested model, and the fit entry point."""

import json
from dataclasses import dataclass
from typing import ClassVar

import numpy as np

from .calls import region_table, runs, tracks
from .config import DEFAULT_NESTED, NestedConfig
from .io import Serializable, tuple_str
from .nested import NestedModel


@dataclass(frozen=True, slots=True)
class NestedFit(Serializable):
    save_attrs: ClassVar[tuple[str, ...]] = (
        "group_names", "depth_mean", "depth_lo", "depth_hi", "p_any", "delta_mean",
        "p_more", "p_less", "p_outer", "p_inner", "icc", "log_evidence", "config_json",
    )
    optional_save_attrs: ClassVar[tuple[str, ...]] = ("draws",)

    group_names: tuple[str, ...]
    depth_mean: np.ndarray     # n x G, log2
    depth_lo: np.ndarray       # n x G, 5% quantile
    depth_hi: np.ndarray       # n x G, 95% quantile
    p_any: np.ndarray          # n, P(min_g depth < 0)
    delta_mean: np.ndarray     # n x G, E[depth_g - median of other groups]
    p_more: np.ndarray         # n x G, P(Delta < -delta)
    p_less: np.ndarray         # n x G, P(Delta > +delta)
    p_outer: np.ndarray        # n - 1, P(outer boundary between j and j + 1)
    p_inner: np.ndarray        # n - 1, P(inner-only boundary between j and j + 1)
    icc: np.ndarray            # n, E[omega^2 / (omega^2 + sigma2_bar)]
    log_evidence: float
    config_json: str
    draws: np.ndarray | None = None   # draws x n x G, float16

    @property
    def config(self) -> NestedConfig:
        return NestedConfig.from_dict(json.loads(self.config_json))

    @property
    def footprints(self) -> list[tuple[int, int]]:
        return runs(self.p_any >= self.config.prob)

    @property
    def calls(self) -> list[tuple[str, str, int, int]]:
        prob = self.config.prob
        return [
            (name, kind, s, e)
            for g, name in enumerate(self.group_names)
            for kind, p in (("more", self.p_more), ("less", self.p_less))
            for s, e in runs(p[:, g] >= prob)
        ]

    def region_table(self, prob: float | None = None):
        return region_table(self, prob)

    def to_dict(self) -> dict[str, object]:
        out = Serializable.to_dict(self)
        out["group_names"] = np.asarray(self.group_names, dtype=str)
        out["log_evidence"] = np.asarray(self.log_evidence)
        out["config_json"] = np.asarray(self.config_json)
        return out

    @classmethod
    def from_dict(cls, data: dict[str, object]) -> "NestedFit":
        values = {name: data[name] for name in cls.save_attrs}
        values["group_names"] = tuple_str(data["group_names"])
        values["log_evidence"] = float(np.asarray(data["log_evidence"]))
        values["config_json"] = str(np.asarray(data["config_json"]))
        return cls(**values, draws=data.get("draws"))


def fit_nested(x, v, group_names, config: NestedConfig = DEFAULT_NESTED) -> NestedFit:
    """Fit one DHS from group summaries x, v (groups x positions, log2 units)."""
    model = NestedModel(x, v, config)
    draws = model.sample(config.n_draws, config.seed)
    tr = {k: v.astype(np.float32) for k, v in tracks(draws.depth, config.delta).items()}
    return NestedFit(
        group_names=tuple(map(str, group_names)),
        p_outer=draws.p_outer.astype(np.float32),
        p_inner=draws.p_inner.astype(np.float32),
        icc=draws.icc.astype(np.float32),
        log_evidence=model.log_evidence(),
        config_json=json.dumps(config.to_dict()),
        draws=draws.depth.astype(np.float16) if config.store_draws else None,
        **tr,
    )
