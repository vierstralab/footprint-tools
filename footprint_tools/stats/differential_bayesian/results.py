"""Per-base tracks and calls of one DHS, computed from posterior draws."""

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
import pandas as pd

from .config import DEFAULT_DF_MODEL, DifferentialFootprintModelConfig
from .model import DifferentialFootprintModel


@dataclass(frozen=True)
class DifferentialFootprintFit:
    group_names: tuple[str, ...]
    depth_mean: np.ndarray     # n x G, log2
    depth_lo: np.ndarray       # n x G, 5% quantile
    depth_hi: np.ndarray       # n x G, 95% quantile
    p_any: np.ndarray          # n, P(min_g depth < 0)
    delta_mean: np.ndarray     # n x G, E[Delta_g], Delta_g = u_g - median of the other groups' u, u = observable(depth)
    p_more: np.ndarray         # n x G, P(Delta < -delta)
    p_less: np.ndarray         # n x G, P(Delta > +delta)
    p_pattern: np.ndarray      # n - 1, P(group-pattern boundary between j and j + 1)
    p_consensus: np.ndarray    # n - 1, P(consensus boundary between j and j + 1)
    icc: np.ndarray            # n, E[omega^2 / (omega^2 + within-group variance)]
    log_evidence: float
    config: str                # ModelConfig as JSON

    @property
    def footprints(self) -> list[tuple[int, int]]:
        return runs(self.p_any >= DifferentialFootprintModelConfig.from_json(self.config).prob)

    def called(self, kind: str) -> np.ndarray:
        """n x G mask of "more" or "less" calls: runs with probability >= prob and at least min_len bases."""
        cfg = DifferentialFootprintModelConfig.from_json(self.config)
        p = self.p_more if kind == "more" else self.p_less
        return np.stack([long_runs(p[:, g] >= cfg.prob, cfg.min_len) for g in range(p.shape[1])], axis=1)

    @property
    def calls(self) -> list[tuple[str, str, int, int]]:
        return [(name, kind, s, e) for kind in ("more", "less")
                for name, col in zip(self.group_names, self.called(kind).T) for s, e in runs(col)]

    def region_table(self) -> pd.DataFrame:
        """One row per run of footprinted / non-footprinted bases, with differential base counts."""
        fp, more, less = self.p_any >= DifferentialFootprintModelConfig.from_json(self.config).prob, self.called("more"), self.called("less")
        cuts = np.r_[0, np.flatnonzero(np.diff(fp)) + 1, fp.size]
        rows = []
        for i, (s, e) in enumerate(zip(cuts[:-1], cuts[1:])):
            m, l = more[s:e], less[s:e]
            row = dict(region_id=i, start=s, end=e, width=e - s,
                       status="footprinted" if fp[s] else "not_footprinted",
                       p_any_mean=self.p_any[s:e].mean(), icc_mean=self.icc[s:e].mean(),
                       diff_more_n_bases=m.any(1).sum(), diff_less_n_bases=l.any(1).sum(),
                       diff_n_bases=(m | l).any(1).sum())
            for g, name in enumerate(self.group_names):
                row[f"{name}__more_n_bases"], row[f"{name}__less_n_bases"] = m[:, g].sum(), l[:, g].sum()
            rows.append(row)
        return pd.DataFrame(rows)


def fit_differential_footprints(x: np.ndarray, v: np.ndarray, sigma2: np.ndarray, group_names: Sequence[str],
                   config: DifferentialFootprintModelConfig = DEFAULT_DF_MODEL) -> DifferentialFootprintFit:
    """Fit one DHS from group summaries: x, v, sigma2 (groups x positions, log2); v is widened by the group nugget."""
    model = DifferentialFootprintModel(x, np.asarray(v) + config.nugget(group_names)[:, None], sigma2, config)
    depth, p_pattern, p_consensus, icc = model.sample(config.n_draws, config.seed)
    delta = delta_from_others(observable(depth, config.call_floor))
    lo, hi = np.quantile(depth, [0.05, 0.95], axis=0)
    f32 = lambda a: np.asarray(a, np.float32)
    return DifferentialFootprintFit(
        tuple(map(str, group_names)), f32(depth.mean(0)), f32(lo), f32(hi),
        f32((depth.min(2) < 0).mean(0)), f32(delta.mean(0)),
        f32((delta < -config.delta).mean(0)), f32((delta > config.delta).mean(0)),
        f32(p_pattern), f32(p_consensus), f32(icc),
        float(sum(model.log_evidence()) + model.const), config.to_json(),
    )


def observable(depth: np.ndarray, c: float) -> np.ndarray:
    """Depth as a sample with floor c sees it: groups are compared where the counts can tell them apart."""
    return np.log2(c + (1 - c) * np.exp2(depth))


def delta_from_others(depth: np.ndarray) -> np.ndarray:
    """Delta_g = depth_g - median depth of the other groups; depth is draws x n x G."""
    G = depth.shape[2]
    return depth - np.stack([np.median(np.delete(depth, g, axis=2), axis=2) for g in range(G)], axis=2)


def long_runs(mask: np.ndarray, min_len: int) -> np.ndarray:
    out = np.zeros(len(mask), bool)
    for s, e in runs(mask):
        out[s:e] = e - s >= min_len
    return out


def runs(mask: np.ndarray) -> list[tuple[int, int]]:
    """Half-open (start, end) runs of True in a 1-D boolean array."""
    edges = np.flatnonzero(np.diff(np.r_[0, np.asarray(mask, int), 0]))
    return list(zip(edges[::2].tolist(), edges[1::2].tolist()))
