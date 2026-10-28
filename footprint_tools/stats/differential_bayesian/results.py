"""Per-base tracks and calls of one DHS, computed from nested-model posterior draws."""

from dataclasses import dataclass

import numpy as np
import pandas as pd

from .config import DEFAULT_NESTED, NestedConfig
from .nested import NestedModel


@dataclass(frozen=True)
class NestedFit:
    group_names: tuple[str, ...]
    depth_mean: np.ndarray     # n x G, log2
    depth_lo: np.ndarray       # n x G, 5% quantile
    depth_hi: np.ndarray       # n x G, 95% quantile
    p_any: np.ndarray          # n, P(min_g depth < 0)
    delta_mean: np.ndarray     # n x G, E[depth_g - median of the other groups]
    p_more: np.ndarray         # n x G, P(Delta < -delta)
    p_less: np.ndarray         # n x G, P(Delta > +delta)
    p_outer: np.ndarray        # n - 1, P(outer boundary between j and j + 1)
    p_inner: np.ndarray        # n - 1, P(inner-only boundary between j and j + 1)
    icc: np.ndarray            # n, E[omega^2 / (omega^2 + sigma2_bar)]
    log_evidence: float
    config: str                # NestedConfig as JSON

    @property
    def footprints(self):
        return runs(self.p_any >= NestedConfig.from_json(self.config).prob)

    @property
    def calls(self):
        prob = NestedConfig.from_json(self.config).prob
        return [(name, kind, s, e)
                for g, name in enumerate(self.group_names)
                for kind, p in (("more", self.p_more), ("less", self.p_less))
                for s, e in runs(p[:, g] >= prob)]

    def region_table(self) -> pd.DataFrame:
        """One row per run of footprinted / non-footprinted bases, with differential base counts."""
        prob = NestedConfig.from_json(self.config).prob
        fp, more, less = self.p_any >= prob, self.p_more >= prob, self.p_less >= prob
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


def fit_nested(x, v, group_names, config: NestedConfig = DEFAULT_NESTED) -> NestedFit:
    """Fit one DHS from group summaries x, v (groups x positions, log2)."""
    model = NestedModel(x, v, config)
    depth, p_outer, p_inner, icc = model.sample(config.n_draws, config.seed)
    delta = delta_from_others(depth)
    lo, hi = np.quantile(depth, [0.05, 0.95], axis=0)
    f32 = lambda a: np.asarray(a, np.float32)
    return NestedFit(
        tuple(map(str, group_names)), f32(depth.mean(0)), f32(lo), f32(hi),
        f32((depth.min(2) < 0).mean(0)), f32(delta.mean(0)),
        f32((delta < -config.delta).mean(0)), f32((delta > config.delta).mean(0)),
        f32(p_outer), f32(p_inner), f32(icc), model.log_evidence, config.to_json(),
    )


def delta_from_others(depth):
    """Delta_g = depth_g - median depth of the other groups; depth is draws x n x G."""
    G = depth.shape[2]
    return depth - np.stack([np.median(np.delete(depth, g, axis=2), axis=2) for g in range(G)], axis=2)


def runs(mask):
    """Half-open (start, end) runs of True in a 1-D boolean array."""
    edges = np.flatnonzero(np.diff(np.r_[0, np.asarray(mask, int), 0]))
    return list(zip(edges[::2].tolist(), edges[1::2].tolist()))
