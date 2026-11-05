"""One DHS: stage-1 posteriors per sample, stage-2 tracks per group, and the calls derived from them."""

from dataclasses import dataclass

import numpy as np
import pandas as pd

from .calls import discriminability, group_tracks, p_footprint, runs
from .config import DEFAULT_DF_CONFIG, DifferentialFootprintConfig
from .segmentation import segment_samples


@dataclass(frozen=True)
class DifferentialFootprintFit:
    group_names: tuple[str, ...]
    sample_names: tuple[str, ...]
    sample_groups: tuple[str, ...]
    theta_mean: np.ndarray     # S x n, stage-1 posterior mean depth (log2)
    theta_var: np.ndarray      # S x n, stage-1 posterior variance
    mean: np.ndarray           # G x n, group mean depth
    se: np.ndarray             # G x n, its standard error
    delta: np.ndarray          # G x n, mean - median of the other groups' means
    a_more: np.ndarray         # G x n, discriminability, more depleted than the other groups
    a_less: np.ndarray         # G x n, less depleted
    lo_more: np.ndarray        # G x n, 95% lower bounds
    lo_less: np.ndarray
    p_any: np.ndarray          # n, P(some group's mean depth < -fold)
    config: str                # DifferentialFootprintConfig as JSON

    @property
    def cfg(self) -> DifferentialFootprintConfig:
        return DifferentialFootprintConfig.from_json(self.config)

    def called(self, kind: str) -> np.ndarray:
        """G x n mask of "more" or "less" calls: lower bound of A >= a_min and |Delta| > fold, runs >= min_len."""
        cfg, sign = self.cfg, -1 if kind == "more" else 1
        ok = (getattr(self, f"lo_{kind}") >= cfg.a_min) & (sign * self.delta > cfg.fold)
        out = np.zeros(ok.shape, bool)
        for g, row in enumerate(ok):
            for s, e in runs(row, cfg.min_len):
                out[g, s:e] = True
        return out

    @property
    def calls(self) -> pd.DataFrame:
        """One row per call, ranked by the mean lower bound of A."""
        rows = [dict(group=name, kind=kind, start=s, end=e, effect=abs(self.delta[g, s:e].mean()),
                     a=getattr(self, f"a_{kind}")[g, s:e].mean(), a_lo=getattr(self, f"lo_{kind}")[g, s:e].mean())
                for kind in ("more", "less") for g, name in enumerate(self.group_names) for s, e in runs(self.called(kind)[g])]
        return pd.DataFrame(rows, columns=["group", "kind", "start", "end", "effect", "a", "a_lo"]).sort_values("a_lo", ascending=False)

    @property
    def footprinted(self) -> np.ndarray:
        cfg = self.cfg; out = np.zeros(self.p_any.shape, bool)
        for s, e in runs(self.p_any >= cfg.p_footprint, cfg.min_len):
            out[s:e] = True
        return out

    def region_table(self) -> pd.DataFrame:
        """One row per run of footprinted / non-footprinted bases, with differential base counts."""
        fp, more, less = self.footprinted, self.called("more"), self.called("less")
        cuts = np.r_[0, np.flatnonzero(np.diff(fp)) + 1, fp.size]
        rows = []
        for i, (s, e) in enumerate(zip(cuts[:-1], cuts[1:])):
            m, l = more[:, s:e], less[:, s:e]
            row = dict(region_id=i, start=s, end=e, width=e - s, status="footprinted" if fp[s] else "not_footprinted",
                       p_any_mean=self.p_any[s:e].mean(), diff_more_n_bases=m.any(0).sum(), diff_less_n_bases=l.any(0).sum(),
                       diff_n_bases=(m | l).any(0).sum())
            for g, name in enumerate(self.group_names):
                row[f"{name}__more_n_bases"], row[f"{name}__less_n_bases"] = m[g].sum(), l[g].sum()
            rows.append(row)
        return pd.DataFrame(rows)


def fit_differential_footprints(groups: pd.Series, obs: pd.DataFrame, exp: pd.DataFrame, disp_models: pd.Series,
                                config: DifferentialFootprintConfig = DEFAULT_DF_CONFIG) -> DifferentialFootprintFit:
    """groups: sample -> group label; obs, exp: samples x positions; disp_models: sample -> dispersion model."""
    samples, labels = groups.index, groups.to_numpy().astype(str)
    names = tuple(pd.unique(labels)); idx = [np.flatnonzero(labels == g) for g in names]
    if min(map(len, idx)) < 2:
        raise ValueError("every group needs at least two samples")
    m, v = segment_samples(obs.loc[samples].to_numpy(float), exp.loc[samples].to_numpy(float), disp_models.loc[samples].to_numpy(), config)
    mean, se, delta = group_tracks(m, v, idx)
    disc = discriminability(m, v, idx, config.n_boot, config.seed)
    f32 = lambda a: np.asarray(a, np.float32)
    return DifferentialFootprintFit(names, tuple(map(str, samples)), tuple(labels), f32(m), f32(v), f32(mean), f32(se), f32(delta),
                                    *(f32(disc[k]) for k in ("a_more", "a_less", "lo_more", "lo_less")),
                                    f32(p_footprint(mean, se, config.fold)), config.to_json())
