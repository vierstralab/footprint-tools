"""Per-base tracks from posterior draws, runs, and the per-DHS region table."""

import numpy as np
import pandas as pd


def runs(mask: np.ndarray) -> list[tuple[int, int]]:
    """Half-open (start, end) runs of True in a 1-D boolean array."""
    edges = np.flatnonzero(np.diff(np.r_[0, np.asarray(mask, int), 0]))
    return list(zip(edges[::2].tolist(), edges[1::2].tolist()))


def delta_from_others(depth: np.ndarray) -> np.ndarray:
    """Delta_g = depth_g - median of the other groups; depth is draws x n x G."""
    G = depth.shape[2]
    if G == 1:
        return np.zeros_like(depth)
    ref = np.stack([np.median(np.delete(depth, g, axis=2), axis=2) for g in range(G)], axis=2)
    return depth - ref


def tracks(depth: np.ndarray, delta: float) -> dict[str, np.ndarray]:
    dl = delta_from_others(depth)
    lo, hi = np.quantile(depth, [0.05, 0.95], axis=0)
    return dict(
        depth_mean=depth.mean(axis=0),
        depth_lo=lo,
        depth_hi=hi,
        p_any=(depth.min(axis=2) < 0).mean(axis=0),
        delta_mean=dl.mean(axis=0),
        p_more=(dl < -delta).mean(axis=0),
        p_less=(dl > delta).mean(axis=0),
    )


def region_table(fit, prob: float | None = None) -> pd.DataFrame:
    """One row per run of footprinted / non-footprinted bases, with differential counts."""
    prob = fit.config.prob if prob is None else prob
    more, less = (fit.p_more >= prob).T, (fit.p_less >= prob).T
    footprinted = fit.p_any >= prob
    bounds = np.flatnonzero(np.diff(footprinted.astype(int))) + 1
    starts, ends = np.r_[0, bounds], np.r_[bounds, footprinted.size]
    rows = []
    for i, (s, e) in enumerate(zip(starts.tolist(), ends.tolist())):
        sl = slice(s, e)
        any_more, any_less = more[:, sl].any(axis=0), less[:, sl].any(axis=0)
        row = dict(
            region_id=i, start=s, end=e, width=e - s,
            status="footprinted" if footprinted[s] else "not_footprinted",
            p_any_mean=float(fit.p_any[sl].mean()),
            icc_mean=float(np.nanmean(fit.icc[sl])) if np.any(np.isfinite(fit.icc[sl])) else np.nan,
            diff_more_n_bases=int(any_more.sum()),
            diff_less_n_bases=int(any_less.sum()),
            diff_n_bases=int((any_more | any_less).sum()),
            diff_more_n_regions=len(runs(any_more)),
            diff_less_n_regions=len(runs(any_less)),
            diff_n_regions=len(runs(any_more | any_less)),
        )
        for g, name in enumerate(fit.group_names):
            row[f"{name}__more_n_bases"] = int(more[g, sl].sum())
            row[f"{name}__less_n_bases"] = int(less[g, sl].sum())
        rows.append(row)
    return pd.DataFrame(rows)
