"""One DHS: stage-1 posteriors per sample, footprints, and the per-position differential p-values; calls once the
genome-wide cutoff is known."""

from dataclasses import dataclass

import numpy as np
import pandas as pd

from .calls import differential_calls, p_footprint, runs
from .config import DEFAULT_DF_CONFIG, DifferentialFootprintConfig
from .pvalues import differential_pvalues, set_weights
from .segmentation import segment_samples


@dataclass(frozen=True)
class DifferentialFootprintFit:
    group_names: tuple[str, ...]
    sample_names: tuple[str, ...]
    sample_groups: tuple[str, ...]
    theta_mean: np.ndarray     # S x n, stage-1 posterior mean depth (log2)
    theta_var: np.ndarray      # S x n, stage-1 posterior variance
    mean: np.ndarray           # G x n, group mean depth
    p_any: np.ndarray          # n, P(some group's mean depth < -fold)
    lp_more: np.ndarray        # n, log10 p of the best more-depleted set at W >= w_min
    set_more: np.ndarray       # n, its index in set_weights(G, k_max)
    w_more: np.ndarray         # n, its W on the called side
    lp_less: np.ndarray
    set_less: np.ndarray
    w_less: np.ndarray
    lp: np.ndarray             # len(report_w) x n, log10 p of the best hypothesis per position, per threshold
    lp_dhs: np.ndarray         # len(report_w), log10 p per DHS (best window of min_len positions)
    xi: np.ndarray             # n, generalised Pareto shape fitted to the permutation maxima
    c95: np.ndarray            # n, 95% point of the permutation maximum
    config: str                # DifferentialFootprintConfig as JSON

    @property
    def cfg(self) -> DifferentialFootprintConfig:
        return DifferentialFootprintConfig.from_json(self.config)

    @property
    def footprinted(self) -> np.ndarray:
        out = np.zeros(self.p_any.shape, bool)
        for s, e in runs(self.p_any >= self.cfg.p_footprint, self.cfg.min_len):
            out[s:e] = True
        return out

    def calls(self, cutoff: float) -> pd.DataFrame:
        """Differential calls at the genome-wide log10 p cutoff, one row per call."""
        cfg, names = self.cfg, self.group_names
        set_list = set_weights(len(names), cfg.k_max)[0]
        C = differential_calls(dict(more=self.lp_more, less=self.lp_less), dict(more=self.set_more, less=self.set_less),
                               dict(more=self.w_more, less=self.w_less), set_list, cutoff, cfg.min_len)
        fp = self.footprinted
        rows = [dict(kind=c['kind'], start=c['start'], end=c['end'],
                     groups=' + '.join(names[k] for k, sp in c['spans'].items() if sp == [(c['start'], c['end'])]),
                     partial='; '.join(f"{names[k]} {','.join(f'{a}-{b}' for a, b in sp)}" for k, sp in c['spans'].items() if sp != [(c['start'], c['end'])]),
                     w=round(c['w'], 4), min_log10p=round(c['min_log10p'], 3), in_footprint=round(float(fp[c['start']:c['end']].mean()), 3))
                for c in C]
        return pd.DataFrame(rows, columns=['kind', 'start', 'end', 'groups', 'partial', 'w', 'min_log10p', 'in_footprint'])


def fit_differential_footprints(groups: pd.Series, obs: pd.DataFrame, exp: pd.DataFrame, disp_models: pd.Series,
                                config: DifferentialFootprintConfig = DEFAULT_DF_CONFIG) -> DifferentialFootprintFit:
    """groups: sample -> group label; obs, exp: samples x positions; disp_models: sample -> dispersion model."""
    samples, labels = groups.index, groups.to_numpy().astype(str)
    names = tuple(pd.unique(labels)); g = np.array([names.index(x) for x in labels])
    if min(np.bincount(g)) < 2:
        raise ValueError("every group needs at least two samples")
    m, v = segment_samples(obs.loc[samples].to_numpy(float), exp.loc[samples].to_numpy(float), disp_models.loc[samples].to_numpy(), config)
    mean, p_any = p_footprint(m, v, [np.flatnonzero(g == k) for k in range(len(names))], config.fold)
    p = differential_pvalues(m, v, g, len(names), config)
    f32 = lambda a: np.asarray(a, np.float32)
    return DifferentialFootprintFit(
        names, tuple(map(str, samples)), tuple(labels), f32(m), f32(v), f32(mean), f32(p_any),
        f32(p['lp_more']), p['set_more'].astype(np.int16), f32(p['w_more']), f32(p['lp_less']), p['set_less'].astype(np.int16), f32(p['w_less']),
        f32(p['lp']), f32(p['lp_dhs']), f32(p['xi']), f32(p['c95']), config.to_json())
