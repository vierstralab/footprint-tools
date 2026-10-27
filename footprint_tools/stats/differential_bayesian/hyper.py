"""Offline fit of the global hyperparameters of the nested model.

1. kappa: split-half calibration on counts (each group's samples split in two;
   kappa chosen so that disagreement between halves / posterior SD has SD 1).
2. tau2, omega_scale (optionally outer_mean): grid search on the log evidence summed
   over DHSs. Block-length means are resolution choices: the evidence favours blocks of
   a few bp at both levels, so they are fitted only when asked.

python -m footprint_tools.stats.differential_bayesian.hyper \
    --summaries results/*.npz [--counts counts/*.npz] --out hyper.json
"""

import argparse
import itertools
import json

import numpy as np
import pandas as pd

from .calls import delta_from_others
from .config import DEFAULT_NESTED, NestedConfig
from .nested import NestedModel
from .priors import gamma_length_log_mass, outer_length_log_weight


def summed_log_evidence(summaries, base: NestedConfig, grid: dict) -> pd.DataFrame:
    """Summed log evidence over DHSs for every grid point.

    `summaries` is a list of (x, v); `grid` maps NestedConfig field names to values.
    Models are rebuilt only when a field other than the length means changes.
    """
    names = list(grid)
    length_keys = [k for k in names if k in ("inner_mean", "outer_mean")]
    model_keys = [k for k in names if k not in length_keys]
    rows = []
    for model_values in itertools.product(*(grid[k] for k in model_keys)):
        cfg = base.with_(**dict(zip(model_keys, model_values)))
        models = [NestedModel(x, v, cfg) for x, v in summaries]
        for length_values in itertools.product(*(grid[k] for k in length_keys)):
            lcfg = cfg.with_(**dict(zip(length_keys, length_values)))
            total = 0.0
            for M in models:
                M.lq = gamma_length_log_mass(M.l_max, lcfg.inner_mean, lcfg.len_shape, lcfg.len_min)
                lq_out = gamma_length_log_mass(M.l_max, lcfg.outer_mean, lcfg.len_shape, lcfg.len_min)
                M.lo = outer_length_log_weight(M.lq, lq_out)
                M._part = None
                total += M.log_evidence()
            rows.append({**dict(zip(model_keys, model_values)), **dict(zip(length_keys, length_values)), "log_evidence": total})
    return pd.DataFrame(rows).sort_values("log_evidence", ascending=False, ignore_index=True)


def log_ratio_summaries(obs, exp, sample_groups, group_names):
    """Plain log2-ratio group means and variances, used for split-half calibration."""
    r = np.log2((obs + 0.5) / (exp + 0.5))
    x = np.stack([r[sample_groups == g].mean(0) for g in group_names])
    v = np.stack([r[sample_groups == g].var(0, ddof=1) / (sample_groups == g).sum() for g in group_names])
    return x, np.maximum(v, 0.25 * np.median(v))


def split_half_sd(obs, exp, sample_groups, config: NestedConfig, kappas, n_rep=2, n_draws=500):
    """SD of (Delta_A - Delta_B) / sqrt(sd_A^2 + sd_B^2) between two sample halves, per kappa."""
    sample_groups = np.asarray(sample_groups).astype(str)
    names = list(dict.fromkeys(sample_groups))
    small = [g for g in names if (sample_groups == g).sum() < 4]
    if small:
        raise ValueError(f"split-half calibration needs at least 4 samples per group: {small}")
    out = {k: [] for k in kappas}
    for rep in range(n_rep):
        rng = np.random.default_rng(rep)
        halves = [[], []]
        for g in names:
            idx = rng.permutation(np.flatnonzero(sample_groups == g))
            h = idx.size // 2
            halves[0] += idx[:h].tolist()
            halves[1] += idx[h:2 * h].tolist()
        for kappa in kappas:
            cfg = config.with_(kappa=kappa)
            stats = []
            for half in halves:
                half = np.asarray(half)
                x, v = log_ratio_summaries(obs[half], exp[half], sample_groups[half], names)
                dl = delta_from_others(NestedModel(x, v, cfg).sample(n_draws, rep).depth)
                stats.append((dl.mean(0), dl.std(0)))
            (ma, sa), (mb, sb) = stats
            out[kappa].append(float(np.std((ma - mb) / np.sqrt(sa ** 2 + sb ** 2))))
    return {k: float(np.median(v)) for k, v in out.items()}


def calibrate_kappa(counts, config: NestedConfig, kappas=(1.0, 1.5, 2.0, 2.5, 3.0, 3.5)):
    """kappa whose median split-half SD over DHSs is closest to 1 (log-linear interpolation)."""
    sds = np.array([[split_half_sd(o, e, g, config, kappas)[k] for k in kappas] for o, e, g in counts])
    med = np.median(sds, axis=0)
    lk = np.log(kappas)
    order = np.argsort(med)
    return float(np.exp(np.interp(0.0, np.log(med[order]), lk[order]))), dict(zip(kappas, med.tolist()))


def _load_summaries(path):
    with np.load(path, allow_pickle=False) as h:
        return h["differential.mu_mean"].astype(float), h["differential.mu_var"].astype(float)


def _load_counts(path):
    with np.load(path, allow_pickle=False) as h:
        return h["obs"].astype(float), h["exp"].astype(float), h["sample_groups"].astype(str)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--summaries", nargs="+", required=True, help="NPZ files written by save_data_results")
    p.add_argument("--counts", nargs="*", default=(), help="NPZ files with obs, exp, sample_groups (for kappa)")
    p.add_argument("--base", help="JSON NestedConfig to start from")
    p.add_argument("--tau2", type=float, nargs="+", default=[0.25, 0.3, 0.35, 0.4, 0.5])
    p.add_argument("--omega-scale", type=float, nargs="+", default=[0.05, 0.1, 0.15, 0.2, 0.3])
    p.add_argument("--outer-mean", type=float, nargs="*", default=[], help="fit outer_mean on this grid")
    p.add_argument("--out", required=True)
    a = p.parse_args(argv)

    cfg = NestedConfig.from_json(a.base) if a.base else DEFAULT_NESTED
    report = {}
    if a.counts:
        kappa, curve = calibrate_kappa([_load_counts(f) for f in a.counts], cfg)
        cfg = cfg.with_(kappa=round(kappa, 2))
        report["kappa_curve"] = curve
    grid = {"tau2": a.tau2, "omega_scale": a.omega_scale}
    if a.outer_mean:
        grid["outer_mean"] = a.outer_mean
    table = summed_log_evidence([_load_summaries(f) for f in a.summaries], cfg, grid)
    best = table.iloc[0]
    cfg = cfg.with_(**{k: float(best[k]) for k in grid})
    cfg.to_json(a.out)
    table.to_csv(a.out.rsplit(".", 1)[0] + ".evidence.tsv", sep="\t", index=False)
    report["best"] = cfg.to_dict()
    report["top"] = table.head(10).to_dict("records")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
