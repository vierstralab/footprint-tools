"""Offline fit of the global constants, once per data release.

kappa: split-half calibration on counts. Each group's samples are split in two halves;
kappa is chosen so that (Delta_A - Delta_B) / sqrt(sd_A^2 + sd_B^2) has SD 1.
tau2, omega_scale: maximum of the log evidence summed over DHSs. Block-length means are
resolution choices (the evidence favours blocks of a few bp) and are not fitted.

python -m footprint_tools.stats.differential_bayesian.hyper \
    --summaries results/*.npz --counts counts/*.npz --out nested.json
"""

import argparse
import itertools

import numpy as np

from .config import DEFAULT_NESTED, NestedConfig
from .nested import NestedModel
from .results import delta_from_others


def summed_log_evidence(summaries, config: NestedConfig, grid: dict) -> dict:
    """{grid point: log evidence summed over (x, v) pairs}; grid maps NestedConfig fields to values."""
    return {point: sum(NestedModel(x, v, config.with_(**dict(zip(grid, point)))).log_evidence for x, v in summaries)
            for point in itertools.product(*grid.values())}


def split_half_sd(obs, exp, sample_groups, config: NestedConfig, n_draws=500, seed=0):
    """SD of the standardized disagreement in Delta between two random halves of each group."""
    rng = np.random.default_rng(seed)
    halves = [[], []]
    for g in np.unique(sample_groups):
        idx = rng.permutation(np.flatnonzero(sample_groups == g))
        if idx.size < 4:
            raise ValueError(f"group {g} needs at least 4 samples")
        halves[0] += idx[: idx.size // 2].tolist()
        halves[1] += idx[idx.size // 2: 2 * (idx.size // 2)].tolist()
    r = np.log2((obs + 0.5) / (exp + 0.5))
    stats = []
    for half in map(np.array, halves):
        grp = sample_groups[half]
        x = np.stack([r[half][grp == g].mean(0) for g in np.unique(grp)])
        v = np.stack([r[half][grp == g].var(0, ddof=1) / (grp == g).sum() for g in np.unique(grp)])
        dl = delta_from_others(NestedModel(x, np.maximum(v, 0.25 * np.median(v)), config).sample(n_draws, seed)[0])
        stats.append((dl.mean(0), dl.std(0)))
    (ma, sa), (mb, sb) = stats
    return float(np.std((ma - mb) / np.hypot(sa, sb)))


def calibrate_kappa(counts, config: NestedConfig, kappas=(1.0, 1.5, 2.0, 2.5, 3.0, 3.5)):
    """kappa at which the median split-half SD over DHSs crosses 1 (interpolated in log kappa)."""
    sd = [np.median([split_half_sd(*c, config.with_(kappa=k)) for c in counts]) for k in kappas]
    return float(np.exp(np.interp(0.0, np.log(sd[::-1]), np.log(kappas[::-1])))), dict(zip(kappas, sd))


def main(argv=None):
    p = argparse.ArgumentParser(description="Fit global constants of the nested model.")
    p.add_argument("--summaries", nargs="+", required=True, help="NPZ files from save_data_results")
    p.add_argument("--counts", nargs="*", default=[], help="NPZ files with obs, exp, sample_groups")
    p.add_argument("--tau2", type=float, nargs="+", default=[0.3, 0.35, 0.4, 0.5])
    p.add_argument("--omega-scale", type=float, nargs="+", default=[0.1, 0.15, 0.2, 0.3])
    p.add_argument("--out", required=True)
    a = p.parse_args(argv)

    config = DEFAULT_NESTED
    if a.counts:
        counts = [tuple(np.load(f)[k] for k in ("obs", "exp", "sample_groups")) for f in a.counts]
        kappa, curve = calibrate_kappa(counts, config)
        config = config.with_(kappa=round(kappa, 2))
        print("split-half SD by kappa:", curve)
    summaries = [(np.load(f)["summaries.x"], np.load(f)["summaries.v"]) for f in a.summaries]
    evidence = summed_log_evidence(summaries, config, {"tau2": a.tau2, "omega_scale": a.omega_scale})
    best = max(evidence, key=evidence.get)
    for point, z in sorted(evidence.items(), key=lambda kv: -kv[1])[:5]:
        print(dict(tau2=point[0], omega_scale=point[1]), round(z - evidence[best], 1))
    config.with_(tau2=best[0], omega_scale=best[1]).to_json(a.out)


if __name__ == "__main__":
    main()
