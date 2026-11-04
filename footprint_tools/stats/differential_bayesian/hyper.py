"""Global constants over DHSs: tau2, b_var, omega_scale by summed log evidence, group nuggets t_g by moments.

The evidence factorizes into a consensus part (tau2, b_var, consensus lengths) and a group-pattern
part (omega_scale, group-pattern lengths), so each part is maximized on its own grid.
The nugget t_g is variance of group g's summaries beyond v that the shared tau2 does not cover: it is
fitted from base-to-base differences of the contrasts, where block levels and tau2 cancel.

Per DHS (inside the pipeline):  stats(x, v, group_names, config), stored in stats.<dhs>.npz
Once, over all DHSs:            python -m footprint_tools.stats.differential_bayesian.hyper stats/*.npz --out model.json
"""

import argparse
from collections.abc import Sequence

import numpy as np

from .config import DEFAULT_DF_MODEL, DifferentialFootprintModelConfig
from .model import DifferentialFootprintModel, _forward, gamma_lengths, half_cauchy

Grid = dict[str, list[float]]

# Block lengths are resolution choices (the evidence prefers blocks of a few bp in both layers), so the
# default grids hold them at the config's values; add cons_* / pattern_* keys to explore them.
GRID_CONSENSUS = dict(tau2=[0.2, 0.3, 0.4, 0.5, 0.6, 0.8], b_var=[0.25, 0.5, 1.0, 2.0])
GRID_PATTERN = dict(omega_scale=[0.05, 0.1, 0.15, 0.2, 0.3, 0.5])
REBUILD = ("tau2", "b_var")   # the other constants only change the forward pass


def stats(x: np.ndarray, v: np.ndarray, group_names: Sequence[str],
          base: DifferentialFootprintModelConfig = DEFAULT_DF_MODEL) -> dict[str, np.ndarray]:
    """Everything one DHS contributes to the global fit, as NPZ fields."""
    z_b, z_d = evidence_tables(x, np.asarray(v) + base.nugget(group_names)[:, None], base)
    num, A = nugget_stats(x, v)
    return {"evidence.consensus": z_b, "evidence.pattern": z_d, "nugget.num": num, "nugget.A": A,
            "group_names": np.asarray(group_names, str)}


def nugget_stats(x: np.ndarray, v: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Moment equations A t = num for the nuggets. c_gj = d_gj - sum_h w_h d_hj, d = base-to-base difference of x, has
    Var c_g = k_g (V_g + 2 t_g) + sum_h w_h^2 (V_h + 2 t_h), k = 1 - 2w, V = v_j + v_j+1; block edges are trimmed at 5 SD."""
    w = 1 / v.mean(1)
    w /= w.sum()
    d, V, k = np.diff(x, axis=1), v[:, 1:] + v[:, :-1], 1 - 2 * w
    c, pred = d - w @ d, k[:, None] * V + (w ** 2) @ V
    wt = np.where(np.abs(c) < 5 * np.sqrt(pred), 1 / pred ** 2, 0.0)
    W = wt.sum(1)
    return (wt * (c ** 2 - pred)).sum(1), np.diag(2 * k * W) + 2 * np.outer(W, w ** 2)


def fit_nugget(dhs_stats: Sequence[dict[str, np.ndarray]]) -> dict[str, float]:
    """Sum per-DHS moment equations over DHSs, aligning groups by name; t_g >= 0."""
    groups = sorted(set().union(*(h["group_names"].tolist() for h in dhs_stats)))
    num, A = np.zeros(len(groups)), np.zeros((len(groups),) * 2)
    for h in dhs_stats:
        i = [groups.index(g) for g in h["group_names"].tolist()]
        num[i] += h["nugget.num"]
        A[np.ix_(i, i)] += h["nugget.A"]
    return dict(zip(groups, np.maximum(np.linalg.solve(A, num), 0.0).round(4).tolist()))


def evidence_tables(x: np.ndarray, v: np.ndarray, base: DifferentialFootprintModelConfig = DEFAULT_DF_MODEL, grid_b: Grid = GRID_CONSENSUS,
                    grid_d: Grid = GRID_PATTERN) -> tuple[np.ndarray, np.ndarray]:
    """Log evidence of one DHS at every grid point: (consensus array, group-pattern array), shaped like the grids."""
    lb, ld = list(grid_b), list(grid_d)
    z_b = np.empty([len(grid_b[k]) for k in lb])
    model = None
    for ib in np.ndindex(z_b.shape):
        cfg = base.with_(**{k: grid_b[k][i] for k, i in zip(lb, ib)})
        if model is None or any(getattr(model.config, k) != getattr(cfg, k) for k in REBUILD):
            model = DifferentialFootprintModel(x, v, np.ones_like(x), cfg)
        z_b[ib] = _forward(model.E_b, gamma_lengths(model.Lm, cfg.cons_mean, cfg.cons_shape, cfg.cons_min))[-1] + model.const
    z_d = np.empty([len(grid_d[k]) for k in ld])
    E_d = {}
    for idd in np.ndindex(z_d.shape):
        cfg = base.with_(**{k: grid_d[k][i] for k, i in zip(ld, idd)})
        if cfg.omega_scale not in E_d:
            E_d = {cfg.omega_scale: np.logaddexp.reduce(model.E_d_omega + half_cauchy(model.omega2, cfg.omega_scale), axis=2)}
        z_d[idd] = _forward(E_d[cfg.omega_scale], gamma_lengths(model.Lm, cfg.pattern_mean, cfg.pattern_shape, cfg.pattern_min))[-1]
    return z_b, z_d


def fit(tables: Sequence[tuple[np.ndarray, np.ndarray]], base: DifferentialFootprintModelConfig = DEFAULT_DF_MODEL, grid_b: Grid = GRID_CONSENSUS,
        grid_d: Grid = GRID_PATTERN) -> tuple[DifferentialFootprintModelConfig, np.ndarray, np.ndarray]:
    """Arg-max of the summed tables; returns the fitted config and the summed tables."""
    z_b, z_d = sum(t[0] for t in tables), sum(t[1] for t in tables)
    best = {**{k: grid_b[k][i] for k, i in zip(grid_b, np.unravel_index(z_b.argmax(), z_b.shape))},
            **{k: grid_d[k][i] for k, i in zip(grid_d, np.unravel_index(z_d.argmax(), z_d.shape))}}
    return base.with_(**best), z_b, z_d


def main(argv: Sequence[str] | None = None) -> None:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("stats", nargs="+")
    p.add_argument("--out", required=True)
    p.add_argument("--config", default="", help="ModelConfig JSON the tables were computed with (default: package defaults)")
    a = p.parse_args(argv)
    dhs_stats = []
    for f in a.stats:
        with np.load(f) as h:
            dhs_stats.append({k: h[k] for k in h.files})
    tables = [(h["evidence.consensus"], h["evidence.pattern"]) for h in dhs_stats]
    config, _, _ = fit(tables, DifferentialFootprintModelConfig.from_json(a.config) if a.config else DEFAULT_DF_MODEL)
    config = config.with_(group_nugget=fit_nugget(dhs_stats))
    config.to_json(a.out)
    print(config)


if __name__ == "__main__":
    main()
