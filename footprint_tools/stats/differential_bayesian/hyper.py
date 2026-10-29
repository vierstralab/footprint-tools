"""Global constants by summed log evidence over DHSs.

The evidence factorizes into a consensus part (tau2, b_var, consensus lengths) and a group-pattern
part (omega_scale, group-pattern lengths), so each part is maximized on its own grid.

Per DHS (inside the pipeline):  evidence_tables(x, v, config), stored as evidence.consensus / evidence.pattern
Once, over all DHSs:            python -m footprint_tools.stats.differential_bayesian.hyper stats/*.npz --out model.json
"""

import argparse

import numpy as np

from .config import DEFAULT_MODEL, ModelConfig
from .model import FootprintModel, _forward, gamma_lengths, half_cauchy

# Block lengths are resolution choices (the evidence prefers blocks of a few bp in both layers), so the
# default grids hold them at the config's values; add cons_* / pattern_* keys to explore them.
GRID_CONSENSUS = dict(tau2=[0.2, 0.3, 0.4, 0.5, 0.7], b_var=[0.25, 0.5, 1.0, 2.0])
GRID_PATTERN = dict(omega_scale=[0.05, 0.1, 0.15, 0.2, 0.3, 0.5])
REBUILD = ("tau2", "b_var")   # the other constants only change the forward pass


def evidence_tables(x, v, base: ModelConfig = DEFAULT_MODEL, grid_b=GRID_CONSENSUS, grid_d=GRID_PATTERN):
    """Log evidence of one DHS at every grid point: (consensus array, group-pattern array), shaped like the grids."""
    lb, ld = list(grid_b), list(grid_d)
    z_b = np.empty([len(grid_b[k]) for k in lb])
    model = None
    for ib in np.ndindex(z_b.shape):
        cfg = base.with_(**{k: grid_b[k][i] for k, i in zip(lb, ib)})
        if model is None or any(getattr(model.config, k) != getattr(cfg, k) for k in REBUILD):
            model = FootprintModel(x, v, np.ones_like(x), cfg)
        z_b[ib] = _forward(model.E_b, gamma_lengths(model.Lm, cfg.cons_mean, cfg.cons_shape, cfg.cons_min))[-1] + model.const
    z_d = np.empty([len(grid_d[k]) for k in ld])
    E_d = {}
    for idd in np.ndindex(z_d.shape):
        cfg = base.with_(**{k: grid_d[k][i] for k, i in zip(ld, idd)})
        if cfg.omega_scale not in E_d:
            E_d = {cfg.omega_scale: np.logaddexp.reduce(model.E_d_omega + half_cauchy(model.omega2, cfg.omega_scale), axis=2)}
        z_d[idd] = _forward(E_d[cfg.omega_scale], gamma_lengths(model.Lm, cfg.pattern_mean, cfg.pattern_shape, cfg.pattern_min))[-1]
    return z_b, z_d


def fit(tables, base: ModelConfig = DEFAULT_MODEL, grid_b=GRID_CONSENSUS, grid_d=GRID_PATTERN):
    """Arg-max of the summed tables; returns the fitted config and the summed tables."""
    z_b, z_d = sum(t[0] for t in tables), sum(t[1] for t in tables)
    best = {**{k: grid_b[k][i] for k, i in zip(grid_b, np.unravel_index(z_b.argmax(), z_b.shape))},
            **{k: grid_d[k][i] for k, i in zip(grid_d, np.unravel_index(z_d.argmax(), z_d.shape))}}
    return base.with_(**best), z_b, z_d


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("stats", nargs="+")
    p.add_argument("--out", required=True)
    p.add_argument("--config", default="", help="ModelConfig JSON the tables were computed with (default: package defaults)")
    a = p.parse_args(argv)
    tables = []
    for f in a.stats:
        with np.load(f) as h:
            tables.append((h["evidence.consensus"], h["evidence.pattern"]))
    config, _, _ = fit(tables, ModelConfig.from_json(a.config) if a.config else DEFAULT_MODEL)
    config.to_json(a.out)
    print(config)


if __name__ == "__main__":
    main()
