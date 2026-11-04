"""Per-sample observation map rho_s(theta) = f_s + (1 - f_s) 2^(k_s theta): contrast k_s and floor f_s.

A sample sees the common depletion theta stretched by k_s and never below its floor f_s, the share of its
expected cleavage that binding does not prevent. Both are fitted once from binned profiles of many DHSs.

Per DHS (inside the pipeline):  profile_stats(obs, exp, sample_ids, labels), stored as O, E, sample_ids
Once, over all DHSs:            python -m footprint_tools.stats.differential_bayesian.samples profile.*.npz --out sample_maps.tsv
"""

import argparse
from collections.abc import Sequence

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

BINS = np.arange(-6.0, 3.01, 0.25)  # reference-depth bins, log2
MIN_EXP = 4                         # a sample's bases enter its profile only where it expects >= 4 cuts
DEEP = -2.0                         # the floor is identified in bins below this depth


def rho(theta: np.ndarray | float, k: np.ndarray | float, f: np.ndarray | float) -> np.ndarray:
    """Share of the expected cleavage a sample shows at depletion theta (log2)."""
    return f + (1 - f) * np.exp2(k * theta)


def halves(sample_ids: np.ndarray, labels: np.ndarray) -> np.ndarray:
    """Each group split in two, alternating over sorted sample ids (the same split in every DHS)."""
    half = np.zeros(len(sample_ids), int)
    for g in np.unique(labels):
        i = np.flatnonzero(labels == g)
        half[i[np.argsort(sample_ids[i])][1::2]] = 1
    return half


def profile_stats(obs: np.ndarray, exp: np.ndarray, sample_ids: Sequence[str],
                  labels: Sequence[str]) -> tuple[np.ndarray, np.ndarray]:
    """Observed and expected cuts per sample and reference-depth bin (samples x bins). A base's reference depth is
    log2 obs/exp in a 3-bp window summed over the other half of the samples, so a sample's noise does not pick its bin."""
    obs, exp = np.asarray(obs, float), np.asarray(exp, float)
    half, nb, win = halves(np.asarray(sample_ids, str), np.asarray(labels, str)), BINS.size - 1, np.ones(3)
    O, E = np.zeros((len(obs), nb)), np.zeros((len(obs), nb))
    for h in (0, 1):
        ro, re = (np.convolve(a[half != h].sum(0), win, "same") for a in (obs, exp))
        b = np.digitize(np.log2((ro + 0.5) / (re + 0.5)), BINS) - 1
        ok = (re > 5) & (b >= 0) & (b < nb)
        for s in np.flatnonzero(half == h):
            m = ok & (exp[s] >= MIN_EXP)
            O[s] += np.bincount(b[m], obs[s, m], nb)
            E[s] += np.bincount(b[m], exp[s, m], nb)
    return O, E


def fit_maps(O: np.ndarray, E: np.ndarray, n_iter: int = 12) -> tuple[np.ndarray, np.ndarray]:
    """k_s, f_s from profiles summed over DHSs: y_sb = log2(O/E) ~ log2 rho(d_b; k_s, f_s), with a latent depth d_b per
    bin; alternating weighted least squares, median k_s = 1. f_s is relative to the least-floored samples."""
    keep = (E.sum(0) > 5e3) & ((O > 0).mean(0) > 0.9)
    O, E = O[:, keep], E[:, keep]
    w = np.where(O > 0, 1 / (1 / ((O + 0.5) * np.log(2) ** 2) + 0.005), 0.0)
    y = np.where(w > 0, np.log2(np.maximum(O, 1e-9) / E), 0.0)
    d, k, f = np.log2(O.sum(0) / E.sum(0)), np.ones(len(O)), np.full(len(O), 0.02)
    curve = lambda d, k, f: np.log2(rho(d, k, f))
    for _ in range(n_iter):
        for s in range(len(O)):
            k[s], f[s] = least_squares(lambda p: np.sqrt(w[s]) * (curve(d, *p) - y[s]), [k[s], f[s]],
                                       bounds=([0.2, 0.0], [3.0, 0.6])).x
        m = np.median(k)
        k, d = k / m, d * m
        for b in range(d.size):
            d[b] = least_squares(lambda x: np.sqrt(w[:, b]) * (curve(x, k, f) - y[:, b]), [d[b]], bounds=(-9, 4)).x[0]
    return k, f


def load_maps(path: str) -> pd.DataFrame:
    """TSV with sample_id, k, f; indexed by sample_id."""
    return pd.read_table(path, index_col="sample_id")


def main(argv: Sequence[str] | None = None) -> None:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("profiles", nargs="+")
    p.add_argument("--out", required=True)
    p.add_argument("--deep_quantile", type=float, default=0.5,
                   help="use DHSs above this quantile of deep-footprint content (the floor is seen only in deep bins)")
    a = p.parse_args(argv)
    rows = []
    for path in a.profiles:
        with np.load(path) as h:
            E = h["E"]
            rows.append((pd.DataFrame(h["O"], index=h["sample_ids"]), pd.DataFrame(E, index=h["sample_ids"]),
                         E[:, BINS[:-1] < DEEP].sum() / max(E.sum(), 1)))
    cut = np.quantile([r[2] for r in rows], a.deep_quantile)
    O, E = (pd.concat([r[i] for r in rows if r[2] >= cut]).groupby(level=0).sum() for i in (0, 1))
    k, f = fit_maps(O.to_numpy(), E.loc[O.index].to_numpy())
    pd.DataFrame(dict(k=k, f=f), index=pd.Index(O.index, name="sample_id")).to_csv(a.out, sep="\t")
    print(f"{sum(r[2] >= cut for r in rows)} of {len(rows)} DHSs; median k = {np.median(k):.2f}, "
          f"median f = {np.median(f):.3f}, 95th percentile f = {np.quantile(f, 0.95):.3f}")


if __name__ == "__main__":
    main()
