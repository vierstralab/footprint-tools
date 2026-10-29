"""Sequence correction of the 6-mer expected cleavage, learned across DHSs.

The base-level residual r_j = high-pass of log2((sum obs + .5) / (sum exp + .5)), sums over all samples, is shared by
all groups and predicted by local sequence. A ridge regression on one-hot bases and dinucleotides at offsets -W..W is
fitted over all DHSs; the corrected expected is exp * 2^correction(seq).

Per DHS (inside the pipeline):  obs and exp sums over samples and the sequence of [start - W, end + W + 1)
Once, over all DHSs:            python -m footprint_tools.stats.differential_bayesian.seqcorr stats/*.npz --out seqcorr.json
"""

import argparse
import json

import numpy as np
from scipy.ndimage import uniform_filter1d

W, HIGHPASS, LAM = 6, 11, 30.0


def features(seq):
    """One-hot bases and dinucleotides at offsets -W..W, one row per base; seq spans [start - W, end + W + 1)."""
    c = np.array(["ACGT".find(b) for b in seq.upper()])
    mono = c[:-1, None] == np.arange(4)
    di = np.where((c[:-1] >= 0) & (c[1:] >= 0), 4 * c[:-1] + c[1:], -1)[:, None] == np.arange(16)
    f = np.hstack([mono, di])
    n = len(f) - 2 * W
    return np.hstack([f[o:o + n] for o in range(2 * W + 1)]).astype(float)


def residual(obs_sum, exp_sum):
    lr = np.log2((obs_sum + 0.5) / (exp_sum + 0.5))
    return lr - uniform_filter1d(lr, HIGHPASS, mode="nearest")


def correction(seq, coef):
    """log2 multiplier of the expected at each base."""
    return coef["intercept"] + features(seq) @ np.array(coef["beta"])


def fit(stats, lam=LAM):
    """Ridge fit from (chrom, obs_sum, exp_sum, seq) per DHS. Returns coefficients and leave-one-chromosome-out R^2."""
    acc = {}
    for chrom, o, e, seq in stats:
        X = np.hstack([np.ones((len(o), 1)), features(seq)])
        y = residual(o, e)
        a = acc.setdefault(chrom, [0.0, 0.0, 0.0])
        a[0], a[1], a[2] = a[0] + X.T @ X, a[1] + X.T @ y, a[2] + y @ y
    pen = lam * np.diag(np.r_[0.0, np.ones(len(next(iter(acc.values()))[1]) - 1)])
    solve = lambda A, b: np.linalg.lstsq(A + pen, b, rcond=None)[0]
    A, b = sum(a[0] for a in acc.values()), sum(a[1] for a in acc.values())
    sse = sst = 0.0
    for Ac, bc, yc in acc.values():
        beta = solve(A - Ac, b - bc)
        sse += yc - 2 * beta @ bc + beta @ Ac @ beta
        sst += yc - bc[0] ** 2 / Ac[0, 0]
    beta = solve(A, b)
    return dict(intercept=float(beta[0]), beta=beta[1:].tolist()), 1 - sse / sst


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("stats", nargs="+")
    p.add_argument("--out", required=True)
    a = p.parse_args(argv)
    stats = []
    for f in a.stats:
        with np.load(f) as h:
            stats.append((str(h["chrom"]), h["seqcorr.obs"], h["seqcorr.exp"], str(h["seqcorr.seq"])))
    coef, r2 = fit(stats)
    with open(a.out, "w") as handle:
        json.dump(coef, handle)
    print(f"{len(stats)} DHSs, leave-one-chromosome-out R^2 = {r2:.3f}")


if __name__ == "__main__":
    main()
