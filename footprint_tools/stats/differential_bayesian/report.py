"""Genome-wide report from per-DHS results: Benjamini-Hochberg over all positions and over DHSs, Storey's pi0, calls,
tables and a summary page.  python -m footprint_tools.stats.differential_bayesian.report FILES_OR_FOLDERS... --prefix diff
Writes <prefix>.calls.tsv (one row per call), <prefix>.regions.tsv (one row per footprinted / not footprinted run of a
DHS), <prefix>.dhs.tsv, <prefix>.cutoffs.tsv, <prefix>.summary.pdf."""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .api import load_fit
from .calls import bh_cutoff, bh_q, storey_pi0

MORE, LESS = "#1f5bff", "#ff2d2d"


def dhs_id(path: Path) -> str:
    return path.name.removeprefix("diff_data.").removesuffix(".npz")


def report(paths: list[Path], prefix: str) -> None:
    rows, lp, xi, c95 = [], [], [], []
    for p in paths:                                                    # pass 1: everything except the calls
        interval, f = load_fit(p)
        rows.append(dict(dhs_id=dhs_id(p), interval=interval.to_ucsc(), width=f.p_any.size, footprint_bases=int(f.footprinted.sum()),
                         **{f"log10p_d{e}": float(x) for e, x in zip(f.cfg.report_d, f.lp_dhs)}))
        lp.append(f.lp); xi.append(f.xi); c95.append(f.c95)
    cfg = f.cfg; E = cfg.report_d
    if cfg.d_min not in E:
        raise ValueError("d_min must be one of report_d")
    lp = np.concatenate(lp, 1); dhs = pd.DataFrame(rows)
    cut = {e: bh_cutoff(lp[i], cfg.fdr) for i, e in enumerate(E)}
    calls, regions = [], []
    for p in paths:                                                    # pass 2: calls at the genome-wide cutoff
        interval, f = load_fit(p)
        for out, t in ((calls, f.calls(cut[cfg.d_min])), (regions, f.region_table(cut[cfg.d_min]))):
            t.insert(0, "interval", interval.to_ucsc()); t.insert(0, "dhs_id", dhs_id(p)); out.append(t)
    calls, regions = pd.concat(calls, ignore_index=True), pd.concat(regions, ignore_index=True)
    for e in E:
        dhs[f"log10q_d{e}"] = bh_q(dhs[f"log10p_d{e}"].to_numpy())
    dhs["n_calls"] = dhs.dhs_id.map(calls.dhs_id.value_counts()).fillna(0).astype(int)
    cutoffs = pd.DataFrame([dict(d=e, base_log10p_cutoff=cut[e], base_frac=(lp[i] <= cut[e]).mean(), base_pi0=storey_pi0(lp[i]),
                                 dhs_frac=(dhs[f"log10q_d{e}"] <= np.log10(cfg.fdr)).mean(), dhs_pi0=storey_pi0(dhs[f"log10p_d{e}"].to_numpy()))
                            for i, e in enumerate(E)])
    for name, t in (("calls", calls), ("regions", regions), ("dhs", dhs), ("cutoffs", cutoffs)):
        t.to_csv(f"{prefix}.{name}.tsv", sep="\t", index=False, float_format="%.4g")
    summary_page(cfg, lp, dhs, calls, cutoffs, np.concatenate(xi), np.concatenate(c95), f.group_names, f"{prefix}.summary.pdf")


def summary_page(cfg, lp, dhs, calls, cutoffs, xi, c95, groups, path):
    plt.rcParams.update({"font.size": 8.5, "axes.titlelocation": "left", "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42})
    E, i_min = cfg.report_d, cfg.report_d.index(cfg.d_min)
    fig, axs = plt.subplots(3, 3, figsize=(15, 13), gridspec_kw=dict(hspace=.55, wspace=.3, left=.06, right=.97, top=.92, bottom=.06))
    fig.suptitle(f"Differential footprints: {len(dhs):,} DHSs, {lp.shape[1]:,} positions; calls at d ≥ {cfg.d_min}, BH {cfg.fdr:g} over all positions "
                 f"(p ≤ 1e{cutoffs.base_log10p_cutoff[i_min]:.2f})", x=.06, ha="left", fontsize=12, weight="bold")
    ax = axs[0, 0]; lab = ["any" if e == 0 else f"d ≥ {e}" for e in E]
    ax.plot(lab, cutoffs.dhs_frac, "o-", c="k", label="DHSs: BH discoveries"); ax.plot(lab, 1 - cutoffs.dhs_pi0, "o:", c="k", label="DHSs: 1 − π0")
    ax.plot(lab, cutoffs.base_frac, "s-", c="0.55", label="positions: BH discoveries"); ax.plot(lab, 1 - cutoffs.base_pi0, "s:", c="0.55", label="positions: 1 − π0")
    ax.set_ylim(0, 1.03); ax.legend(fontsize=7, frameon=False); ax.set_title("Fraction differential by the effect required\n(above 'any', 1 − π0 is a lower bound)")
    bins = np.linspace(0, 1, 21)
    for ax, i, t in ((axs[0, 1], 0, "any difference"), (axs[0, 2], i_min, f"d ≥ {cfg.d_min}")):
        ax.hist(10 ** lp[i], bins=bins, density=True, color="0.6"); ax.axhline(cutoffs.base_pi0[i], c="k", ls="--", lw=1)
        ax.set_xlabel("p"); ax.set_title(f"Per-position p, null: {t} not reached\nπ0 = {cutoffs.base_pi0[i]:.2f}")
    ax = axs[1, 0]
    for e, col in ((E[0], "0.5"), (cfg.d_min, MORE)):
        ax.hist(10 ** dhs[f"log10p_d{e}"], bins=bins, density=True, histtype="step", lw=1.6, color=col, label="any" if e == 0 else f"d ≥ {e}")
    ax.legend(fontsize=7.5, frameon=False); ax.set_xlabel("p"); ax.set_title("Per-DHS p (best 5-bp window)")
    ax = axs[1, 1]; n = np.minimum(dhs.n_calls, 10); ax.bar(*np.unique(n, return_counts=True), color="0.6")
    ax.set_xlabel("calls per DHS (10 = 10+)"); ax.set_title(f"{(dhs.n_calls > 0).sum():,} DHSs with ≥ 1 call; {len(calls):,} calls")
    ax = axs[1, 2]; L = calls.end - calls.start
    ax.hist(np.minimum(L, 80), bins=np.arange(cfg.min_len, 82, 2), color="0.6"); ax.set_xlabel("call length (bp; 80 = 80+)")
    ax.set_title(f"Median length {L.median():.0f} bp; {(L * calls.in_footprint).sum() / L.sum():.0%} of called bases in footprints")
    ax = axs[2, 0]; y = np.arange(len(groups))
    for kind, col, dy in (("more", MORE, .2), ("less", LESS, -.2)):
        c = calls[calls.kind == kind]; ax.barh(y + dy, [(c.groups + ";" + c.partial.fillna("")).str.contains(g, regex=False).sum() for g in groups], .4, color=col, label=f"{kind} depleted")
    ax.set_yticks(y, [g.split("/")[0] for g in groups], fontsize=7); ax.invert_yaxis(); ax.legend(fontsize=7.5, frameon=False); ax.set_title("Groups in the called sets")
    ax = axs[2, 1]; ax.hist(np.clip(xi, -1, 1), bins=50, color="0.6"); ax.axvline(0, c="k", lw=.8)
    ax.set_xlabel("ξ"); ax.set_title(f"Generalised Pareto shape of the permutation maxima\n{(xi > 0).mean():.0%} > 0; ξ < 0 is used as 0")
    ax = axs[2, 2]; ax.hist(np.clip(c95, 0, 8), bins=50, color="0.6"); ax.set_xlabel("95% point of the permutation maximum")
    ax.set_title(f"Critical value per position: median {np.median(c95):.2f}")
    fig.savefig(path)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("inputs", nargs="+", type=Path); ap.add_argument("--prefix", default="diff")
    a = ap.parse_args()
    report(sorted(f for x in a.inputs for f in (x.glob("diff_data.*.npz") if x.is_dir() else [x])), a.prefix)
