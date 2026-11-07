"""Footprints and differential calls from per-position statistics, and the genome-wide multiple-testing helpers."""

import numpy as np
from scipy.special import ndtr


def runs(mask: np.ndarray, min_len: int = 1) -> list[tuple[int, int]]:
    """Half-open (start, end) runs of True of at least min_len in a 1-D boolean array."""
    edges = np.flatnonzero(np.diff(np.r_[0, np.asarray(mask, int), 0]))
    return [(s, e) for s, e in zip(edges[::2].tolist(), edges[1::2].tolist()) if e - s >= min_len]


def p_footprint(m: np.ndarray, v: np.ndarray, idx: list[np.ndarray], fold: float) -> tuple[np.ndarray, np.ndarray]:
    """Group mean depth (groups x positions) and P(at least one group's mean depth < -fold), groups as independent
    normals with standard error sqrt((var of m + mean of v) / n_g)."""
    mean = np.stack([m[i].mean(0) for i in idx])
    se = np.stack([np.sqrt((m[i].var(0) + v[i].mean(0)) / len(i)) for i in idx])
    return mean, 1 - np.prod(ndtr((mean + fold) / se), 0)


def bh_cutoff(lp: np.ndarray, q: float) -> float:
    """log10 of the largest p-value Benjamini-Hochberg rejects at level q (lp: log10 p-values)."""
    s = np.sort(lp); k = np.flatnonzero(s <= np.log10(q * np.arange(1, s.size + 1) / s.size))
    return float(s[k[-1]]) if k.size else -np.inf


def bh_q(lp: np.ndarray) -> np.ndarray:
    """log10 Benjamini-Hochberg q-values."""
    o = np.argsort(lp); q = lp[o] + np.log10(lp.size / np.arange(1, lp.size + 1))
    out = np.empty(lp.size); out[o] = np.minimum(np.minimum.accumulate(q[::-1])[::-1], 0)
    return out


def storey_pi0(lp: np.ndarray, lam: float = .5) -> float:
    return min(1., float((lp > np.log10(lam)).mean() / (1 - lam)))


def differential_calls(lp: dict, sets: dict, w: dict, set_list: list[tuple[int, ...]], cutoff: float, min_len: int) -> list[dict]:
    """Calls from the best set per direction ('more', 'less') at every position: runs of >= min_len positions passing
    the cutoff. A group belongs to a call where it is in the best set for >= min_len consecutive positions. Where both
    directions pass and their sets share a group, the weaker direction gives way."""
    ok = {d: lp[d] <= cutoff for d in lp}
    clash = ok['more'] & ok['less'] & np.array([bool(set(set_list[a]) & set(set_list[b])) for a, b in zip(sets['more'], sets['less'])])
    weaker_more = lp['more'] > lp['less']
    ok['more'] &= ~(clash & weaker_more); ok['less'] &= ~(clash & ~weaker_more)
    out = []
    for d in ('more', 'less'):
        for s, e in runs(ok[d], min_len):
            members = [set_list[i] for i in sets[d][s:e]]
            spans = {k: [(s + a, s + b) for a, b in runs([k in x for x in members], min_len)] for k in sorted(set().union(*members))}
            spans = {k: v for k, v in spans.items() if v} or {k: [(s, e)] for k in set_list[np.bincount(sets[d][s:e]).argmax()]}
            out.append(dict(kind=d, start=s, end=e, spans=spans, w=float(w[d][s:e].mean()), min_log10p=float(lp[d][s:e].min())))
    return sorted(out, key=lambda c: c['start'])
