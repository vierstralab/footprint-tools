"""DataBundle loaders and NPZ storage of per-DHS results."""

from dataclasses import fields

import numpy as np
from genome_tools import GenomicInterval
from genome_tools.data.utils import realign_matrix
from genome_tools.plotting.modular_plot.api import DataBundle, PlotDataLoader

from .config import DEFAULT_NESTED, DEFAULT_SUMMARY
from .results import NestedFit, fit_nested
from .summaries import Summaries, fit_summaries

RESULTS = {"summaries": Summaries, "nested": NestedFit}
POSITION_AXIS = {"x": 1, "v": 1, "theta_mean": 1, "theta_sd": 1, "p_outer": 0, "p_inner": 0,
                 **{k: 0 for k in ("depth_mean", "depth_lo", "depth_hi", "p_any", "delta_mean", "p_more", "p_less", "icc")}}
BOUNDARY = {"p_outer", "p_inner"}


class SummaryLoader(PlotDataLoader):
    """Needs data.groups_data, data.obs, data.exp, data.disp_models."""

    def _load(self, data, config=DEFAULT_SUMMARY):
        data.summaries = fit_summaries(data.groups_data, data.obs, data.exp, data.disp_models, config)
        return data


class NestedLoader(PlotDataLoader):
    def _load(self, data, config=DEFAULT_NESTED):
        s = data.summaries
        data.nested = fit_nested(s.x, s.v, s.group_names, config)
        return data


def save_data_results(data, path):
    payload = {"__interval__": np.asarray(data.interval.to_ucsc())}
    for name in RESULTS:
        if getattr(data, name, None) is not None:
            for f in fields(RESULTS[name]):
                payload[f"{name}.{f.name}"] = np.asarray(getattr(getattr(data, name), f.name))
    np.savez_compressed(path, **payload)


def load_data_results(path, data=None):
    """Load results saved by `save_data_results`; arrays are realigned to data.interval when given."""
    with np.load(path, allow_pickle=False) as h:
        source = GenomicInterval.from_ucsc(h["__interval__"].item())
        data = DataBundle(interval=source) if data is None else data
        for name, cls in RESULTS.items():
            if f"{name}.{fields(cls)[0].name}" in h.files:
                values = {f.name: _value(h[f"{name}.{f.name}"], f.name, source, data.interval) for f in fields(cls)}
                setattr(data, name, cls(**values))
    return data


def _value(x, key, source, target):
    if x.dtype.kind == "U":
        return str(x) if x.ndim == 0 else tuple(x.tolist())
    if x.ndim == 0:
        return float(x)
    if key not in POSITION_AXIS:
        return x
    if key in BOUNDARY:
        source, target = source.widen(left=0, right=-1), target.widen(left=0, right=-1)
    axis = POSITION_AXIS[key]
    return np.moveaxis(realign_matrix(np.moveaxis(x, axis, -1), source, target), -1, axis)
