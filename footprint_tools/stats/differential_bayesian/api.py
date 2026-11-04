"""DataBundle loaders and NPZ storage of per-DHS results."""

from dataclasses import fields
from pathlib import Path

import numpy as np
import pandas as pd
from genome_tools import GenomicInterval
from genome_tools.data.utils import realign_matrix
from genome_tools.plotting.modular_plot.api import DataBundle, PlotDataLoader

from .config import DEFAULT_DF_MODEL, DEFAULT_DF_SUMMARY, DifferentialFootprintModelConfig, DifferentialFootprintSummaryConfig
from .results import DifferentialFootprintFit, fit_differential_footprints
from .summaries import Summaries, fit_summaries

RESULTS = {"summaries": Summaries, "footprints": DifferentialFootprintFit}
POSITION_AXIS = {**{k: 1 for k in ("x", "v", "sigma2", "theta_mean", "theta_sd")},
                 **{k: 0 for k in ("depth_mean", "depth_lo", "depth_hi", "p_any", "delta_mean", "p_more", "p_less", "icc",
                                   "p_pattern", "p_consensus")}}
BOUNDARY = {"p_pattern", "p_consensus"}


class DifferentialFootprintSummaryLoader(PlotDataLoader):
    """Needs data.groups_data, data.obs, data.exp, data.disp_models; maps: sample -> k, f (samples.load_maps)."""

    def _load(self, data: DataBundle, config: DifferentialFootprintSummaryConfig = DEFAULT_DF_SUMMARY,
              maps: pd.DataFrame | None = None) -> DataBundle:
        data.summaries = fit_summaries(data.groups_data, data.obs, data.exp, data.disp_models, config, maps)
        return data


class DifferentialFootprintLoader(PlotDataLoader):
    def _load(self, data: DataBundle, config: DifferentialFootprintModelConfig = DEFAULT_DF_MODEL) -> DataBundle:
        s = data.summaries
        data.footprints = fit_differential_footprints(s.x, s.v, s.sigma2, s.group_names, config)
        return data


def save_data_results(data: DataBundle, path: str | Path) -> None:
    payload = {"__interval__": np.asarray(data.interval.to_ucsc())}
    for name in RESULTS:
        if getattr(data, name, None) is not None:
            for f in fields(RESULTS[name]):
                payload[f"{name}.{f.name}"] = np.asarray(getattr(getattr(data, name), f.name))
    np.savez_compressed(path, **payload)


def load_data_results(path: str | Path, data: DataBundle | None = None) -> DataBundle:
    """Load results saved by `save_data_results`; arrays are realigned to data.interval when given."""
    with np.load(path, allow_pickle=False) as h:
        source = GenomicInterval.from_ucsc(h["__interval__"].item())
        data = DataBundle(interval=source) if data is None else data
        for name, cls in RESULTS.items():
            if f"{name}.{fields(cls)[0].name}" in h.files:
                values = {f.name: _value(h[f"{name}.{f.name}"], f.name, source, data.interval) for f in fields(cls)}
                setattr(data, name, cls(**values))
    return data


def _value(x: np.ndarray, key: str, source: GenomicInterval, target: GenomicInterval) -> np.ndarray | float | str | tuple:
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
