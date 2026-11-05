"""DataBundle loader and NPZ storage of per-DHS results."""

from dataclasses import fields
from pathlib import Path

import numpy as np
from genome_tools import GenomicInterval
from genome_tools.data.utils import realign_matrix
from genome_tools.plotting.modular_plot.api import DataBundle, PlotDataLoader

from .config import DEFAULT_DF_CONFIG, DifferentialFootprintConfig
from .results import DifferentialFootprintFit, fit_differential_footprints


class DifferentialFootprintLoader(PlotDataLoader):
    """Needs data.groups_data, data.obs, data.exp, data.disp_models."""

    def _load(self, data: DataBundle, config: DifferentialFootprintConfig = DEFAULT_DF_CONFIG) -> DataBundle:
        data.footprints = fit_differential_footprints(data.groups_data, data.obs, data.exp, data.disp_models, config)
        return data


def save_data_results(data: DataBundle, path: str | Path) -> None:
    payload = {f.name: np.asarray(getattr(data.footprints, f.name)) for f in fields(DifferentialFootprintFit)}
    np.savez_compressed(path, __interval__=np.asarray(data.interval.to_ucsc()), **payload)


def load_data_results(path: str | Path, data: DataBundle | None = None) -> DataBundle:
    """Load results saved by `save_data_results`; per-position arrays are realigned to data.interval when given."""
    with np.load(path, allow_pickle=False) as h:
        source = GenomicInterval.from_ucsc(h["__interval__"].item())
        data = DataBundle(interval=source) if data is None else data
        data.footprints = DifferentialFootprintFit(**{f.name: _value(h[f.name], source, data.interval) for f in fields(DifferentialFootprintFit)})
    return data


def _value(x: np.ndarray, source: GenomicInterval, target: GenomicInterval) -> np.ndarray | str | tuple:
    if x.dtype.kind == "U":
        return str(x) if x.ndim == 0 else tuple(x.tolist())
    return realign_matrix(x, source, target)       # positions are the last axis
