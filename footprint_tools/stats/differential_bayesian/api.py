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


def load_fit(path: str | Path) -> tuple[GenomicInterval, DifferentialFootprintFit]:
    with np.load(path, allow_pickle=False) as h:
        interval = GenomicInterval.from_ucsc(h["__interval__"].item())
        values = {f.name: h[f.name] for f in fields(DifferentialFootprintFit)}
    text = lambda x: str(x) if x.ndim == 0 else tuple(x.tolist())
    return interval, DifferentialFootprintFit(**{k: text(x) if x.dtype.kind == "U" else x for k, x in values.items()})


def load_data_results(path: str | Path, data: DataBundle | None = None) -> DataBundle:
    """Load results saved by `save_data_results`; per-position arrays are realigned to data.interval when given."""
    source, fit = load_fit(path)
    data = DataBundle(interval=source) if data is None else data
    n = len(source)
    realign = lambda x: realign_matrix(x, source, data.interval) if isinstance(x, np.ndarray) and x.ndim and x.shape[-1] == n else x
    data.footprints = DifferentialFootprintFit(**{f.name: realign(getattr(fit, f.name)) for f in fields(DifferentialFootprintFit)})
    return data
