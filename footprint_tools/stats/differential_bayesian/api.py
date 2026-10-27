import warnings

import numpy as np
from genome_tools import GenomicInterval
from genome_tools.data.utils import realign_matrix
from genome_tools.plotting.modular_plot.api import DataBundle, PlotDataLoader

from .base_differential import Differential, DifferentialModel
from .config import DEFAULT_DIFFERENTIAL, DEFAULT_NESTED, DEFAULT_THETA, DEFAULT_THETA_SEGMENTATION
from .results import NestedFit, fit_nested
from .segmentation import Segmentation
from .theta import ThetaLikelihood, ThetaModel, fit_theta_segmentation


_DATA_RESULT_CLASSES = {
    "differential": Differential,
    "nested": NestedFit,
    "theta": ThetaLikelihood,
    "theta_segmentation": Segmentation,
}

_GENOMIC_AXES = {
    "differential": {"loglik_mu": 1, "mu_mean": 1, "mu_var": 1},
    "nested": {
        "depth_mean": 0, "depth_lo": 0, "depth_hi": 0, "p_any": 0, "delta_mean": 0,
        "p_more": 0, "p_less": 0, "icc": 0, "draws": 1, "p_outer": 0, "p_inner": 0,
    },
    "theta": {"loglik": 1},
    "theta_segmentation": {"log_posterior": 1, "boundary": 1},
}

_BOUNDARY_KEYS = {"boundary", "p_outer", "p_inner"}

_REMOVED_LOADERS = {
    "GroupMeanSegmentationLoader": "NestedSegmentationLoader",
    "VarianceRatioLoader": "NestedSegmentationLoader (ICC is NestedFit.icc)",
    "EtaSegmentationLoader": "NestedSegmentationLoader (ICC is NestedFit.icc)",
    "Mu0SegmentationLoader": "NestedSegmentationLoader",
    "CommonCoefficientLikelihoodLoader": "NestedSegmentationLoader",
    "CommonCoefficientSegmentationLoader": "NestedSegmentationLoader (calls are NestedFit.calls)",
    "ZeroCoefficientLikelihoodLoader": "NestedSegmentationLoader",
    "ZeroCoefficientSegmentationLoader": "NestedSegmentationLoader",
    "ZeroFootprintCountLoader": "NestedFit.p_any / NestedFit.footprints",
    "DeviantFootprintCountLoader": "NestedFit.p_more, NestedFit.p_less / NestedFit.calls",
}


def __getattr__(name):
    if name in _REMOVED_LOADERS:
        raise ImportError(f"{name} was removed with the nested model; use {_REMOVED_LOADERS[name]}")
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def data_results_to_dict(data, fields=None):
    """Serializable loader-produced fields present on ``data``."""
    out = {"__interval__": data.interval.to_ucsc()}
    for name in _selected_result_fields(fields):
        obj = getattr(data, name, None)
        if obj is not None:
            out[name] = obj.to_dict()
    return out


def data_results_from_dict(values, data=None, fields=None):
    """Attach serialized loader-produced fields to ``data`` and return it."""
    data = DataBundle(interval=values["__interval__"]) if data is None else data
    selected = set(values) if fields is None else set(_selected_result_fields(fields))
    for name, cls in _DATA_RESULT_CLASSES.items():
        if name in selected and name in values:
            setattr(data, name, cls.from_dict(values[name]))
    return data


def save_data_results(data, path, fields=None):
    """Save all present loader-produced result fields from ``data`` to one NPZ."""
    values = data_results_to_dict(data, fields)
    result_fields = tuple(k for k in values if k != "__interval__")
    payload = {
        "__fields__": np.asarray(result_fields, dtype=str),
        "__interval__": values["__interval__"],
    }
    for name in result_fields:
        payload[f"{name}.__class__"] = np.asarray(_DATA_RESULT_CLASSES[name].__name__)
        for key, value in values[name].items():
            payload[f"{name}.{key}"] = value
    np.savez_compressed(path, **payload)


def _slice_if_genomic(x, result_name, key, source_interval, data):
    if data is None:
        return x
    axis = _GENOMIC_AXES.get(result_name, {}).get(key)
    if axis is None:
        return x
    x = np.asarray(x)
    axis %= x.ndim
    if key in _BOUNDARY_KEYS:
        source_interval = source_interval.widen(left=0, right=-1)
        target_interval = data.interval.widen(left=0, right=-1)
    else:
        target_interval = data.interval
    assert x.shape[axis] == len(source_interval), (
        f"{result_name}.{key}: axis {axis} has length {x.shape[axis]}, "
        f"expected {len(source_interval)} for {source_interval}"
    )
    x = np.moveaxis(x, axis, -1)
    x = realign_matrix(x, source_interval, target_interval)
    return np.moveaxis(x, -1, axis)


def load_data_results(path, data=None, fields=None):
    """Load an NPZ written by `save_data_results` into `data` and return it."""
    with np.load(path, allow_pickle=False) as handle:
        interval = GenomicInterval.from_ucsc(handle["__interval__"].item())
        saved = tuple(str(x) for x in handle["__fields__"].tolist())
        unknown = [name for name in saved if name not in _DATA_RESULT_CLASSES]
        if unknown:
            warnings.warn(f"{path}: skipping results of removed stages {unknown}")
        selected = set(saved) - set(unknown) if fields is None else set(_selected_result_fields(fields))
        values = {"__interval__": interval}
        for name in saved:
            if name not in selected:
                continue
            prefix = f"{name}."
            values[name] = {
                key[len(prefix):]: _slice_if_genomic(handle[key], name, key[len(prefix):], interval, data)
                for key in handle.files
                if key.startswith(prefix) and not key.endswith(".__class__")
            }
        return data_results_from_dict(values, data)


def _selected_result_fields(fields):
    if fields is None:
        return tuple(_DATA_RESULT_CLASSES)
    unknown = set(fields) - set(_DATA_RESULT_CLASSES)
    if unknown:
        raise ValueError(f"unknown data result fields: {sorted(unknown)}")
    return tuple(fields)


class DifferentialLoader(PlotDataLoader):
    def _load(self, data, selected_groups=None, config=DEFAULT_DIFFERENTIAL):
        data.differential = DifferentialModel(config).fit(
            data.groups_data,
            data.obs,
            data.exp,
            data.disp_models,
            selected_groups,
        )
        return data


class NestedSegmentationLoader(PlotDataLoader):
    def _load(self, data, config=DEFAULT_NESTED):
        x, v = data.differential.moments()
        data.nested = fit_nested(x, v, data.differential.group_names, config)
        return data


class ThetaLoader(PlotDataLoader):
    def _load(self, data, config=DEFAULT_THETA, log_mu_prior=None):
        data.theta = ThetaModel(config).fit(
            data.groups_data,
            data.obs,
            data.exp,
            data.disp_models,
            data.differential,
            log_mu_prior,
        )
        return data


class ThetaSegmentationLoader(PlotDataLoader):
    def _load(self, data, length_prior, config=DEFAULT_THETA_SEGMENTATION, log_theta_prior=None):
        data.theta_segmentation = fit_theta_segmentation(
            data.theta,
            length_prior,
            config,
            log_theta_prior,
        )
        return data
