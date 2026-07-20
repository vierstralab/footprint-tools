from types import SimpleNamespace

import numpy as np
from genome_tools.plotting.modular_plot.api import PlotDataLoader

from .coefficients import (
    CoefficientLikelihood,
    CommonCoefficientModel,
    ZeroCoefficientModel,
    fit_coefficient_segmentation,
    infer_kfp_dev,
    infer_kfp_zero,
)
from .config import (
    DEFAULT_COEFFICIENT,
    DEFAULT_COEFFICIENT_SEGMENTATION,
    DEFAULT_DIFFERENTIAL,
    DEFAULT_ETA_SEGMENTATION,
    DEFAULT_MEAN_SEGMENTATION,
    DEFAULT_THETA,
    DEFAULT_THETA_SEGMENTATION,
    DEFAULT_VARIANCE_RATIO,
)
from .differential import Differential, DifferentialModel, fit_group_means_segmentation
from .eta import EtaSegmentation, fit_eta_segmentation
from .posterior import GridPosterior
from .segmentation import Segmentation
from .theta import ThetaLikelihood, ThetaModel, fit_theta_segmentation
from .variance_ratio import VarianceRatioLikelihood, VarianceRatioModel, fit_mu0_segmentation


_DATA_RESULT_CLASSES = {
    "differential": Differential,
    "segmentation": Segmentation,
    "variance_ratio": VarianceRatioLikelihood,
    "eta_segmentation": EtaSegmentation,
    "mu0_segmentation": Segmentation,
    "common_coefficient_likelihood": CoefficientLikelihood,
    "common_coefficient_segmentation": Segmentation,
    "zero_coefficient_likelihood": CoefficientLikelihood,
    "zero_coefficient_segmentation": Segmentation,
    "kfp_zero": GridPosterior,
    "kfp_dev": GridPosterior,
    "theta": ThetaLikelihood,
    "theta_segmentation": Segmentation,
}


def data_results_to_dict(data, fields=None):
    """Return serializable loader-produced fields from a data object.

    Only known fields that are present on ``data`` are included.  Each value is
    produced by the object's own ``to_dict()`` method, so likelihood/result
    classes retain ownership of their serialization.
    """
    selected = _selected_result_fields(fields)
    out = {}
    for name in selected:
        obj = getattr(data, name, None)
        if obj is not None:
            out[name] = obj.to_dict()
    return out


def data_results_from_dict(values, data=None, fields=None):
    """Attach serialized loader-produced fields to ``data`` and return it."""
    data = SimpleNamespace() if data is None else data
    selected = set(values) if fields is None else set(_selected_result_fields(fields))
    for name in _DATA_RESULT_CLASSES:
        if name in selected and name in values:
            setattr(data, name, _DATA_RESULT_CLASSES[name].from_dict(values[name]))
    return data


def save_data_results(data, path, fields=None):
    """Save all present loader-produced result fields from ``data`` to one NPZ."""
    values = data_results_to_dict(data, fields)
    payload = {"__fields__": np.asarray(tuple(values), dtype=str)}
    for name, item in values.items():
        payload[f"{name}.__class__"] = np.asarray(_DATA_RESULT_CLASSES[name].__name__)
        for key, value in item.items():
            payload[f"{name}.{key}"] = value
    np.savez_compressed(path, **payload)


def load_data_results(path, data=None, fields=None):
    """Load an NPZ written by ``save_data_results`` into ``data`` and return it."""
    with np.load(path, allow_pickle=False) as handle:
        saved = tuple(str(x) for x in handle["__fields__"].tolist())
        selected = set(saved) if fields is None else set(_selected_result_fields(fields))
        values = {}
        for name in saved:
            if name not in selected:
                continue
            prefix = f"{name}."
            values[name] = {
                key[len(prefix):]: handle[key]
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


class GroupMeanSegmentationLoader(PlotDataLoader):
    def _load(
        self,
        data,
        length_prior,
        config=DEFAULT_MEAN_SEGMENTATION,
        log_mu_prior=None,
    ):
        data.segmentation = fit_group_means_segmentation(
            data.differential,
            length_prior,
            config,
            log_mu_prior,
        )
        return data


class VarianceRatioLoader(PlotDataLoader):
    def _load(
        self,
        data,
        config=DEFAULT_VARIANCE_RATIO,
        log_mu0_prior=None,
        log_eta_prior=None,
    ):
        data.variance_ratio = VarianceRatioModel(config).fit(
            data.differential,
            log_mu0_prior,
            log_eta_prior,
        )
        return data


class EtaSegmentationLoader(PlotDataLoader):
    def _load(
        self,
        data,
        length_prior,
        config=DEFAULT_ETA_SEGMENTATION,
        log_mu0_prior=None,
        log_eta_prior=None,
    ):
        data.eta_segmentation = fit_eta_segmentation(
            data.variance_ratio,
            length_prior,
            config,
            log_mu0_prior,
            log_eta_prior,
        )
        return data


class Mu0SegmentationLoader(PlotDataLoader):
    def _load(
        self,
        data,
        length_prior,
        config=DEFAULT_MEAN_SEGMENTATION,
        log_mu0_prior=None,
        log_eta_prior=None,
    ):
        data.mu0_segmentation = fit_mu0_segmentation(
            data.variance_ratio,
            length_prior,
            config,
            log_mu0_prior,
            log_eta_prior,
        )
        return data


class CommonCoefficientLikelihoodLoader(PlotDataLoader):
    def _load(
        self,
        data,
        config=DEFAULT_COEFFICIENT,
        log_mu0_prior=None,
        log_eta_prior=None,
    ):
        data.common_coefficient_likelihood = CommonCoefficientModel(config).fit(
            data.differential,
            data.variance_ratio,
            log_mu0_prior,
            log_eta_prior,
        )
        return data


class CommonCoefficientSegmentationLoader(PlotDataLoader):
    def _load(
        self,
        data,
        length_prior,
        config=DEFAULT_COEFFICIENT_SEGMENTATION,
    ):
        data.common_coefficient_segmentation = fit_coefficient_segmentation(
            data.common_coefficient_likelihood,
            length_prior,
            config,
        )
        return data


class ZeroCoefficientLikelihoodLoader(PlotDataLoader):
    def _load(
        self,
        data,
        config=DEFAULT_COEFFICIENT,
        log_z_prior=None,
    ):
        data.zero_coefficient_likelihood = ZeroCoefficientModel(config).fit(
            data.differential,
            log_z_prior,
        )
        return data


class ZeroCoefficientSegmentationLoader(PlotDataLoader):
    def _load(
        self,
        data,
        length_prior,
        config=DEFAULT_COEFFICIENT_SEGMENTATION,
    ):
        data.zero_coefficient_segmentation = fit_coefficient_segmentation(
            data.zero_coefficient_likelihood,
            length_prior,
            config,
        )
        return data


class ZeroFootprintCountLoader(PlotDataLoader):
    def _load(self, data, threshold=1.0, source="segmented"):
        if source == "segmented":
            posterior = data.zero_coefficient_segmentation.posterior
        elif source == "pointwise":
            posterior = data.zero_coefficient_likelihood.posterior()
        else:
            raise ValueError("source must be 'segmented' or 'pointwise'")
        data.kfp_zero = infer_kfp_zero(posterior, threshold)
        return data


class DeviantFootprintCountLoader(PlotDataLoader):
    def _load(self, data, threshold=1.0, source="segmented"):
        if source == "segmented":
            posterior = data.common_coefficient_segmentation.posterior
        elif source == "pointwise":
            posterior = data.common_coefficient_likelihood.posterior()
        else:
            raise ValueError("source must be 'segmented' or 'pointwise'")
        data.kfp_dev = infer_kfp_dev(posterior, threshold)
        return data


class ThetaLoader(PlotDataLoader):
    def _load(
        self,
        data,
        config=DEFAULT_THETA,
        log_mu_prior=None,
    ):
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
    def _load(
        self,
        data,
        length_prior,
        config=DEFAULT_THETA_SEGMENTATION,
        log_theta_prior=None,
    ):
        data.theta_segmentation = fit_theta_segmentation(
            data.theta,
            length_prior,
            config,
            log_theta_prior,
        )
        return data
