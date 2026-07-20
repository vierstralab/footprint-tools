
from dataclasses import dataclass

import numpy as np
from numba import njit, prange
from scipy.special import logsumexp, ndtr

from .differential import Differential, _quadrature_weights


@dataclass(frozen=True, slots=True)
class GaussianSummary:
    mu: np.ndarray
    mu_var: np.ndarray
    sig2: np.ndarray
    log_evidence: np.ndarray


def gaussian_summary(data: Differential, variance_floor: float | None = None) -> GaussianSummary:
    log_width = np.log(_quadrature_weights(data.mu_x))
    joint = data.loglik_mu + log_width[None, None]
    evidence = logsumexp(joint, axis=2)
    mass = np.exp(joint - evidence[:, :, None])
    mu = np.sum(mass * data.mu_x[None, None], axis=2)
    mu_var = np.sum(mass * (data.mu_x[None, None] - mu[:, :, None]) ** 2, axis=2)
    floor = variance_floor
    if floor is None:
        floor = np.median(np.diff(data.mu_x)) ** 2 / 12
    mu_var = np.maximum(mu_var, floor)

    log_sig2 = (
        logsumexp(data.loglik_mu_sig2 + log_width[None, None, :, None], axis=2)
        + data.log_sig2_prior[:, None]
    )
    log_sig2 -= logsumexp(log_sig2, axis=2, keepdims=True)
    sig2 = np.sum(np.exp(log_sig2) * data.sig2_x[None, None], axis=2)
    return GaussianSummary(mu, mu_var, sig2, evidence)


def variance_ratio_group_terms(
    data: Differential,
    mu0_x: np.ndarray,
    eta_x: np.ndarray,
    method: str,
    variance_floor: float | None = None,
) -> np.ndarray:
    ratio_x = np.where(np.isneginf(eta_x), 0.0, np.exp(eta_x))
    if method == "gaussian":
        return _variance_ratio_gaussian(data, mu0_x, ratio_x, variance_floor)
    if method == "exact":
        return _variance_ratio_exact(data, mu0_x, ratio_x)
    raise ValueError("method must be 'gaussian' or 'exact'")


def coefficient_target_terms(
    data: Differential,
    mu0_x: np.ndarray,
    z_x: np.ndarray,
    method: str,
    variance_floor: float | None = None,
) -> np.ndarray:
    if method == "gaussian":
        s = gaussian_summary(data, variance_floor)
        mean = mu0_x[None, None, :, None] + z_x[None, None, None] * np.sqrt(s.sig2)[:, :, None, None]
        return (
            s.log_evidence[:, :, None, None]
            - 0.5 * (
                np.log(2 * np.pi * s.mu_var)[:, :, None, None]
                + (s.mu[:, :, None, None] - mean) ** 2 / s.mu_var[:, :, None, None]
            )
        )
    if method == "exact":
        out = np.empty((len(data.group_names), data.loglik_mu.shape[1], len(mu0_x), len(z_x)))
        for g in range(len(data.group_names)):
            out[g] = _target_exact(
                np.ascontiguousarray(data.loglik_mu_sig2[g]),
                np.ascontiguousarray(data.log_sig2_prior[g]),
                np.ascontiguousarray(data.mu_x),
                np.ascontiguousarray(data.sig2_x),
                np.ascontiguousarray(mu0_x),
                np.ascontiguousarray(z_x),
            )
        return out
    raise ValueError("method must be 'gaussian' or 'exact'")


def _variance_ratio_gaussian(data, mu0_x, ratio_x, variance_floor):
    s = gaussian_summary(data, variance_floor)
    g, n = s.mu.shape
    out = np.empty((g, n, len(mu0_x), len(ratio_x)))
    for k, ratio in enumerate(ratio_x):
        if ratio == 0:
            for i in range(g):
                out[i, :, :, k] = _interpolate_grid(data.loglik_mu[i], data.mu_x, mu0_x)
        else:
            variance = s.mu_var + ratio * s.sig2
            out[:, :, :, k] = (
                s.log_evidence[:, :, None]
                - 0.5 * (
                    np.log(2 * np.pi * variance)[:, :, None]
                    + (mu0_x[None, None] - s.mu[:, :, None]) ** 2 / variance[:, :, None]
                )
            )
    return out


def _variance_ratio_exact(data, mu0_x, ratio_x):
    g, n, _, _ = data.loglik_mu_sig2.shape
    out = np.empty((g, n, len(mu0_x), len(ratio_x)))
    for k, ratio in enumerate(ratio_x):
        if ratio == 0:
            for i in range(g):
                out[i, :, :, k] = _interpolate_grid(data.loglik_mu[i], data.mu_x, mu0_x)
            continue
        kernel = _normal_grid_mass(data.mu_x, mu0_x, ratio * data.sig2_x)
        with np.errstate(divide="ignore"):
            log_kernel = np.log(kernel)
        out[:, :, :, k] = _quadrature_state(
            np.ascontiguousarray(data.loglik_mu_sig2),
            np.ascontiguousarray(data.log_sig2_prior),
            np.ascontiguousarray(log_kernel),
        )
    return out


def _normal_grid_mass(support_x, mean_x, variance_x):
    bounds = np.r_[-np.inf, (support_x[:-1] + support_x[1:]) / 2, np.inf]
    z = (bounds[None, None] - mean_x[None, :, None]) / np.sqrt(variance_x)[:, None, None]
    mass = np.maximum(np.diff(ndtr(z), axis=2), 0.0)
    return mass / mass.sum(axis=2, keepdims=True)


def _interpolate_grid(loglik, source_x, target_x):
    if np.array_equal(source_x, target_x):
        return loglik.copy()
    return np.vstack([
        np.interp(target_x, source_x, row, left=-np.inf, right=-np.inf)
        for row in loglik
    ])


@njit(cache=True, parallel=True)
def _quadrature_state(raw, log_sig2_prior, log_kernel):
    g, n, m, s = raw.shape
    u = log_kernel.shape[1]
    out = np.empty((g, n, u))
    for flat in prange(g * n * u):
        ui = flat % u
        q = flat // u
        p = q % n
        gi = q // n
        maximum = -np.inf
        for si in range(s):
            for mi in range(m):
                v = raw[gi, p, mi, si] + log_sig2_prior[gi, si] + log_kernel[si, ui, mi]
                maximum = max(maximum, v)
        total = 0.0
        for si in range(s):
            for mi in range(m):
                v = raw[gi, p, mi, si] + log_sig2_prior[gi, si] + log_kernel[si, ui, mi]
                if np.isfinite(v):
                    total += np.exp(v - maximum)
        out[gi, p, ui] = maximum + np.log(total)
    return out


@njit(cache=True, parallel=True)
def _target_exact(raw, log_sig2_prior, mu_x, sig2_x, mu0_x, z_x):
    n, _, s = raw.shape
    out = np.empty((n, len(mu0_x), len(z_x)))
    for flat in prange(n * len(mu0_x) * len(z_x)):
        zi = flat % len(z_x)
        q = flat // len(z_x)
        ui = q % len(mu0_x)
        p = q // len(mu0_x)
        maximum = -np.inf
        values = np.empty(s)
        for si in range(s):
            x = mu0_x[ui] + z_x[zi] * np.sqrt(sig2_x[si])
            values[si] = _interp_1d(raw[p, :, si], mu_x, x) + log_sig2_prior[si]
            maximum = max(maximum, values[si])
        total = 0.0
        for si in range(s):
            if np.isfinite(values[si]):
                total += np.exp(values[si] - maximum)
        out[p, ui, zi] = maximum + np.log(total) if total > 0 else -np.inf
    return out


@njit(cache=True)
def _interp_1d(y, x, value):
    if value < x[0] or value > x[-1]:
        return -np.inf
    if value == x[-1]:
        return y[-1]
    lo, hi = 0, len(x) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if x[mid] <= value:
            lo = mid
        else:
            hi = mid
    if not np.isfinite(y[lo]) or not np.isfinite(y[hi]):
        return y[lo] if value == x[lo] else y[hi] if value == x[hi] else -np.inf
    w = (value - x[lo]) / (x[hi] - x[lo])
    return y[lo] + w * (y[hi] - y[lo])
