from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from typing import ClassVar

import numpy as np
from numba import njit
from scipy.special import logsumexp

from .io import Serializable, tuple_str
from .posterior import GridPosterior


_INDEPENDENT = 0
_INDEPENDENT_NO_SELF = 1
_SPARSE = 2


@dataclass(frozen=True, slots=True)
class LengthPrior(Serializable):
    """Explicit length weights with an optional inferred geometric tail."""

    save_attrs: ClassVar[tuple[str, ...]] = ("log_mass", "log_tail_ratio")

    log_mass: np.ndarray
    log_tail_ratio: np.ndarray | float = -np.inf

    def __post_init__(self) -> None:
        value = np.asarray(self.log_mass, dtype=float)
        ratio = np.asarray(self.log_tail_ratio, dtype=float)
        if value.ndim == 1:
            value = value[None]
        if ratio.ndim == 0:
            ratio = np.full(value.shape[0], ratio)
        weight = np.exp(value - np.max(value, axis=1, keepdims=True))
        r = np.exp(ratio)
        if np.any(r >= 1):
            raise ValueError("tail ratios must be below one")
        total = weight.sum(axis=1) + weight[:, -1] * r / (1 - r)
        with np.errstate(divide="ignore"):
            value = np.log(weight / total[:, None])
        object.__setattr__(self, "log_mass", np.ascontiguousarray(value))
        object.__setattr__(self, "log_tail_ratio", np.ascontiguousarray(ratio))

    @classmethod
    def from_log_mass(
        cls,
        log_mass: np.ndarray,
        infer_tail: bool = False,
        n_tail: int = 2,
    ) -> "LengthPrior":
        value = np.asarray(log_mass, dtype=float)
        rows = value[None] if value.ndim == 1 else value

        ratio = np.full(rows.shape[0], -np.inf)

        if infer_tail:
            if not 2 <= n_tail <= rows.shape[1]:
                raise ValueError(
                    f"n_tail must be between 2 and {rows.shape[1]}"
                )

            ratio = np.diff(rows[:, -n_tail:], axis=1).mean(axis=1)

            if np.any(ratio >= 0):
                raise ValueError("the inferred tail log-ratio must be below zero")

        return cls(value, ratio)


    def for_states(self, n_state: int) -> tuple[np.ndarray, np.ndarray]:
        mass = self.log_mass
        ratio = self.log_tail_ratio
        if mass.shape[0] == 1:
            mass = np.broadcast_to(mass, (n_state, mass.shape[1]))
            ratio = np.broadcast_to(ratio, n_state)
        if mass.shape[0] != n_state:
            raise ValueError("length-prior state dimension does not match")
        return np.ascontiguousarray(mass), np.ascontiguousarray(ratio)


@dataclass(frozen=True, slots=True)
class Segmentation(Serializable):
    save_attrs: ClassVar[tuple[str, ...]] = (
        "names", "state_x", "log_posterior", "boundary", "log_partition"
    )

    names: tuple[str, ...]
    posterior: GridPosterior
    boundary: np.ndarray
    log_partition: np.ndarray
    _sampler: "_SegmentationSampler | None" = field(
        default=None, repr=False, compare=False
    )

    def sample(
        self,
        n_draws: int = 1,
        rng: np.random.Generator | None = None,
    ) -> dict[str, tuple["SegmentationPath", ...]]:
        """Draw independent exact paths from the fitted posterior."""
        if self._sampler is None:
            raise RuntimeError(
                "sampling is unavailable after loading a summary-only NPZ"
            )
        return self._sampler.sample(n_draws, rng)

    def sample_prior(
        self,
        n_draws: int = 1,
        rng: np.random.Generator | None = None,
    ) -> dict[str, tuple["SegmentationPath", ...]]:
        """Draw stationary finite-window paths from the fitted prior."""
        if self._sampler is None:
            raise RuntimeError(
                "sampling is unavailable after loading a summary-only NPZ"
            )
        return self._sampler.sample_prior(n_draws, rng)

    def to_dict(self) -> dict[str, object]:
        return {
            "names": np.asarray(self.names, dtype=str),
            "state_x": self.posterior.x,
            "log_posterior": self.posterior.log_mass,
            "boundary": self.boundary,
            "log_partition": self.log_partition,
        }

    @classmethod
    def from_dict(cls, data: dict[str, object]) -> "Segmentation":
        return cls(
            tuple_str(data["names"]),
            GridPosterior(data["state_x"], data["log_posterior"]),
            data["boundary"], data["log_partition"],
        )


@dataclass(frozen=True, slots=True)
class SegmentationPath:
    starts: np.ndarray
    ends: np.ndarray
    state_index: np.ndarray
    state_value: np.ndarray

    @property
    def n_segments(self) -> int:
        return self.starts.size

    def dense(self) -> np.ndarray:
        return np.repeat(self.state_value, self.ends - self.starts)

    def dense_index(self) -> np.ndarray:
        return np.repeat(self.state_index, self.ends - self.starts)


def segment(
    loglik: np.ndarray,
    state_x: np.ndarray,
    names: tuple[str, ...],
    length_prior: LengthPrior,
    log_state_prior: np.ndarray,
    transition_sd: float | None = None,
    forbid_same_state: bool = True,
) -> Segmentation:
    """Exactly marginalize segmentations and retain traceback messages."""
    loglik, state_x = _track_inputs(loglik, state_x, names)
    g, n, m = loglik.shape
    prior = _normalize_track_prior(log_state_prior, g, m)
    log_q, log_tail = length_prior.for_states(m)
    mean, survival, equilibrium = _length_terms(log_q, log_tail, n)

    def run_track(i):
        model = _transition_model(
            state_x, prior[i], transition_sd, forbid_same_state
        )
        stationary_mean = float(np.exp(model[-1]) @ mean)
        z, posterior, boundary, prefix, reverse = _partition(
            np.ascontiguousarray(loglik[i]), log_q, log_tail,
            survival, equilibrium, prior[i], *model[:-1], model[-1],
            stationary_mean,
        )
        return z, posterior, boundary, (model, prefix, reverse)

    results = [run_track(0)]
    if g > 1:
        with ThreadPoolExecutor(max_workers=min(g - 1, 32)) as pool:
            results.extend(pool.map(run_track, range(1, g)))
    z = np.asarray([x[0] for x in results])
    posterior = np.asarray([x[1] for x in results])
    boundary = np.asarray([x[2] for x in results])
    sampler = _SegmentationSampler(
        state_x, names, prior, log_q, log_tail, survival, equilibrium, mean,
        tuple(x[3] for x in results),
    )
    return Segmentation(
        names, GridPosterior(state_x, posterior), boundary, z, sampler
    )


@dataclass(frozen=True, slots=True)
class _SegmentationSampler:
    state_x: np.ndarray
    save_attrs: ClassVar[tuple[str, ...]] = (
        "names", "state_x", "log_posterior", "boundary", "log_partition"
    )

    names: tuple[str, ...]
    prior: np.ndarray
    log_q: np.ndarray
    log_tail: np.ndarray
    survival: np.ndarray
    equilibrium: np.ndarray
    mean: np.ndarray
    tracks: tuple[tuple, ...]

    def sample(self, n_draws=1, rng=None):
        if n_draws < 1:
            raise ValueError("n_draws must be positive")
        rng = np.random.default_rng() if rng is None else rng
        result = {}
        for i, name in enumerate(self.names):
            model, prefix, reverse = self.tracks[i]
            outgoing = _outgoing_messages(reverse, self.prior[i], *model[:-1])
            result[name] = _posterior_paths(
                self.state_x, self.prior[i], self.log_q, self.log_tail,
                self.survival, self.equilibrium, self.mean, model, prefix,
                reverse, outgoing, n_draws, rng,
            )
        return result

    def sample_prior(self, n_draws=1, rng=None):
        if n_draws < 1:
            raise ValueError("n_draws must be positive")
        rng = np.random.default_rng() if rng is None else rng
        n = self.tracks[0][1].shape[0] - 1
        return {
            name: _prior_paths(
                self.state_x, self.prior[i], self.log_q, self.log_tail,
                self.survival, self.equilibrium, self.mean,
                self.tracks[i][0], n, n_draws, rng,
            )
            for i, name in enumerate(self.names)
        }

def _track_inputs(loglik, state_x, names):
    loglik = np.asarray(loglik, dtype=float)
    state_x = np.asarray(state_x, dtype=float)
    if loglik.ndim == 2:
        loglik = loglik[None]
    if loglik.ndim != 3 or loglik.shape[0] < 1:
        raise ValueError("loglik must contain at least one track")
    if state_x.shape != (loglik.shape[2],) or len(names) != loglik.shape[0]:
        raise ValueError("state grid or names do not match the likelihood")
    if loglik.shape[1] < 1:
        raise ValueError("at least one position is required")
    return loglik, state_x

def _normalize_track_prior(value: np.ndarray, n_track: int, n_state: int) -> np.ndarray:
    value = np.asarray(value, dtype=float)
    if value.ndim == 1:
        if value.shape != (n_state,):
            raise ValueError(f"prior must have shape ({n_state},) or ({n_track}, {n_state})")
        value = np.broadcast_to(value, (n_track, n_state)).copy()
    elif value.shape != (n_track, n_state):
        raise ValueError(f"prior must have shape ({n_state},) or ({n_track}, {n_state})")
    value -= logsumexp(value, axis=1, keepdims=True)
    return np.ascontiguousarray(value)


def _transition_model(x, log_prior, sd, forbid_same):
    m = len(x)
    if sd is not None:
        return (_SPARSE, *_sparse_transition(x, log_prior, sd, forbid_same))
    mode = _INDEPENDENT_NO_SELF if forbid_same else _INDEPENDENT
    if mode == _INDEPENDENT_NO_SELF:
        p = np.exp(log_prior)
        if m < 2 or np.any(p >= 1):
            raise ValueError("forbid_same_state requires at least two positive states")
        log_nu = log_prior + np.log1p(-p)
        log_nu -= logsumexp(log_nu)
    else:
        log_nu = log_prior.copy()
    empty_i = np.empty(0, dtype=np.int64)
    empty_f = np.empty(0)
    return mode, empty_i, empty_i, empty_f, empty_i, empty_i, empty_f, log_nu


def _sparse_transition(x, log_prior, sd, forbid_same):
    """Return row/column sparse forms and the stationary state mass."""
    m = len(x)
    if m == 1:
        if forbid_same:
            raise ValueError("forbid_same_state requires at least two states")
        ptr = np.array([0, 1], dtype=np.int64)
        idx = np.array([0], dtype=np.int64)
        zero = np.zeros(1)
        return ptr, idx, zero, ptr, idx, zero, zero

    log_t = np.broadcast_to(log_prior, (m, m)).copy()
    log_t -= (x[:, None] - x[None]) ** 2 / (2 * sd * sd)
    if forbid_same:
        np.fill_diagonal(log_t, -np.inf)
    log_t -= logsumexp(log_t, axis=1, keepdims=True)

    # Conversion to probabilities removes only transitions below float64's
    # representable range and gives a compact band for narrow attraction.
    t = np.exp(log_t)
    t /= t.sum(axis=1, keepdims=True)
    stationary = np.full(m, 1 / m)
    for _ in range(10000):
        new = stationary @ t
        if np.max(np.abs(new - stationary)) < 1e-14:
            stationary = new
            break
        stationary = new
    log_nu = np.log(stationary / stationary.sum())

    rows, cols = np.nonzero(t)
    values = np.log(t[rows, cols])
    row_ptr = np.zeros(m + 1, dtype=np.int64)
    np.add.at(row_ptr, rows + 1, 1)
    np.cumsum(row_ptr, out=row_ptr)

    order = np.lexsort((rows, cols))
    c_rows = rows[order].astype(np.int64)
    c_cols = cols[order]
    c_values = values[order]
    col_ptr = np.zeros(m + 1, dtype=np.int64)
    np.add.at(col_ptr, c_cols + 1, 1)
    np.cumsum(col_ptr, out=col_ptr)

    return (
        np.ascontiguousarray(row_ptr),
        np.ascontiguousarray(cols.astype(np.int64)),
        np.ascontiguousarray(values),
        np.ascontiguousarray(col_ptr),
        np.ascontiguousarray(c_rows),
        np.ascontiguousarray(c_values),
        np.ascontiguousarray(log_nu),
    )


def _length_terms(log_mass, log_tail, n):
    m, l0 = log_mass.shape
    survival = np.full((m, n + 1), -np.inf)
    equilibrium = np.full((m, n + 1), -np.inf)
    mean = np.empty(m)
    for state in range(m):
        explicit = np.exp(log_mass[state])
        r = np.exp(log_tail[state])
        for visible in range(1, n + 1):
            if visible <= l0:
                s = explicit[visible - 1:].sum() + explicit[-1] * r / (1 - r)
                e = np.sum((np.arange(visible, l0 + 1) - visible + 1) * explicit[visible - 1:])
                if r > 0:
                    e += explicit[-1] * (
                        (l0 - visible + 1) * r / (1 - r)
                        + r / (1 - r) ** 2
                    )
            else:
                s = explicit[-1] * r ** (visible - l0) / (1 - r) if r > 0 else 0.0
                e = explicit[-1] * r ** (visible - l0) / (1 - r) ** 2 if r > 0 else 0.0
            survival[state, visible] = np.log(s) if s > 0 else -np.inf
            equilibrium[state, visible] = np.log(e) if e > 0 else -np.inf
        base = np.sum(np.arange(1, l0 + 1) * explicit)
        if r > 0:
            base += explicit[-1] * (
                l0 * r / (1 - r) + r / (1 - r) ** 2
            )
        mean[state] = base
    return mean, survival, equilibrium


def _probabilities(log_weight):
    log_weight = np.asarray(log_weight, dtype=float)
    peak = np.max(log_weight)
    if not np.isfinite(peak):
        raise RuntimeError("no segmentation has positive probability")
    weight = np.exp(log_weight - peak)
    return weight / weight.sum()


def _duration_log_mass(log_q, log_tail, state, lengths):
    lengths = np.asarray(lengths)
    l0 = log_q.shape[1]
    out = np.full(lengths.shape, -np.inf)
    explicit = lengths <= l0
    out[explicit] = log_q[state, lengths[explicit] - 1]
    tail = ~explicit
    if np.isfinite(log_tail[state]):
        out[tail] = log_q[state, -1] + (lengths[tail] - l0) * log_tail[state]
    return out


def _transition_candidates(previous, prior, model, reverse=None):
    mode, row_ptr, col_idx, row_logp, *_ = model
    if mode == _SPARSE:
        sl = slice(row_ptr[previous], row_ptr[previous + 1])
        state, logp = col_idx[sl], row_logp[sl]
    else:
        state = np.arange(prior.size)
        if mode == _INDEPENDENT_NO_SELF:
            state = state[state != previous]
            logp = prior[state] - np.log1p(-np.exp(prior[previous]))
        else:
            logp = prior
    if reverse is not None:
        logp = logp + reverse[state]
    return state, _probabilities(logp)


def _path(starts, ends, states, state_x):
    states = np.asarray(states, dtype=np.int64)
    return SegmentationPath(
        np.asarray(starts, dtype=np.int64),
        np.asarray(ends, dtype=np.int64),
        states,
        state_x[states],
    )


def _posterior_paths(state_x, prior, log_q, log_tail, survival, equilibrium,
                     mean, model, prefix, reverse, outgoing, n_draws, rng):
    n, m = prefix.shape[0] - 1, state_x.size
    log_nu = model[-1]
    first = np.empty((m, n))
    for state in range(m):
        if n > 1:
            first[state, :-1] = (
                log_nu[state] + survival[state, 1:n]
                + prefix[1:n, state] + outgoing[1:n, state]
            )
        first[state, -1] = (
            log_nu[state] + equilibrium[state, n] + prefix[n, state]
        )
    first_state = _probabilities(logsumexp(first, axis=1))
    first_end = [_probabilities(row) for row in first]
    state_cache, end_cache = {}, {}
    paths = []
    for _ in range(n_draws):
        state = rng.choice(m, p=first_state)
        end = rng.choice(np.arange(1, n + 1), p=first_end[state])
        starts, ends, states = [0], [end], [state]
        while end < n:
            key = (end, state)
            if key not in state_cache:
                state_cache[key] = _transition_candidates(
                    state, prior, model, reverse[end]
                )
            candidates, probability = state_cache[key]
            state = rng.choice(candidates, p=probability)
            start = end
            key = (start, state)
            if key not in end_cache:
                possible = np.arange(start + 1, n + 1)
                weight = np.empty(n - start)
                if possible.size > 1:
                    lengths = possible[:-1] - start
                    weight[:-1] = (
                        _duration_log_mass(log_q, log_tail, state, lengths)
                        + prefix[possible[:-1], state] - prefix[start, state]
                        + outgoing[possible[:-1], state]
                    )
                weight[-1] = (
                    survival[state, n - start]
                    + prefix[n, state] - prefix[start, state]
                )
                end_cache[key] = (possible, _probabilities(weight))
            possible, probability = end_cache[key]
            end = rng.choice(possible, p=probability)
            starts.append(start); ends.append(end); states.append(state)
        paths.append(_path(starts, ends, states, state_x))
    return tuple(paths)


def _prior_paths(state_x, prior, log_q, log_tail, survival, equilibrium,
                 mean, model, n, n_draws, rng):
    m = state_x.size
    first_state = _probabilities(model[-1] + np.log(mean))
    first_end = []
    for state in range(m):
        weight = np.empty(n)
        if n > 1:
            weight[:-1] = survival[state, 1:n] - np.log(mean[state])
        weight[-1] = equilibrium[state, n] - np.log(mean[state])
        first_end.append(_probabilities(weight))
    transition_cache, length_cache = {}, {}
    paths = []
    for _ in range(n_draws):
        state = rng.choice(m, p=first_state)
        end = rng.choice(np.arange(1, n + 1), p=first_end[state])
        starts, ends, states = [0], [end], [state]
        while end < n:
            if state not in transition_cache:
                transition_cache[state] = _transition_candidates(state, prior, model)
            candidates, probability = transition_cache[state]
            state = rng.choice(candidates, p=probability)
            start = end
            key = (n - start, state)
            if key not in length_cache:
                remaining = n - start
                possible = np.arange(1, remaining + 1)
                weight = np.empty(remaining)
                if remaining > 1:
                    weight[:-1] = _duration_log_mass(
                        log_q, log_tail, state, possible[:-1]
                    )
                weight[-1] = survival[state, remaining]
                length_cache[key] = (possible, _probabilities(weight))
            lengths, probability = length_cache[key]
            end = start + rng.choice(lengths, p=probability)
            starts.append(start); ends.append(end); states.append(state)
        paths.append(_path(starts, ends, states, state_x))
    return tuple(paths)


@njit(cache=True, inline="always")
def _logadd(a, b):
    if a == -np.inf:
        return b
    if b == -np.inf:
        return a
    x = max(a, b)
    return x + np.log(np.exp(a - x) + np.exp(b - x))


@njit(cache=True, inline="always")
def _logsum(values):
    result = -np.inf
    for i in range(values.size):
        result = _logadd(result, values[i])
    return result


@njit(cache=True)
def _excluding_logsum(values, out, prefix, suffix):
    m = values.size
    prefix[0] = -np.inf
    for i in range(m):
        prefix[i + 1] = _logadd(prefix[i], values[i])
    suffix[m] = -np.inf
    for i in range(m - 1, -1, -1):
        suffix[i] = _logadd(values[i], suffix[i + 1])
    for i in range(m):
        out[i] = _logadd(prefix[i], suffix[i + 1])


@njit(cache=True)
def _forward_transition(values, out, mode, log_prior, log_denom,
                        col_ptr, row_idx, col_logp, scratch, prefix, suffix):
    m = values.size
    if mode == _INDEPENDENT:
        total = _logsum(values)
        for target in range(m):
            out[target] = total + log_prior[target]
    elif mode == _INDEPENDENT_NO_SELF:
        for source in range(m):
            scratch[source] = values[source] - log_denom[source]
        _excluding_logsum(scratch, out, prefix, suffix)
        for target in range(m):
            out[target] += log_prior[target]
    else:
        for target in range(m):
            value = -np.inf
            for j in range(col_ptr[target], col_ptr[target + 1]):
                value = _logadd(value, values[row_idx[j]] + col_logp[j])
            out[target] = value


@njit(cache=True)
def _reverse_transition(values, out, mode, log_prior, log_denom,
                        row_ptr, col_idx, row_logp, scratch, prefix, suffix):
    m = values.size
    if mode == _INDEPENDENT:
        for target in range(m):
            scratch[target] = log_prior[target] + values[target]
        total = _logsum(scratch)
        for source in range(m):
            out[source] = total
    elif mode == _INDEPENDENT_NO_SELF:
        for target in range(m):
            scratch[target] = log_prior[target] + values[target]
        _excluding_logsum(scratch, out, prefix, suffix)
        for source in range(m):
            out[source] -= log_denom[source]
    else:
        for source in range(m):
            value = -np.inf
            for j in range(row_ptr[source], row_ptr[source + 1]):
                value = _logadd(value, row_logp[j] + values[col_idx[j]])
            out[source] = value


@njit(cache=True, nogil=True)
def _outgoing_messages(reverse, log_prior, mode, row_ptr, col_idx, row_logp,
                       col_ptr, row_idx, col_logp):
    n, m = reverse.shape
    log_denom = np.zeros(m)
    if mode == _INDEPENDENT_NO_SELF:
        for state in range(m):
            log_denom[state] = np.log1p(-np.exp(log_prior[state]))
    outgoing = np.full((n, m), -np.inf)
    scratch = np.empty(m)
    work_prefix = np.empty(m + 1)
    work_suffix = np.empty(m + 1)
    for position in range(n):
        _reverse_transition(
            reverse[position], outgoing[position], mode, log_prior, log_denom,
            row_ptr, col_idx, row_logp, scratch, work_prefix, work_suffix,
        )
    return outgoing


@njit(cache=True)
def _boundary_transition(forward, reverse, mode, log_prior, log_denom,
                         row_ptr, col_idx, row_logp, scratch, excluded,
                         prefix, suffix):
    m = forward.size
    if mode == _INDEPENDENT:
        for target in range(m):
            scratch[target] = log_prior[target] + reverse[target]
        return _logsum(forward) + _logsum(scratch)
    if mode == _INDEPENDENT_NO_SELF:
        for target in range(m):
            scratch[target] = log_prior[target] + reverse[target]
        _excluding_logsum(scratch, excluded, prefix, suffix)
        value = -np.inf
        for source in range(m):
            value = _logadd(
                value, forward[source] - log_denom[source] + excluded[source]
            )
        return value
    value = -np.inf
    for source in range(m):
        for j in range(row_ptr[source], row_ptr[source + 1]):
            value = _logadd(
                value,
                forward[source] + row_logp[j] + reverse[col_idx[j]],
            )
    return value


@njit(cache=True, nogil=True)
def _prefix_sum(emission):
    n, m = emission.shape
    prefix = np.zeros((n + 1, m))
    for position in range(n):
        for state in range(m):
            prefix[position + 1, state] = (
                prefix[position, state] + emission[position, state]
            )
    return prefix


@njit(cache=True, nogil=True)
def _forward_messages(prefix, log_q, log_tail, log_survival, log_prior,
                      mode, row_ptr, col_idx, row_logp,
                      col_ptr, row_idx, col_logp, log_nu, mean_length):
    n, m = prefix.shape[0] - 1, prefix.shape[1]
    l0, log_mean = log_q.shape[1], np.log(mean_length)
    log_denom = np.empty(m)
    for state in range(m):
        log_denom[state] = (
            np.log1p(-np.exp(log_prior[state]))
            if mode == _INDEPENDENT_NO_SELF else 0.0
        )
    scratch = np.empty(m)
    work_prefix = np.empty(m + 1)
    work_suffix = np.empty(m + 1)
    forward = np.full((n + 1, m), -np.inf)
    incoming = np.full((n + 1, m), -np.inf)
    tail_acc = np.full(m, -np.inf)
    for end in range(1, n + 1):
        eligible = end - l0 - 1
        for state in range(m):
            if eligible >= 1 and log_tail[state] > -np.inf:
                tail_acc[state] = _logadd(
                    tail_acc[state],
                    incoming[eligible, state] - prefix[eligible, state]
                    - eligible * log_tail[state],
                )
            value = (
                log_nu[state] + log_survival[state, end]
                - log_mean + prefix[end, state]
            )
            for length in range(1, min(l0, end - 1) + 1):
                start = end - length
                if log_q[state, length - 1] > -np.inf:
                    value = _logadd(
                        value,
                        incoming[start, state] + log_q[state, length - 1]
                        + prefix[end, state] - prefix[start, state],
                    )
            if tail_acc[state] > -np.inf:
                value = _logadd(
                    value,
                    prefix[end, state] + end * log_tail[state]
                    + log_q[state, -1] - l0 * log_tail[state]
                    + tail_acc[state],
                )
            forward[end, state] = value
        _forward_transition(
            forward[end], incoming[end], mode, log_prior, log_denom,
            col_ptr, row_idx, col_logp, scratch, work_prefix, work_suffix,
        )
    return forward, incoming


@njit(cache=True, nogil=True)
def _reverse_messages(prefix, log_q, log_tail, log_survival, log_prior,
                      mode, row_ptr, col_idx, row_logp,
                      col_ptr, row_idx, col_logp):
    n, m = prefix.shape[0] - 1, prefix.shape[1]
    l0 = log_q.shape[1]
    log_denom = np.empty(m)
    for state in range(m):
        log_denom[state] = (
            np.log1p(-np.exp(log_prior[state]))
            if mode == _INDEPENDENT_NO_SELF else 0.0
        )
    scratch = np.empty(m)
    work_prefix = np.empty(m + 1)
    work_suffix = np.empty(m + 1)
    reverse = np.full((n + 1, m), -np.inf)
    outgoing = np.full((n + 1, m), -np.inf)
    tail_acc = np.full(m, -np.inf)
    for start in range(n - 1, -1, -1):
        eligible = start + l0 + 1
        for state in range(m):
            if eligible < n and log_tail[state] > -np.inf:
                tail_acc[state] = _logadd(
                    tail_acc[state],
                    prefix[eligible, state] + outgoing[eligible, state]
                    + eligible * log_tail[state],
                )
            remaining = n - start
            value = (
                log_survival[state, remaining]
                + prefix[n, state] - prefix[start, state]
            )
            for length in range(1, min(l0, remaining - 1) + 1):
                end = start + length
                if log_q[state, length - 1] > -np.inf:
                    value = _logadd(
                        value,
                        log_q[state, length - 1]
                        + prefix[end, state] - prefix[start, state]
                        + outgoing[end, state],
                    )
            if tail_acc[state] > -np.inf:
                value = _logadd(
                    value,
                    log_q[state, -1] - l0 * log_tail[state]
                    - prefix[start, state] - start * log_tail[state]
                    + tail_acc[state],
                )
            reverse[start, state] = value
        _reverse_transition(
            reverse[start], outgoing[start], mode, log_prior, log_denom,
            row_ptr, col_idx, row_logp, scratch, work_prefix, work_suffix,
        )
    return reverse, outgoing


@njit(cache=True, nogil=True)
def _start_log_partition(prefix, outgoing, log_survival, log_equilibrium,
                         log_nu, mean_length):
    n, m = prefix.shape[0] - 1, prefix.shape[1]
    value = -np.inf
    log_mean = np.log(mean_length)
    for state in range(m):
        value = _logadd(
            value,
            log_nu[state] + log_equilibrium[state, n]
            - log_mean + prefix[n, state],
        )
        for end in range(1, n):
            value = _logadd(
                value,
                log_nu[state] + log_survival[state, end]
                - log_mean + prefix[end, state] + outgoing[end, state],
            )
    return value


@njit(cache=True, nogil=True)
def _partition(emission, log_q, log_tail, log_survival, log_equilibrium,
               log_prior, mode, row_ptr, col_idx, row_logp,
               col_ptr, row_idx, col_logp, log_nu, mean_length):
    n, m = emission.shape
    l0, log_mean = log_q.shape[1], np.log(mean_length)
    prefix = _prefix_sum(emission)
    forward, incoming = _forward_messages(
        prefix, log_q, log_tail, log_survival, log_prior, mode,
        row_ptr, col_idx, row_logp, col_ptr, row_idx, col_logp,
        log_nu, mean_length,
    )
    reverse, outgoing = _reverse_messages(
        prefix, log_q, log_tail, log_survival, log_prior, mode,
        row_ptr, col_idx, row_logp, col_ptr, row_idx, col_logp,
    )
    log_z = _start_log_partition(
        prefix, outgoing, log_survival, log_equilibrium, log_nu, mean_length
    )
    log_denom = np.empty(m)
    for state in range(m):
        log_denom[state] = (
            np.log1p(-np.exp(log_prior[state]))
            if mode == _INDEPENDENT_NO_SELF else 0.0
        )
    scratch = np.empty(m)
    excluded = np.empty(m)
    work_prefix = np.empty(m + 1)
    work_suffix = np.empty(m + 1)

    scale = np.full(m, -np.inf)
    for state in range(m):
        scale[state] = (
            log_nu[state] + log_equilibrium[state, n]
            - log_mean + prefix[n, state] - log_z
        )
        for end in range(1, n):
            scale[state] = max(
                scale[state],
                log_nu[state] + log_survival[state, end]
                - log_mean + prefix[end, state]
                + outgoing[end, state] - log_z,
            )
        for start in range(1, n):
            scale[state] = max(
                scale[state],
                incoming[start, state] + log_survival[state, n - start]
                + prefix[n, state] - prefix[start, state] - log_z,
            )
        for start in range(1, n - 1):
            max_length = min(l0, n - start - 1)
            for length in range(1, max_length + 1):
                end = start + length
                if log_q[state, length - 1] > -np.inf:
                    scale[state] = max(
                        scale[state],
                        incoming[start, state] + log_q[state, length - 1]
                        + prefix[end, state] - prefix[start, state]
                        + outgoing[end, state] - log_z,
                    )

        if log_tail[state] > -np.inf:
            best_start = -np.inf
            tail_const = log_q[state, l0 - 1] - l0 * log_tail[state]
            for end in range(1, n):
                new_start = end - l0 - 1
                if new_start >= 1:
                    best_start = max(
                        best_start,
                        incoming[new_start, state]
                        - prefix[new_start, state]
                        - new_start * log_tail[state],
                    )
                if best_start > -np.inf:
                    scale[state] = max(
                        scale[state],
                        tail_const + best_start + prefix[end, state]
                        + outgoing[end, state] + end * log_tail[state]
                        - log_z,
                    )

    difference = np.zeros((n + 1, m))
    suffix_tail = np.full(n + 1, -np.inf)
    for state in range(m):
        weight = np.exp(
            log_nu[state] + log_equilibrium[state, n]
            - log_mean + prefix[n, state] - log_z - scale[state]
        )
        difference[0, state] += weight
        difference[n, state] -= weight

        for end in range(1, n):
            weight = np.exp(
                log_nu[state] + log_survival[state, end]
                - log_mean + prefix[end, state] + outgoing[end, state]
                - log_z - scale[state]
            )
            difference[0, state] += weight
            difference[end, state] -= weight

        for start in range(1, n):
            weight = np.exp(
                incoming[start, state] + log_survival[state, n - start]
                + prefix[n, state] - prefix[start, state]
                - log_z - scale[state]
            )
            difference[start, state] += weight
            difference[n, state] -= weight

        for start in range(1, n - 1):
            max_length = min(l0, n - start - 1)
            for length in range(1, max_length + 1):
                end = start + length
                if log_q[state, length - 1] > -np.inf:
                    value = (
                        incoming[start, state] + log_q[state, length - 1]
                        + prefix[end, state] - prefix[start, state]
                        + outgoing[end, state] - log_z
                    )
                    weight = np.exp(value - scale[state])
                    difference[start, state] += weight
                    difference[end, state] -= weight

        if log_tail[state] > -np.inf:
            tail_const = log_q[state, l0 - 1] - l0 * log_tail[state] - log_z
            suffix_tail[:] = -np.inf
            for end in range(n - 1, 0, -1):
                value = (
                    prefix[end, state] + outgoing[end, state]
                    + end * log_tail[state]
                )
                suffix_tail[end] = _logadd(value, suffix_tail[end + 1])

            for start in range(1, n - l0 - 1):
                first_end = start + l0 + 1
                value = (
                    tail_const + incoming[start, state]
                    - prefix[start, state] - start * log_tail[state]
                    + suffix_tail[first_end]
                )
                difference[start, state] += np.exp(value - scale[state])

            prefix_tail = -np.inf
            for end in range(1, n):
                new_start = end - l0 - 1
                if new_start >= 1:
                    value = (
                        incoming[new_start, state]
                        - prefix[new_start, state]
                        - new_start * log_tail[state]
                    )
                    prefix_tail = _logadd(prefix_tail, value)
                if prefix_tail > -np.inf:
                    value = (
                        tail_const + prefix[end, state]
                        + outgoing[end, state] + end * log_tail[state]
                        + prefix_tail
                    )
                    difference[end, state] -= np.exp(value - scale[state])

    log_posterior = np.full((n, m), -np.inf)
    running = np.zeros(m)
    for p in range(n):
        normalizer = -np.inf
        for state in range(m):
            running[state] += difference[p, state]
            if running[state] > 0:
                log_posterior[p, state] = np.log(running[state]) + scale[state]
            normalizer = _logadd(normalizer, log_posterior[p, state])
        for state in range(m):
            log_posterior[p, state] -= normalizer

    boundary = np.empty(max(n - 1, 0))
    for b in range(1, n):
        value = _boundary_transition(
            forward[b], reverse[b], mode, log_prior, log_denom,
            row_ptr, col_idx, row_logp, scratch, excluded,
            work_prefix, work_suffix,
        )
        boundary[b - 1] = min(1.0, np.exp(value - log_z))
    return log_z, log_posterior, boundary, prefix, reverse
