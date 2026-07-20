"""Small serialization helpers for likelihood and result objects."""

from __future__ import annotations

from typing import ClassVar

import numpy as np


class Serializable:
    """Mixin for compact dict/NPZ serialization.

    Subclasses either set ``save_attrs`` or override ``to_dict``/``from_dict``.
    Values must be arrays, scalars, or tuples coercible by ``np.savez``.
    """

    save_attrs: ClassVar[tuple[str, ...]] = ()
    optional_save_attrs: ClassVar[tuple[str, ...]] = ()

    def to_dict(self) -> dict[str, object]:
        out = {name: getattr(self, name) for name in self.save_attrs}
        for name in self.optional_save_attrs:
            value = getattr(self, name)
            if value is not None:
                out[name] = value
        return out

    @classmethod
    def from_dict(cls, data: dict[str, object]):
        return cls(**{name: data[name] for name in cls.save_attrs})

    def to_npz(self, path) -> None:
        np.savez_compressed(path, **self.to_dict())

    @classmethod
    def from_npz(cls, path):
        with np.load(path, allow_pickle=False) as x:
            return cls.from_dict({name: x[name] for name in x.files})


def tuple_str(value) -> tuple[str, ...]:
    return tuple(str(x) for x in np.asarray(value).tolist())


def optional_float(value) -> float | None:
    if value is None:
        return None
    value = float(np.asarray(value).item())
    return None if np.isnan(value) else value
