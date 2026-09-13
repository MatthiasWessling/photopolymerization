"""Irradiance protocols ``I(t)`` in W/m^2."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass


IntensityFn = Callable[[float], float]


@dataclass(frozen=True)
class ConstantIntensity:
    I0: float

    def __call__(self, t: float) -> float:
        return float(self.I0)


@dataclass(frozen=True)
class SquarePulse:
    I0: float
    t_on: float
    t_off: float

    def __call__(self, t: float) -> float:
        if self.t_on <= t < self.t_off:
            return float(self.I0)
        return 0.0
