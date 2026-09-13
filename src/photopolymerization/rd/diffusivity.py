"""Conversion-dependent diffusivity, Montgomery Eq. 28."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from photopolymerization.parameters import DiffusionParams


def harmonic_D(
    p: NDArray[np.floating] | float,
    D_liquid: float,
    D_solid: float,
) -> NDArray[np.floating] | float:
    """``1/D = p/D_solid + (1-p)/D_liquid``."""
    p = np.clip(np.asarray(p, dtype=float), 0.0, 1.0)
    return 1.0 / (p / D_solid + (1.0 - p) / D_liquid)


def species_diffusivities(
    p: NDArray[np.floating],
    diff: DiffusionParams,
    *,
    enabled: bool = True,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """``(D_I, D_R, D_O, D_M)`` at each node."""
    p = np.asarray(p, dtype=float)
    if not enabled:
        z = np.zeros_like(p)
        return z, z.copy(), z.copy(), z.copy()
    return (
        np.asarray(harmonic_D(p, diff.D_I_liquid, diff.D_I_solid), dtype=float),
        np.asarray(harmonic_D(p, diff.D_R_liquid, diff.D_R_solid), dtype=float),
        np.asarray(harmonic_D(p, diff.D_O_liquid, diff.D_O_solid), dtype=float),
        np.asarray(harmonic_D(p, diff.D_M_liquid, diff.D_M_solid), dtype=float),
    )
