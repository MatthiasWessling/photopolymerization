"""Conservative 1-D diffusion ``d/dx (D dC/dx)`` with Neumann or Dirichlet."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray


def divergence_diffusion_1d(
    C: NDArray[np.floating],
    D: NDArray[np.floating],
    dx: float,
    *,
    left: str = "neumann",
    right: str = "neumann",
    left_value: float = 0.0,
    right_value: float = 0.0,
) -> NDArray[np.floating]:
    """Finite-volume Laplacian with harmonic face diffusivities.

    ``left``/``right`` are ``'neumann'`` (zero flux) or ``'dirichlet'``.
    Dirichlet nodes are held by returning zero time derivative there; the
    caller must set the nodal value to the boundary datum.
    """
    C = np.asarray(C, dtype=float)
    D = np.asarray(D, dtype=float)
    n = C.size
    out = np.zeros(n, dtype=float)
    if n < 2:
        return out

    D_face = 2.0 * D[:-1] * D[1:] / (D[:-1] + D[1:] + 1e-30)
    flux = D_face * (C[1:] - C[:-1]) / dx

    flux_left = 0.0
    flux_right = 0.0
    if left == "dirichlet":
        D_g = D[0]
        C_ghost = 2.0 * left_value - C[0]
        flux_left = D_g * (C[0] - C_ghost) / dx
    if right == "dirichlet":
        D_g = D[-1]
        C_ghost = 2.0 * right_value - C[-1]
        flux_right = D_g * (C_ghost - C[-1]) / dx

    out[0] = (flux[0] - flux_left) / dx
    out[-1] = (flux_right - flux[-1]) / dx
    if n > 2:
        out[1:-1] = (flux[1:] - flux[:-1]) / dx
    if left == "dirichlet":
        out[0] = 0.0
    if right == "dirichlet":
        out[-1] = 0.0
    return out


def fourier_explicit_limit(D_max: float, dx: float) -> float:
    """1-D explicit Fourier limit ``dt <= 0.5 dx^2 / D`` (documented, not used)."""
    if D_max <= 0.0:
        return float("inf")
    return 0.5 * dx * dx / D_max
