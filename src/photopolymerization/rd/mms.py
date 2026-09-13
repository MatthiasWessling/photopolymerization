"""Manufactured 1-D diffusion (no reaction)."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

from photopolymerization.rd.operators import divergence_diffusion_1d


def cosine_mode_exact(
    x: NDArray[np.floating],
    t: float,
    D: float,
    length: float,
) -> NDArray[np.floating]:
    """Neumann eigenmode ``1 + exp(-D (pi/L)^2 t) cos(pi x / L)``."""
    k = np.pi / length
    return 1.0 + np.exp(-D * k * k * t) * np.cos(k * x)


def integrate_cosine_mms(
    x: NDArray[np.floating],
    D: float,
    t_end: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Diffuse the t=0 cosine with constant D and Neumann ends."""
    dx = float(x[1] - x[0])
    y0 = cosine_mode_exact(x, 0.0, D, float(x[-1] - x[0]))

    def rhs(_t: float, y: np.ndarray) -> np.ndarray:
        return divergence_diffusion_1d(y, np.full_like(y, D), dx)

    sol = solve_ivp(rhs, (0.0, t_end), y0, method="BDF", rtol=1e-8, atol=1e-10)
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol.y[:, -1], cosine_mode_exact(x, t_end, D, float(x[-1] - x[0]))
