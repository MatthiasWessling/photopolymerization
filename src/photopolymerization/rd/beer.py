"""Beer--Lambert attenuation, Montgomery Eqs. 4--6."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import cumulative_trapezoid

from photopolymerization.parameters import OpticalParams


def absorption(
    C_I: NDArray[np.floating] | float,
    p: NDArray[np.floating] | float,
    optics: OpticalParams,
) -> NDArray[np.floating] | float:
    """Total absorption ``A`` (1/m)."""
    C_I = np.asarray(C_I, dtype=float)
    p = np.asarray(p, dtype=float)
    return (
        optics.alpha_I * C_I
        + p * optics.A_polymer
        + (1.0 - p) * optics.A_monomer
        + optics.w_absorber * optics.A_absorber
    )


def irradiance_along_z(
    z: NDArray[np.floating],
    A: NDArray[np.floating],
    I0: float,
) -> NDArray[np.floating]:
    """``I(z) = I0 exp(-int_0^z A dz')`` on a 1-D node grid with ``z[0]`` illuminated."""
    z = np.asarray(z, dtype=float)
    A = np.asarray(A, dtype=float)
    if z.size == 1:
        return np.array([I0], dtype=float)
    tau = cumulative_trapezoid(A, z, initial=0.0)
    return I0 * np.exp(-tau)


def optical_depth(A_char: float, H: float) -> float:
    """``beta_opt = A H``."""
    return float(A_char * H)
