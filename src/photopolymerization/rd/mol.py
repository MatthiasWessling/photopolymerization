"""Packed method-of-lines RHS for four species on a 1-D grid."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from photopolymerization.kinetics.closures import conversion
from photopolymerization.kinetics.rhs import N_STATE, reaction_rates
from photopolymerization.parameters import DiffusionParams, KineticParams, OpticalParams
from photopolymerization.rd.beer import absorption, irradiance_along_z
from photopolymerization.rd.diffusivity import species_diffusivities
from photopolymerization.rd.operators import divergence_diffusion_1d


@dataclass(frozen=True)
class Grid1D:
    x: NDArray[np.floating]
    dx: float

    @property
    def n(self) -> int:
        return int(self.x.size)


def make_grid(length: float, n: int) -> Grid1D:
    x = np.linspace(0.0, length, n)
    dx = float(x[1] - x[0]) if n > 1 else length
    return Grid1D(x=x, dx=dx)


def unpack(y: NDArray[np.floating], n: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    y = np.asarray(y, dtype=float)
    return y[0:n], y[n : 2 * n], y[2 * n : 3 * n], y[3 * n : 4 * n]


def pack(
    C_I: NDArray[np.floating],
    C_R: NDArray[np.floating],
    C_O: NDArray[np.floating],
    C_M: NDArray[np.floating],
) -> np.ndarray:
    return np.concatenate([C_I, C_R, C_O, C_M])


def initial_field(params: KineticParams, n: int) -> np.ndarray:
    return pack(
        np.full(n, params.C_I0),
        np.full(n, params.C_R0),
        np.full(n, params.C_O0),
        np.full(n, params.C_M0),
    )


def mol_rhs_1d(
    t: float,
    y: NDArray[np.floating],
    grid: Grid1D,
    params: KineticParams,
    optics: OpticalParams,
    diff: DiffusionParams,
    I0_of_x: NDArray[np.floating] | float,
    *,
    along: str = "z",
    oxygen_bc: str = "sealed",
    diffuse: bool = True,
) -> np.ndarray:
    """Reaction + diffusion. ``along='z'`` applies Beer--Lambert from ``x[0]``.

    ``along='x'`` treats ``I0_of_x`` as the in-plane irradiance (thin-film,
    no depth attenuation). ``oxygen_bc`` is ``'sealed'`` (Neumann) or
    ``'dirichlet'`` (fixed ``C_O`` at the illuminated/left node).
    """
    n = grid.n
    C_I, C_R, C_O, C_M = unpack(y, n)
    if oxygen_bc == "dirichlet":
        C_O = C_O.copy()
        C_O[0] = params.C_O0
    p = np.asarray(conversion(C_M, params.C_M0), dtype=float)
    I0 = np.broadcast_to(np.asarray(I0_of_x, dtype=float), (n,))
    if along == "z":
        A = np.asarray(absorption(C_I, p, optics), dtype=float)
        I = irradiance_along_z(grid.x, A, float(I0[0]))
    elif along == "x":
        I = I0
    else:
        raise ValueError(f"unknown axis {along}")
    dCI_r, dCR_r, dCO_r, dCM_r = reaction_rates(C_I, C_R, C_O, C_M, I, params)
    D_I, D_R, D_O, D_M = species_diffusivities(p, diff, enabled=diffuse)
    o2_left = "dirichlet" if oxygen_bc == "dirichlet" else "neumann"
    dCI = dCI_r + divergence_diffusion_1d(C_I, D_I, grid.dx)
    dCR = dCR_r + divergence_diffusion_1d(C_R, D_R, grid.dx)
    dCO = dCO_r + divergence_diffusion_1d(
        C_O, D_O, grid.dx, left=o2_left, left_value=params.C_O0
    )
    dCM = dCM_r + divergence_diffusion_1d(C_M, D_M, grid.dx)
    if oxygen_bc == "dirichlet":
        dCO[0] = 0.0
    return pack(dCI, dCR, dCO, dCM)


N_STATE_SPATIAL = N_STATE
