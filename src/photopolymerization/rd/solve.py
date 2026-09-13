"""BDF method of lines for Stage 4 1-D fields."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

from photopolymerization.kinetics.closures import conversion
from photopolymerization.parameters import DiffusionParams, KineticParams, OpticalParams
from photopolymerization.rd.beer import absorption, irradiance_along_z
from photopolymerization.rd.mol import Grid1D, initial_field, mol_rhs_1d, unpack


@dataclass(frozen=True)
class SpatialTrajectory:
    t: NDArray[np.floating]
    x: NDArray[np.floating]
    C_I: NDArray[np.floating]
    C_R: NDArray[np.floating]
    C_O: NDArray[np.floating]
    C_M: NDArray[np.floating]
    p: NDArray[np.floating]
    I: NDArray[np.floating]
    message: str
    success: bool


def integrate_1d(
    grid: Grid1D,
    params: KineticParams,
    optics: OpticalParams,
    diff: DiffusionParams,
    I0_of_x: NDArray[np.floating] | float,
    t_end: float,
    *,
    along: str = "z",
    oxygen_bc: str = "sealed",
    diffuse: bool = True,
    rtol: float = 1.0e-6,
    n_eval: int = 25,
    max_step: float | None = 0.25,
) -> SpatialTrajectory:
    n = grid.n
    atol = np.tile(
        np.array([1e-12, 1e-18, 1e-12, 1e-8], dtype=float),
        n,
    )
    t_eval = np.linspace(0.0, t_end, n_eval)
    sol = solve_ivp(
        lambda t, y: mol_rhs_1d(
            t,
            y,
            grid,
            params,
            optics,
            diff,
            I0_of_x,
            along=along,
            oxygen_bc=oxygen_bc,
            diffuse=diffuse,
        ),
        (0.0, t_end),
        initial_field(params, n),
        method="BDF",
        rtol=rtol,
        atol=atol,
        t_eval=t_eval,
        max_step=np.inf if max_step is None else max_step,
    )
    nt = sol.t.size
    C_I = np.empty((nt, n))
    C_R = np.empty((nt, n))
    C_O = np.empty((nt, n))
    C_M = np.empty((nt, n))
    I_hist = np.empty((nt, n))
    for k, yk in enumerate(sol.y.T):
        ci, cr, co, cm = unpack(yk, n)
        if oxygen_bc == "dirichlet":
            co = co.copy()
            co[0] = params.C_O0
        C_I[k], C_R[k], C_O[k], C_M[k] = ci, cr, co, cm
        p = np.asarray(conversion(cm, params.C_M0), dtype=float)
        I0 = np.broadcast_to(np.asarray(I0_of_x, dtype=float), (n,))
        if along == "z":
            A = np.asarray(absorption(ci, p, optics), dtype=float)
            I_hist[k] = irradiance_along_z(grid.x, A, float(I0[0]))
        else:
            I_hist[k] = I0
    p_all = np.asarray(conversion(C_M, params.C_M0), dtype=float)
    return SpatialTrajectory(
        t=sol.t,
        x=grid.x,
        C_I=C_I,
        C_R=C_R,
        C_O=C_O,
        C_M=C_M,
        p=p_all,
        I=I_hist,
        message=sol.message,
        success=bool(sol.success),
    )
