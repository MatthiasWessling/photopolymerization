"""Stiff integration of the Stage 1 local ODE."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp

from photopolymerization.kinetics.closures import conversion, kp
from photopolymerization.kinetics.protocols import IntensityFn
from photopolymerization.kinetics.rhs import CM, initial_state, rhs
from photopolymerization.parameters import KineticParams


@dataclass(frozen=True)
class Trajectory:
    t: np.ndarray
    y: np.ndarray
    I: np.ndarray
    p: np.ndarray
    Rp: np.ndarray
    message: str
    success: bool


def integrate(
    params: KineticParams,
    intensity: IntensityFn,
    t_end: float,
    *,
    constant_rate: bool = False,
    rtol: float = 1.0e-7,
    atol: tuple[float, float, float, float] = (1e-12, 1e-18, 1e-12, 1e-8),
    max_step: float | None = None,
    n_eval: int = 400,
) -> Trajectory:
    """Integrate on ``[0, t_end]`` with BDF (radicals are stiff)."""
    t_eval = np.linspace(0.0, t_end, n_eval)
    sol = solve_ivp(
        lambda t, y: rhs(t, y, params, intensity, constant_rate=constant_rate),
        (0.0, t_end),
        initial_state(params),
        method="BDF",
        rtol=rtol,
        atol=np.asarray(atol, dtype=float),
        t_eval=t_eval,
        max_step=np.inf if max_step is None else max_step,
        dense_output=False,
    )
    y = sol.y
    I = np.array([intensity(float(ti)) for ti in sol.t], dtype=float)
    p = np.asarray(conversion(y[CM], params.C_M0), dtype=float)
    if constant_rate:
        kp_val = np.full_like(p, params.kp0)
    else:
        kp_val = np.asarray(kp(p, params), dtype=float)
    Rp = kp_val * y[1] * y[CM]
    return Trajectory(
        t=sol.t,
        y=y,
        I=I,
        p=p,
        Rp=Rp,
        message=sol.message,
        success=bool(sol.success),
    )
