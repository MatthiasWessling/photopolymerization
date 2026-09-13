"""Steady plug flow: ``U dy/dx = R(y, kd(x))``, Slutzky Eqs. 1a--d."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

from photopolymerization.flow.params import SlutzkyParams
from photopolymerization.kinetics.rhs import CI, CM, CO, CR

PI, R, O2, M = CI, CR, CO, CM


def kd_of_x(x: float, params: SlutzkyParams) -> float:
    """Photolysis only inside the UV window ``[0, L]``."""
    if 0.0 <= x <= params.L:
        return float(params.kd)
    return 0.0


def reaction_slutzky(
    y: NDArray[np.floating],
    kd: float,
    params: SlutzkyParams,
) -> NDArray[np.floating]:
    """Time-like rates (Lagrangian). State ``(Pi, R, y_O2, M)``."""
    Pi, rad, yox, mon = np.asarray(y, dtype=float)
    dPi = -kd * Pi
    dR = kd * Pi - params.kt * rad * rad - params.ko * yox * rad
    dO = -params.ko * yox * rad
    dM = -params.kp * mon * rad
    return np.array([dPi, dR, dO, dM], dtype=float)


def rhs_x(
    x: float,
    y: NDArray[np.floating],
    U: float,
    params: SlutzkyParams,
) -> NDArray[np.floating]:
    if U <= 0.0:
        raise ValueError("plug speed U must be positive")
    return reaction_slutzky(y, kd_of_x(x, params), params) / U


@dataclass(frozen=True)
class PlugTrajectory:
    x: NDArray[np.floating]
    y: NDArray[np.floating]
    kd: NDArray[np.floating]
    U: float
    message: str
    success: bool

    @property
    def M_over_M0(self) -> np.ndarray:
        return self.y[M] / self.y[M, 0]

    @property
    def p_double_bond(self) -> np.ndarray:
        """Slutzky conversion ``1 - M/M0`` (not Zhu P_gel)."""
        return 1.0 - self.M_over_M0


def initial_state(params: SlutzkyParams) -> np.ndarray:
    return np.array([params.Pi0, 0.0, params.y0, params.M0], dtype=float)


def integrate_plug(
    params: SlutzkyParams,
    U: float,
    *,
    x_end: float | None = None,
    n_eval: int = 400,
    rtol: float = 1e-7,
) -> PlugTrajectory:
    """Integrate from the window inlet ``x=0`` to ``x_end`` (default ``2L``)."""
    x_end = params.L * 2.0 if x_end is None else float(x_end)
    t_eval = np.linspace(0.0, x_end, n_eval)
    sol = solve_ivp(
        lambda x, y: rhs_x(x, y, U, params),
        (0.0, x_end),
        initial_state(params),
        method="BDF",
        rtol=rtol,
        atol=np.array([1e-10, 1e-18, 1e-12, 1e-8], dtype=float),
        t_eval=t_eval,
    )
    kd = np.array([kd_of_x(float(xi), params) for xi in sol.t], dtype=float)
    return PlugTrajectory(
        x=sol.t,
        y=sol.y,
        kd=kd,
        U=float(U),
        message=sol.message,
        success=bool(sol.success),
    )


def gelled(traj: PlugTrajectory, params: SlutzkyParams) -> bool:
    """Operational gel: ``min(M/M0) <= 0.98`` (Slutzky/Dendukuri, not Macosko)."""
    return float(np.min(traj.M_over_M0)) <= params.gel_M_over_M0


def critical_speed(
    params: SlutzkyParams,
    *,
    U_lo: float = 1e-5,
    U_hi: float = 1.0,
    n_bisect: int = 28,
) -> float:
    """Largest U (by bisection) at which the 2% double-bond predicate holds."""
    lo, hi = float(U_lo), float(U_hi)
    if not gelled(integrate_plug(params, lo), params):
        return 0.0
    if gelled(integrate_plug(params, hi), params):
        return hi
    for _ in range(n_bisect):
        mid = 0.5 * (lo + hi)
        if gelled(integrate_plug(params, mid), params):
            lo = mid
        else:
            hi = mid
    return float(lo)


def kinematic_fiber_length(U: float, t_uv: float) -> float:
    """Uniform-fiber length ``U t_uv`` (pulse duration; not a conversion field)."""
    return float(U * t_uv)


def residence_time(params: SlutzkyParams, U: float) -> float:
    return params.L / U


def tau_exp(t_uv: float, params: SlutzkyParams, U: float) -> float:
    """``t_uv / t_L = U t_uv / L``."""
    return float(U * t_uv / params.L)
