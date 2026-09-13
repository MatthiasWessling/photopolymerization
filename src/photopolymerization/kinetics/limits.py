"""Analytic limits for Stage 2 (constant-rate fork).

Slutzky, Stone & Nunes, Soft Matter 15, 9553 (2019), Eq. 3:
``C_O / C_{O,0} = (C_M / C_{M,0})^{k_O / k_p}`` when ``k_p`` and ``k_O``
are constant.

Quasi-steady radicals (oxygen-free, constant ``k_t``), cf. Zhu et al.,
J. Mech. Phys. Solids 142, 104041 (2020), Eq. 15 analogue:
``C_R ≈ sqrt(R_i / (2 k_t))``.
"""

from __future__ import annotations

import numpy as np

from photopolymerization.kinetics.closures import kt
from photopolymerization.kinetics.solve import Trajectory
from photopolymerization.parameters import KineticParams


def slutzky_oxygen_monomer_residual(
    traj: Trajectory,
    params: KineticParams,
    *,
    min_oxygen_frac: float = 0.05,
) -> float:
    """Max relative residual of Slutzky Eq. 3 while oxygen is still present.

    With ``k_O / k_p ~ 10^3``, oxygen is exhausted at tiny conversion, so the
    identity is only tested on the induction window (not after both sides
    underflow to zero).
    """
    return slutzky_eq3_residual(
        traj.y[2],
        traj.y[3],
        params.C_O0,
        params.C_M0,
        params.k_oxygen,
        params.kp0,
        min_oxygen_frac=min_oxygen_frac,
    )


def slutzky_eq3_residual(
    C_O: np.ndarray,
    C_M: np.ndarray,
    C_O0: float,
    C_M0: float,
    k_oxygen: float,
    kp: float,
    *,
    min_oxygen_frac: float = 0.05,
) -> float:
    """Max relative residual of ``C_O/C_O0 = (C_M/C_M0)^{k_O/k_p}``."""
    frac_M = np.clip(np.asarray(C_M, dtype=float) / C_M0, 1e-30, 1.0)
    frac_O = np.asarray(C_O, dtype=float) / C_O0
    mask = (frac_O > min_oxygen_frac) & (C_O0 > 0.0)
    if not np.any(mask):
        return float("nan")
    exponent = k_oxygen / kp
    predicted = np.exp(exponent * np.log(frac_M[mask]))
    denom = np.maximum(frac_O[mask], 1e-12)
    return float(np.max(np.abs(predicted - frac_O[mask]) / denom))


def qss_radical(
    C_I: np.ndarray | float,
    I: np.ndarray | float,
    params: KineticParams,
) -> np.ndarray:
    """Oxygen-free QSS ``C_R = sqrt(m beta I C_I / (2 k_t(p=0)))``."""
    kt0 = float(kt(0.0, params))
    Ri = params.m * params.beta * np.asarray(I, dtype=float) * np.asarray(
        C_I, dtype=float
    )
    return np.sqrt(np.maximum(Ri, 0.0) / (2.0 * kt0))


def qss_residual(
    traj: Trajectory,
    params: KineticParams,
    *,
    t_skip: float,
) -> float:
    """Relative QSS residual after the fast radical transient."""
    mask = traj.t >= t_skip
    if not np.any(mask):
        return float("nan")
    C_R = traj.y[1, mask]
    pred = qss_radical(traj.y[0, mask], traj.I[mask], params)
    scale = np.maximum(pred, 1e-18)
    return float(np.max(np.abs(C_R - pred) / scale))
