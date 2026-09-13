"""Local four-species photopolymerization RHS.

Montgomery et al., Extreme Mech. Lett. 53, 101714 (2022), Eqs. 16--19
with spatial diffusion omitted (Stage 1).

State ``y = (C_I, C_R, C_O, C_M)`` in mol/m^3.
"""

from __future__ import annotations

import numpy as np

from photopolymerization.kinetics.closures import conversion, kp, kt
from photopolymerization.kinetics.protocols import IntensityFn
from photopolymerization.parameters import KineticParams

N_STATE = 4
CI, CR, CO, CM = 0, 1, 2, 3


def reaction_rates(
    C_I: np.ndarray,
    C_R: np.ndarray,
    C_O: np.ndarray,
    C_M: np.ndarray,
    I: np.ndarray,
    params: KineticParams,
    *,
    constant_rate: bool = False,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Montgomery Eqs. 16--19 reaction terms (no diffusion).

    Accepts scalars or arrays so the Stage 4 method of lines can call
    the same chemistry as Stage 1.
    """
    C_I = np.asarray(C_I, dtype=float)
    C_R = np.asarray(C_R, dtype=float)
    C_O = np.asarray(C_O, dtype=float)
    C_M = np.asarray(C_M, dtype=float)
    I = np.asarray(I, dtype=float)
    p = conversion(C_M, params.C_M0)
    if constant_rate:
        kp_val = np.full_like(np.asarray(p, dtype=float), params.kp0)
        kt_val = np.full_like(kp_val, float(kt(0.0, params)))
    else:
        kp_val = np.asarray(kp(p, params), dtype=float)
        kt_val = np.asarray(kt(p, params), dtype=float)
    R_init = params.m * params.beta * I * C_I
    dCI = -params.beta * I * C_I
    dCR = R_init - 2.0 * kt_val * C_R * C_R - params.k_oxygen * C_O * C_R
    dCO = -params.k_oxygen * C_O * C_R
    dCM = -kp_val * C_R * C_M
    return dCI, dCR, dCO, dCM


def rhs(
    t: float,
    y: np.ndarray,
    params: KineticParams,
    intensity: IntensityFn,
    *,
    constant_rate: bool = False,
) -> np.ndarray:
    """Time derivatives of ``(C_I, C_R, C_O, C_M)``.

    ``constant_rate=True`` freezes ``k_p = k_{p0}`` and ``k_t = k_t(p=0)``
    so Slutzky's algebraic oxygen--monomer identity can be tested.
    """
    C_I, C_R, C_O, C_M = np.asarray(y, dtype=float)
    dCI, dCR, dCO, dCM = reaction_rates(
        C_I, C_R, C_O, C_M, intensity(t), params, constant_rate=constant_rate
    )
    return np.array([dCI, dCR, dCO, dCM], dtype=float)


def initial_state(params: KineticParams) -> np.ndarray:
    return np.array(
        [params.C_I0, params.C_R0, params.C_O0, params.C_M0],
        dtype=float,
    )
