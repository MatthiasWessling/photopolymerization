"""Conversion-dependent propagation and termination.

Montgomery et al., Extreme Mech. Lett. 53, 101714 (2022), Eqs. 22--27
and Table 2.  The paper writes ``exp(cp)`` for ``exp(c p)`` with
``c`` the relative viscosity coefficient (Table 2).
"""

from __future__ import annotations

import numpy as np

from photopolymerization.parameters import KineticParams


def conversion(C_M: float | np.ndarray, C_M0: float) -> float | np.ndarray:
    """Acrylate conversion ``p = 1 - C_M / C_M0``."""
    return 1.0 - np.asarray(C_M, dtype=float) / C_M0


def kp(p: float | np.ndarray, params: KineticParams) -> float | np.ndarray:
    """Propagation coefficient, Montgomery Eq. 26."""
    p = np.asarray(p, dtype=float)
    return params.kp0 / (
        1.0 + (params.kp0 / params.kp_D0) * np.exp(params.c_p * p)
    )


def kt(p: float | np.ndarray, params: KineticParams) -> float | np.ndarray:
    """Termination coefficient ``k_{t,D} + k_{t,RD}``, Montgomery Eqs. 22--27."""
    p = np.asarray(p, dtype=float)
    kt_D_inv = (1.0 / params.kt_SD) + np.exp(params.c_p * p) / params.kt_TD0
    kt_D = 1.0 / kt_D_inv
    kt_RD = params.C_RD * kp(p, params) * (1.0 - p)
    return kt_D + kt_RD
