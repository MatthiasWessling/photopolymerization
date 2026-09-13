"""Slutzky 2019 plug-flow photopolymerization (PEGDA 575 jet).

Not Montgomery PEGDA-250 / Irgacure 819. Do not mix books.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping

from photopolymerization.parameters import load_parameter_book


@dataclass(frozen=True)
class SlutzkyParams:
    """Steady plug-flow coefficients, Soft Matter Eqs. 1a--d and Table 2.

    Photolysis frequency ``kd`` is already ``phi * epsilon * I`` (1/s).
    Termination is ``k_t R^2`` as written (no extra factor of two).
    One radical per initiator photolysis (their Table 1, reaction I).
    """

    chemistry_key: str
    L: float
    kd: float
    kt: float
    kp: float
    ko: float
    Pi0: float
    M0: float
    y0: float
    gel_M_over_M0: float


def slutzky_params_from_book(
    book: Mapping[str, Any],
    *,
    kd: float | None = None,
    Pi0: float | None = None,
    kt: float | None = None,
) -> SlutzkyParams:
    p = book["parameters"]

    def v(name: str) -> float:
        return float(p[name]["value"])

    return SlutzkyParams(
        chemistry_key=str(book["chemistry_key"]),
        L=v("L"),
        kd=float(kd) if kd is not None else v("kd_ref"),
        kt=float(kt) if kt is not None else v("kt"),
        kp=v("kp"),
        ko=v("ko"),
        Pi0=float(Pi0) if Pi0 is not None else v("Pi0_ref"),
        M0=v("M0"),
        y0=v("y0"),
        gel_M_over_M0=v("gel_M_over_M0"),
    )


def load_slutzky_book(path: str) -> dict[str, Any]:
    return load_parameter_book(path)


def photolysis_frequency(phi: float, epsilon: float, I_einstein: float) -> float:
    """``k_d = phi * epsilon * I`` with ``I`` in Einstein/(m^2 s)."""
    return float(phi * epsilon * I_einstein)


def scaling_Uc(params: SlutzkyParams) -> float:
    """``U_c ~ L k_d Pi0 / y0`` (their oxygen-depletion estimate)."""
    return params.L * params.kd * params.Pi0 / params.y0
