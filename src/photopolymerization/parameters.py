"""Parameter book with paper provenance.

Values for the Type-I PEGDA-250 / Irgacure 819 branch are transcribed from
Montgomery et al., Extreme Mech. Lett. 53, 101714 (2022), Table 2.
They are not universal PEGDA constants.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

import json

ALLOWED_STATUS = frozenset({"transcribed", "calculated", "still_required", "assumed"})
ALLOWED_INITIATOR = frozenset({"type_I", "type_II"})


@dataclass(frozen=True)
class ProvenanceEntry:
    name: str
    value: float
    units: str
    source: str
    status: str
    chemistry_key: str


@dataclass(frozen=True)
class KineticParams:
    """Numeric values used by the Stage 1 ODE (SI except beta, which follows Montgomery).

    ``beta`` multiplies irradiance in W/m^2 and has units s^2/kg so that
    ``beta * I`` is a first-order photolysis frequency (1/s).
    """

    chemistry_key: str
    initiator_class: str
    m: float
    beta: float
    k_oxygen: float
    kp0: float
    kp_D0: float
    c_p: float
    kt_SD: float
    kt_TD0: float
    C_RD: float
    C_I0: float
    C_R0: float
    C_O0: float
    C_M0: float


@dataclass(frozen=True)
class OpticalParams:
    """Beer-Lambert coefficients, Montgomery Eqs. 4--6 and Table 1.

    ``w_absorber`` multiplies ``A_absorber``. The paper calls ``w`` a weight
    percent; the multilayer example uses 0.06 wt% photoabsorber. FTIR-cell
    runs should set ``w_absorber=0``.
    """

    alpha_I: float
    A_monomer: float
    A_polymer: float
    A_absorber: float
    w_absorber: float
    I_ideal: float
    I_pixel: float
    sigma_m: float


@dataclass(frozen=True)
class DiffusionParams:
    """Harmonic liquid--solid diffusivities, Montgomery Eq. 28 and Table 3."""

    D_I_liquid: float
    D_I_solid: float
    D_R_liquid: float
    D_R_solid: float
    D_O_liquid: float
    D_O_solid: float
    D_M_liquid: float
    D_M_solid: float


def load_parameter_book(path: str | Path) -> dict[str, Any]:
    """Load a Stage 0 JSON book and check the provenance schema."""
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    _validate_book(data)
    return data


def _validate_book(data: Mapping[str, Any]) -> None:
    required_top = (
        "chemistry_key",
        "initiator_class",
        "formulation_label",
        "do_not_mix_with",
        "parameters",
    )
    missing = [k for k in required_top if k not in data]
    if missing:
        raise ValueError(f"parameter book missing keys: {missing}")
    if data["initiator_class"] not in ALLOWED_INITIATOR:
        raise ValueError(f"unknown initiator_class: {data['initiator_class']}")
    params = data["parameters"]
    if not isinstance(params, dict) or not params:
        raise ValueError("parameters must be a non-empty object")
    for name, entry in params.items():
        for field in ("value", "units", "source", "status"):
            if field not in entry:
                raise ValueError(f"{name} missing {field}")
        if entry["status"] not in ALLOWED_STATUS:
            raise ValueError(f"{name} has illegal status {entry['status']}")
        if not isinstance(entry["value"], (int, float)):
            raise ValueError(f"{name} value must be numeric")


def kinetic_params_from_book(book: Mapping[str, Any]) -> KineticParams:
    """Extract the ODE vector from a provenance book."""
    p = book["parameters"]

    def v(name: str) -> float:
        return float(p[name]["value"])

    return KineticParams(
        chemistry_key=str(book["chemistry_key"]),
        initiator_class=str(book["initiator_class"]),
        m=v("m"),
        beta=v("beta"),
        k_oxygen=v("k_oxygen"),
        kp0=v("kp0"),
        kp_D0=v("kp_D0"),
        c_p=v("c_p"),
        kt_SD=v("kt_SD"),
        kt_TD0=v("kt_TD0"),
        C_RD=v("C_RD"),
        C_I0=v("C_I0"),
        C_R0=v("C_R0"),
        C_O0=v("C_O0"),
        C_M0=v("C_M0"),
    )


def optical_params_from_book(
    book: Mapping[str, Any],
    *,
    w_absorber: float | None = None,
) -> OpticalParams:
    """Extract Montgomery Table 1 optics. Optional ``w_absorber`` override."""
    p = book["parameters"]

    def v(name: str) -> float:
        return float(p[name]["value"])

    w = float(w_absorber) if w_absorber is not None else v("w_absorber")
    return OpticalParams(
        alpha_I=v("alpha_I"),
        A_monomer=v("A_monomer"),
        A_polymer=v("A_polymer"),
        A_absorber=v("A_absorber"),
        w_absorber=w,
        I_ideal=v("I_ref"),
        I_pixel=v("I_pixel"),
        sigma_m=v("sigma_pixel"),
    )


def diffusion_params_from_book(book: Mapping[str, Any]) -> DiffusionParams:
    """Extract Montgomery Table 3 diffusivities."""
    p = book["parameters"]

    def v(name: str) -> float:
        return float(p[name]["value"])

    return DiffusionParams(
        D_I_liquid=v("D_I_liquid"),
        D_I_solid=v("D_I_solid"),
        D_R_liquid=v("D_R_liquid"),
        D_R_solid=v("D_R_solid"),
        D_O_liquid=v("D_O_liquid"),
        D_O_solid=v("D_O_solid"),
        D_M_liquid=v("D_M_liquid"),
        D_M_solid=v("D_M_solid"),
    )


def still_required_names(book: Mapping[str, Any]) -> list[str]:
    return sorted(
        name
        for name, entry in book["parameters"].items()
        if entry["status"] == "still_required"
    )
