"""Zhu 2020 loop-aware Macosko--Miller network (post-processor only).

Does not advance species. Eosin Y bleaching ODEs are not used.
PEGDA is treated as difunctional in acrylate groups (f=2), with two
carbon ends per reacted acrylate (Zhu Eqs. 32--35).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping

import numpy as np
from numpy.typing import NDArray
from scipy.special import comb

K_BOLTZMANN = 1.380649e-23
N_AVOGADRO = 6.02214076e23
F_ACRYLATE = 2


@dataclass(frozen=True)
class NetworkParams:
    """Zhu Table A2 network constants. Concentrations as volume percent."""

    chemistry_key: str
    C0_volpct: float
    T: float
    f: int = F_ACRYLATE


def network_params_from_book(book: Mapping[str, Any]) -> NetworkParams:
    p = book["parameters"]
    return NetworkParams(
        chemistry_key=str(book["chemistry_key"]),
        C0_volpct=float(p["C0_volpct"]["value"]),
        T=float(p["T"]["value"]),
        f=int(p["f"]["value"]),
    )


def theta(c_pegda_volpct: float, params: NetworkParams) -> float:
    """Probability of intermolecular (not loop) propagation, Zhu Eq. 21."""
    c = float(c_pegda_volpct)
    return c / (c + 2.0 * params.C0_volpct)


def p_gel(c_pegda_volpct: float, params: NetworkParams) -> float:
    """Acrylate conversion at gelation, Zhu Eq. 26. ``>=1`` means never gels."""
    th = theta(c_pegda_volpct, params)
    if th <= 0.0:
        return float("inf")
    return (1.0 - th) / (2.0 * th)


def can_gel(c_pegda_volpct: float, params: NetworkParams) -> bool:
    """False when ``[PEGDA]_0 <= C_0`` (Zhu after Eq. 26)."""
    return float(c_pegda_volpct) > params.C0_volpct


def _phi_F_out_star(p: NDArray[np.floating], th: float) -> NDArray[np.floating]:
    """``φ(F_out^{CC*})`` from Zhu Eqs. 27--31.

    ``F=1`` (all looks finite) is always a root. After gelation the physical
    root is the smaller one in ``[0, 1)``.
    """
    p = np.clip(np.asarray(p, dtype=float), 0.0, 1.0)
    if th >= 1.0 - 1e-15:
        return np.zeros_like(p)
    if th <= 0.0:
        return np.ones_like(p)
    a = th * p * p
    b = -(1.0 - th + th * p * p)
    c = 1.0 - th
    F = np.ones_like(p)
    mask = a > 1e-18
    disc = b * b - 4.0 * a * c
    good = mask & (disc >= 0.0)
    sqrt_d = np.sqrt(np.clip(disc[good], 0.0, None))
    F1 = (-b[good] - sqrt_d) / (2.0 * a[good])
    F2 = (-b[good] + sqrt_d) / (2.0 * a[good])
    cand = np.where(F1 < F2, F1, F2)
    cand = np.clip(cand, 0.0, 1.0)
    F[good] = cand
    return F


def eta_noloop(p: NDArray[np.floating] | float) -> NDArray[np.floating] | float:
    """Zhu Eq. 36: ``η = P_CC^2`` when loops are forbidden."""
    p = np.asarray(p, dtype=float)
    return p * p


def eta(
    p: NDArray[np.floating] | float,
    c_pegda_volpct: float,
    params: NetworkParams,
    *,
    loops: bool = True,
) -> NDArray[np.floating]:
    """Fractional elastically active density, Zhu Eq. 35 (or 36 if ``loops=False``)."""
    p = np.clip(np.asarray(p, dtype=float, copy=True), 0.0, 1.0)
    scalar = p.ndim == 0
    p = np.atleast_1d(p)
    if not loops:
        out = eta_noloop(p)
        return out.reshape(()) if scalar else out
    if not can_gel(c_pegda_volpct, params):
        out = np.zeros_like(p)
        return out.reshape(()) if scalar else out
    th = theta(c_pegda_volpct, params)
    pg = p_gel(c_pegda_volpct, params)
    F = _phi_F_out_star(p, th)
    F = np.where(p + 1e-15 >= pg, F, 1.0)
    out = _eta_from_F(p, F, params.f)
    out = np.where(p + 1e-15 >= pg, out, 0.0)
    out = np.clip(out, 0.0, 1.0)
    return out.reshape(()) if scalar else out


def _eta_from_F(
    p: NDArray[np.floating],
    F: NDArray[np.floating],
    f: int,
) -> NDArray[np.floating]:
    """Zhu Eqs. 32--35 with ``f=2``."""
    eta_val = np.zeros_like(p, dtype=float)
    q = 1.0 - F
    for k in range(0, f + 1):
        pk = comb(f, k) * p**k * (1.0 - p) ** (f - k)
        n_ends = 2 * k
        for m in range(3, n_ends + 1):
            pm = comb(n_ends, m) * (F ** (n_ends - m)) * (q**m)
            eta_val += (2.0 / f) * ((m - 2) / 2.0) * pm * pk
    return eta_val


def c_pegda_mol_m3(c_volpct: float, *, rho: float = 1120.0, M: float = 0.575) -> float:
    """Volume-percent PEGDA-575 to mol/m^3 (density assumed, labelled)."""
    return (c_volpct / 100.0) * rho / M


def G_shear(
    eta_val: NDArray[np.floating] | float,
    c_pegda_volpct: float,
    params: NetworkParams,
) -> NDArray[np.floating] | float:
    """Phantom modulus ``G = η N_A [PEGDA]_0 k_B T``, Zhu Eq. 38."""
    c = c_pegda_mol_m3(c_pegda_volpct)
    return np.asarray(eta_val, dtype=float) * N_AVOGADRO * c * K_BOLTZMANN * params.T


def G_ideal(c_pegda_volpct: float, params: NetworkParams) -> float:
    """``G_ideal = N_A [PEGDA]_0 k_B T`` (``½ ν_0 kT`` with ``ν_0 = 2 N_A c``)."""
    return float(G_shear(1.0, c_pegda_volpct, params))


@dataclass(frozen=True)
class NetworkMap:
    p: NDArray[np.floating]
    eta: NDArray[np.floating]
    G: NDArray[np.floating]
    p_gel: float
    gelled: NDArray[np.bool_]
    hypothetical: bool


def map_conversion(
    p: NDArray[np.floating] | float,
    c_pegda_volpct: float,
    params: NetworkParams,
    *,
    loops: bool = True,
    hypothetical: bool = False,
) -> NetworkMap:
    """Map a conversion field to ``η`` and ``G``. Does not call kinetics."""
    p = np.asarray(p, dtype=float)
    et = eta(p, c_pegda_volpct, params, loops=loops)
    G = G_shear(et, c_pegda_volpct, params)
    pg = p_gel(c_pegda_volpct, params)
    return NetworkMap(
        p=p,
        eta=np.asarray(et, dtype=float),
        G=np.asarray(G, dtype=float),
        p_gel=pg,
        gelled=p >= pg if np.isfinite(pg) and pg < 1.0 else np.zeros_like(p, dtype=bool),
        hypothetical=hypothetical,
    )
