"""Stage 2 analytic-limit tests."""

from __future__ import annotations

from pathlib import Path

from photopolymerization.kinetics.limits import (
    qss_residual,
    slutzky_oxygen_monomer_residual,
)
from photopolymerization.kinetics.protocols import ConstantIntensity
from photopolymerization.kinetics.solve import integrate
from photopolymerization.parameters import kinetic_params_from_book, load_parameter_book

BOOK = (
    Path(__file__).resolve().parents[1]
    / "configs"
    / "parameter_book.montgomery_2022.json"
)


def _params():
    return kinetic_params_from_book(load_parameter_book(BOOK))


def test_slutzky_identity_constant_rate():
    params = _params()
    traj = integrate(params, ConstantIntensity(64.0), 8.0, constant_rate=True)
    residual = slutzky_oxygen_monomer_residual(traj, params)
    assert residual < 1e-4


def test_qss_oxygen_free_constant_rate():
    from dataclasses import replace

    params = replace(_params(), C_O0=0.0)
    traj = integrate(params, ConstantIntensity(64.0), 8.0, constant_rate=True)
    residual = qss_residual(traj, params, t_skip=0.5)
    assert residual < 5e-3


def test_montgomery_kp_differs_from_constant_fork():
    from photopolymerization.kinetics.closures import kp

    params = _params()
    assert kp(0.8, params) < 0.99 * kp(0.0, params)
