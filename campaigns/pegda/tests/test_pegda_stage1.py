"""Stage 1 local ODE tests."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from numpy.testing import assert_allclose

from photopolymerization.kinetics.closures import kp, kt
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


def test_kp_decreases_with_conversion():
    params = _params()
    assert kp(0.0, params) > kp(0.8, params)


def test_kt_finite_at_p0_and_p1():
    params = _params()
    assert np.isfinite(kt(0.0, params))
    assert np.isfinite(kt(1.0, params))
    assert kt(0.0, params) > 0.0


def test_dark_no_conversion():
    params = _params()
    traj = integrate(params, ConstantIntensity(0.0), t_end=5.0, n_eval=50)
    assert traj.success
    assert_allclose(traj.p, 0.0, atol=1e-12)


def test_illuminated_consumes_initiator():
    params = _params()
    traj = integrate(params, ConstantIntensity(64.0), t_end=5.0, n_eval=80)
    assert traj.success
    assert traj.y[0, -1] < traj.y[0, 0]
