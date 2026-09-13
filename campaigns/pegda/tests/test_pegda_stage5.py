"""Stage 5 Slutzky plug-flow tests."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from numpy.testing import assert_allclose

from photopolymerization.flow.params import (
    load_slutzky_book,
    photolysis_frequency,
    slutzky_params_from_book,
)
from photopolymerization.flow.plug import gelled, integrate_plug, kinematic_fiber_length
from photopolymerization.kinetics.limits import slutzky_eq3_residual

BOOK = (
    Path(__file__).resolve().parents[1]
    / "configs"
    / "parameter_book.slutzky_2019.json"
)


def test_kd_from_einstein_intensity():
    book = load_slutzky_book(BOOK)
    p = book["parameters"]
    kd = photolysis_frequency(p["phi"]["value"], p["epsilon"]["value"], 1.3e-3)
    assert_allclose(kd, 1.3e-3 * 0.6 * 1.6, rtol=1e-12)


def test_fig2_no_gel():
    params = slutzky_params_from_book(
        load_slutzky_book(BOOK), kd=1e-3, Pi0=200.0, kt=1000.0
    )
    traj = integrate_plug(params, 0.003)
    assert traj.success
    assert not gelled(traj, params)


def test_fig3_gel():
    params = slutzky_params_from_book(
        load_slutzky_book(BOOK), kd=5e-3, Pi0=400.0, kt=1000.0
    )
    traj = integrate_plug(params, 0.003)
    assert traj.success
    assert gelled(traj, params)


def test_eq3_along_x():
    params = slutzky_params_from_book(
        load_slutzky_book(BOOK), kd=5e-3, Pi0=400.0, kt=1000.0
    )
    traj = integrate_plug(params, 0.003)
    r = slutzky_eq3_residual(
        traj.y[2], traj.y[3], params.y0, params.M0, params.ko, params.kp
    )
    assert r < 1e-4


def test_kinematic_length_scales_with_pulse():
    assert_allclose(kinematic_fiber_length(0.05, 0.04), 0.002)
    assert kinematic_fiber_length(0.05, 0.08) == 2 * kinematic_fiber_length(0.05, 0.04)
