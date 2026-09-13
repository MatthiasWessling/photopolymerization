"""Stage 4 reaction--diffusion tests."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from numpy.testing import assert_allclose
from scipy.integrate import trapezoid

from photopolymerization.parameters import (
    diffusion_params_from_book,
    kinetic_params_from_book,
    load_parameter_book,
    optical_params_from_book,
)
from photopolymerization.rd.beer import absorption, irradiance_along_z
from photopolymerization.rd.diffusivity import harmonic_D
from photopolymerization.rd.mms import integrate_cosine_mms
from photopolymerization.rd.mol import make_grid
from photopolymerization.rd.solve import integrate_1d

BOOK = (
    Path(__file__).resolve().parents[1]
    / "configs"
    / "parameter_book.montgomery_2022.json"
)


def test_harmonic_limits():
    assert_allclose(harmonic_D(0.0, 1e-10, 1e-14), 1e-10)
    assert_allclose(harmonic_D(1.0, 1e-10, 1e-14), 1e-14)


def test_beer_no_absorption():
    z = np.linspace(0.0, 1e-4, 11)
    I = irradiance_along_z(z, np.zeros_like(z), 64.0)
    assert_allclose(I, 64.0)


def test_mms_cosine():
    x = np.linspace(0.0, 2e-4, 41)
    y_num, y_ex = integrate_cosine_mms(x, 1e-10, 2.0)
    rel = np.max(np.abs(y_num - y_ex)) / np.max(np.abs(y_ex))
    assert rel < 1e-2


def test_dark_mass_conservation():
    book = load_parameter_book(BOOK)
    kin = kinetic_params_from_book(book)
    opt = optical_params_from_book(book, w_absorber=0.0)
    diff = diffusion_params_from_book(book)
    grid = make_grid(2e-4, 17)
    traj = integrate_1d(
        grid, kin, opt, diff, 0.0, 1.0, along="z", oxygen_bc="sealed", max_step=0.2
    )
    assert traj.success
    m0 = trapezoid(np.full(grid.n, kin.C_M0), grid.x)
    m1 = trapezoid(traj.C_M[-1], grid.x)
    assert abs(m1 - m0) / m0 < 1e-6


def test_absorption_increases_with_polymer():
    book = load_parameter_book(BOOK)
    opt = optical_params_from_book(book, w_absorber=0.0)
    a0 = absorption(2.651, 0.0, opt)
    a1 = absorption(2.651, 1.0, opt)
    assert a1 > a0
