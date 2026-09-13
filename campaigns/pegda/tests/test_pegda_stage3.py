"""Stage 3 projection tests (no kinetics)."""

from __future__ import annotations

import numpy as np
from numpy.testing import assert_allclose
from scipy.integrate import trapezoid

from photopolymerization.projection.psf import (
    analytic_integral_1d,
    gaussian_1d,
    grayscale_to_amplitude,
    irradiance_1d,
    scan_dose_localized_1d,
)


def test_analytic_1d_integral():
    sigma = 1.77e-5
    I0 = 41.08
    x = np.linspace(-12 * sigma, 12 * sigma, 2001)
    field = I0 * gaussian_1d(x, 0.0, sigma)
    assert_allclose(trapezoid(field, x), analytic_integral_1d(I0, sigma), rtol=1e-4)


def test_dark_bitmap():
    x = np.linspace(-1e-4, 1e-4, 51)
    field = irradiance_1d(x, np.array([0.0]), np.array([0.0]), 1.77e-5)
    assert_allclose(field, 0.0)


def test_g80_linear_map():
    assert_allclose(grayscale_to_amplitude(80.0, 41.08), 0.2 * 41.08)


def test_scan_peak_dose_inverse_velocity():
    sigma = 1.77e-5
    x = np.linspace(-8e-4, 8e-4, 1601)
    field = 41.08 * gaussian_1d(x, 0.0, sigma)
    d1 = scan_dose_localized_1d(x, field, 0.002)
    d2 = scan_dose_localized_1d(x, field, 0.004)
    assert_allclose(np.max(d1) * 0.002, np.max(d2) * 0.004, rtol=2e-2)
