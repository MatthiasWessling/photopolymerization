"""Gaussian pixel kernel, Montgomery Eq. 7 (fallback PSF).

Montgomery et al., Extreme Mech. Lett. 53, 101714 (2022):

``I_0(x,y) = sum_{i,j} I_pixel(g) exp[ -((x-x0_i)/(2σ))^2 - ((y-y0_j)/(2σ))^2 ]``

This is an optical boundary condition, not a polymer state.
A measured non-Gaussian kernel (Meenakshisundaram et al., 2020) is not used yet.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

SIGMA_MONTGOMERY_M = 1.77e-5
PITCH_FIG3_M = 50e-6
I_PIXEL_W_M2 = 41.08
I_IDEAL_W_M2 = 64.0


def gaussian_1d(x: NDArray[np.floating], x0: float, sigma: float) -> NDArray[np.floating]:
    """Montgomery 1-D factor ``exp(-((x-x0)/(2σ))^2)``."""
    return np.exp(-((x - x0) / (2.0 * sigma)) ** 2)


def irradiance_1d(
    x: NDArray[np.floating],
    centers: NDArray[np.floating],
    amplitudes: NDArray[np.floating],
    sigma: float,
) -> NDArray[np.floating]:
    """Superposition of 1-D Montgomery Gaussians (W/m^2 if amplitudes are)."""
    dx = (x[:, None] - centers[None, :]) / (2.0 * sigma)
    return (amplitudes[None, :] * np.exp(-(dx**2))).sum(axis=1)


def irradiance_2d(
    x: NDArray[np.floating],
    y: NDArray[np.floating],
    x0: NDArray[np.floating],
    y0: NDArray[np.floating],
    amplitudes: NDArray[np.floating],
    sigma: float,
) -> NDArray[np.floating]:
    """``I(x,y)`` on a meshgrid-like pair of 1-D axes, shape ``(ny, nx)``."""
    X, Y = np.meshgrid(x, y)
    field = np.zeros_like(X, dtype=float)
    for xc, yc, amp in zip(x0, y0, amplitudes):
        field += amp * np.exp(-((X - xc) / (2.0 * sigma)) ** 2 - ((Y - yc) / (2.0 * sigma)) ** 2)
    return field


def analytic_integral_1d(amplitude: float, sigma: float) -> float:
    r"""``\int exp(-((x)/(2σ))^2) dx = 2 σ \sqrt{π}``."""
    return amplitude * 2.0 * sigma * np.sqrt(np.pi)


def analytic_integral_2d(amplitude: float, sigma: float) -> float:
    r"""``\iint G\,dx\,dy = 4 π σ^2`` times peak amplitude (W if I in W/m^2)."""
    return amplitude * 4.0 * np.pi * sigma**2


def grayscale_to_amplitude(g: float, I_pixel: float) -> float:
    """Linear map: G0 full, G100 dark. SI calibration is still required."""
    g = float(np.clip(g, 0.0, 100.0))
    return I_pixel * (1.0 - g / 100.0)


def scan_dose_localized_1d(
    x: NDArray[np.floating],
    pattern_I: NDArray[np.floating],
    velocity_m_s: float,
) -> NDArray[np.floating]:
    """Dose from a compact 1-D pattern translating at speed ``v`` (non-periodic).

    A laboratory point at ``x[i]`` sees ``I_pattern(x[i] - v t)``. On a uniform
    grid this is a causal shift sum with ``dt = dx / v``.
    Peak dose scales as ``1/v`` (Meenakshisundaram dwell).
    """
    if velocity_m_s <= 0.0:
        raise ValueError("velocity must be positive")
    dx = float(np.mean(np.diff(x)))
    dt = dx / velocity_m_s
    n = x.size
    dose = np.zeros(n, dtype=float)
    for k in range(n):
        sl = slice(k, n)
        dose[sl] += pattern_I[: n - k] * dt
    return dose
