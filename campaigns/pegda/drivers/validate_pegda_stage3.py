#!/usr/bin/env python3
"""Stage 3 gates: PSF integral, bleed, dark bitmap, scan 1/v. No chemistry."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import trapezoid

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))

from photopolymerization.projection.psf import (
    analytic_integral_1d,
    analytic_integral_2d,
    gaussian_1d,
    grayscale_to_amplitude,
    irradiance_1d,
    irradiance_2d,
    scan_dose_localized_1d,
)

from _paths import config_path


def _gate(name: str, ok: bool, detail: str) -> bool:
    status = "PASS" if ok else "FAIL"
    print(f"  [{status}] {name}: {detail}")
    return ok


def validate(*, smoke: bool) -> int:
    cfg = json.loads(config_path("config.pegda_stage3.json").read_text())
    print("Stage 3 — projection engine (chemistry off)")
    results = []
    results.append(
        _gate(
            "kernel_labelled_fallback",
            cfg["kernel"] == "gaussian_montgomery_eq7",
            cfg["kernel_status"],
        )
    )
    sigma = float(cfg["sigma_m"])
    pitch = float(cfg["pitch_m"])
    I_pixel = float(cfg["I_pixel_W_m2"])
    n = 401 if smoke else 1201
    half = 8 * sigma if smoke else 12 * sigma
    x = np.linspace(-half, half, n)
    one = I_pixel * gaussian_1d(x, 0.0, sigma)
    num = float(trapezoid(one, x))
    ana = analytic_integral_1d(I_pixel, sigma)
    rtol = 5e-3 if smoke else float(cfg["gates"]["integral_rtol"])
    results.append(
        _gate(
            "single_pixel_integral_1d",
            abs(num - ana) / ana < rtol,
            f"numeric={num:.6e} analytic={ana:.6e} rel={abs(num-ana)/ana:.3e}",
        )
    )

    y = x
    field = irradiance_2d(
        x, y, np.array([0.0]), np.array([0.0]), np.array([I_pixel]), sigma
    )
    num2 = float(trapezoid(trapezoid(field, x), y))
    ana2 = analytic_integral_2d(I_pixel, sigma)
    results.append(
        _gate(
            "single_pixel_integral_2d",
            abs(num2 - ana2) / ana2 < rtol,
            f"numeric={num2:.6e} analytic={ana2:.6e} rel={abs(num2-ana2)/ana2:.3e}",
        )
    )

    dark = irradiance_1d(x, np.array([0.0]), np.array([0.0]), sigma)
    results.append(_gate("dark_bitmap", float(np.max(np.abs(dark))) == 0.0, "I=0"))

    x_line = np.linspace(-1.5 * pitch, 10.5 * pitch, 2001)
    centers = np.arange(10, dtype=float) * pitch
    ten = irradiance_1d(x_line, centers, np.full(10, I_pixel), sigma)
    left_edge = -0.5 * pitch
    I_edge = float(np.interp(left_edge, x_line, ten))
    I_mid = float(np.interp(4.5 * pitch, x_line, ten))
    results.append(
        _gate(
            "ten_pixel_bleed",
            I_edge > 0.05 * I_pixel and I_edge < 0.95 * I_mid,
            f"I(outer edge)={I_edge:.3f}, I(interior)={I_mid:.3f} W/m^2",
        )
    )
    g80 = grayscale_to_amplitude(80.0, I_pixel)
    results.append(
        _gate(
            "grayscale_not_chemistry",
            abs(g80 - 0.2 * I_pixel) < 1e-12,
            f"G80 amplitude={g80:.4f} (linear map, assumed; no p reported)",
        )
    )

    if not smoke:
        v1, v2 = 0.002, 0.004
        dose1 = scan_dose_localized_1d(x_line, ten, v1)
        dose2 = scan_dose_localized_1d(x_line, ten, v2)
        q1 = float(np.max(dose1) * v1)
        q2 = float(np.max(dose2) * v2)
        scan_tol = float(cfg["gates"]["scan_vE_rtol"])
        results.append(
            _gate(
                "scan_dose_scales_as_1_over_v",
                abs(q1 - q2) / max(q1, q2) < scan_tol,
                f"max(E)*v = {q1:.6e} vs {q2:.6e}",
            )
        )

    overall = all(results)
    print(f"OVERALL {'PASS' if overall else 'FAIL'}")
    return 0 if overall else 1


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()
    raise SystemExit(validate(smoke=args.smoke))


if __name__ == "__main__":
    main()
