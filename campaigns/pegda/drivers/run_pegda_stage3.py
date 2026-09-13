#!/usr/bin/env python3
"""Stage 3 projection figure (Montgomery Fig. 3 style). Chemistry is not called."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import trapezoid

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))

from photopolymerization.projection.psf import (
    analytic_integral_1d,
    gaussian_1d,
    grayscale_to_amplitude,
    irradiance_1d,
    scan_dose_localized_1d,
)

from _paths import FIGURES, config_path, output_dir


def main() -> None:
    cfg = json.loads(config_path("config.pegda_stage3.json").read_text())
    sigma = float(cfg["sigma_m"])
    pitch = float(cfg["pitch_m"])
    I_pixel = float(cfg["I_pixel_W_m2"])
    x = np.linspace(-1.5 * pitch, 10.5 * pitch, 2001)
    centers = np.arange(10, dtype=float) * pitch

    one = I_pixel * gaussian_1d(x, 0.0, sigma)
    ten_g0 = irradiance_1d(x, centers, np.full(10, I_pixel), sigma)
    amp_g40 = np.array([I_pixel] * 5 + [grayscale_to_amplitude(40.0, I_pixel)] * 5)
    amp_g80 = np.array([I_pixel] * 5 + [grayscale_to_amplitude(80.0, I_pixel)] * 5)
    ten_g40 = irradiance_1d(x, centers, amp_g40, sigma)
    ten_g80 = irradiance_1d(x, centers, amp_g80, sigma)

    fig, axes = plt.subplots(2, 2, figsize=(8.6, 6.0), constrained_layout=True)
    panels = [
        (axes[0, 0], one, "Single pixel (G0)", True),
        (axes[0, 1], ten_g0, "Ten pixels, all G0", True),
        (axes[1, 0], ten_g40, "Five G0 + five G40", True),
        (axes[1, 1], ten_g80, "Five G0 + five G80", True),
    ]
    for ax, field, title, _ in panels:
        ax.plot(x * 1e6, field, color="tab:blue", lw=1.8)
        for i in range(11):
            ax.axvline((i - 0.5) * pitch * 1e6, color="tab:red", lw=0.8, alpha=0.6)
        ax.set_xlabel(r"$x$ / $\mu$m")
        ax.set_ylabel(r"$I_\perp$ / W m$^{-2}$")
        ax.set_title(title)

    fig.suptitle(
        "Stage 3 Gaussian PSF (Montgomery Eq. 7 fallback; not a measured DMD kernel)",
        fontsize=11,
    )
    out = output_dir("pegda_stage3")
    fig.savefig(out / "pegda_stage3_overview.pdf")
    fig.savefig(FIGURES / "pegda_stage3_overview.pdf")
    fig.savefig(FIGURES / "pegda_stage3_overview.png")
    plt.close(fig)

    num = trapezoid(one, x)
    summary = {
        "analytic_1d_integral": analytic_integral_1d(I_pixel, sigma),
        "numeric_1d_integral": float(num),
        "kernel": cfg["kernel"],
        "chemistry": "not_used",
    }
    v1, v2 = 0.002, 0.004
    dose1 = scan_dose_localized_1d(x, one, v1)
    dose2 = scan_dose_localized_1d(x, one, v2)
    summary["scan_peak_dose_times_v"] = {
        "v1": float(np.max(dose1) * v1),
        "v2": float(np.max(dose2) * v2),
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
