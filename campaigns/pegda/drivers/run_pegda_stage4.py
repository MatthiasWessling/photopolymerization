#!/usr/bin/env python3
"""Stage 4 overview: depth front, oxygen wall, grayscale blur. No flow."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))

from photopolymerization.parameters import (
    diffusion_params_from_book,
    kinetic_params_from_book,
    load_parameter_book,
    optical_params_from_book,
)
from photopolymerization.projection.psf import grayscale_to_amplitude, irradiance_1d
from photopolymerization.rd.beer import absorption, optical_depth
from photopolymerization.rd.mol import Grid1D, make_grid
from photopolymerization.rd.operators import fourier_explicit_limit
from photopolymerization.rd.solve import integrate_1d

from _paths import FIGURES, config_path, output_dir


def _grayscale_I(cfg: dict, n: int) -> tuple[np.ndarray, np.ndarray]:
    geo = cfg["geometry_x"]
    pitch = float(geo["pitch_m"])
    n_pix = int(geo["n_pixels"])
    x = np.linspace(-1.5 * pitch, (n_pix - 0.5) * pitch, n)
    centers = np.arange(n_pix, dtype=float) * pitch
    I_pixel = float(geo["I_pixel_W_m2"])
    amp = np.array(
        [I_pixel] * (n_pix // 2)
        + [grayscale_to_amplitude(80.0, I_pixel)] * (n_pix - n_pix // 2)
    )
    field = irradiance_1d(x, centers, amp, float(geo["sigma_m"]))
    return x, field


def main() -> None:
    cfg = json.loads(config_path("config.pegda_stage4.json").read_text())
    book = load_parameter_book(config_path(cfg["parameter_book"]))
    kin = kinetic_params_from_book(book)
    diff = diffusion_params_from_book(book)
    opt_thick = optical_params_from_book(
        book, w_absorber=float(cfg["optics_variants"]["thick_vat"]["w_absorber"])
    )
    opt_thin = optical_params_from_book(
        book, w_absorber=float(cfg["optics_variants"]["ftir_like"]["w_absorber"])
    )
    t_end = float(cfg["solver"]["t_end_s"])
    I0 = float(cfg["I0_W_m2"])
    n_z = int(cfg["geometry_z"]["n"])
    H = float(cfg["geometry_z"]["H_m"])
    grid_z = make_grid(H, n_z)
    max_step = float(cfg["solver"]["max_step_s"])

    sealed = integrate_1d(
        grid_z, kin, opt_thick, diff, I0, t_end, along="z", oxygen_bc="sealed",
        max_step=max_step,
    )
    wall = integrate_1d(
        grid_z, kin, opt_thin, diff, I0, t_end, along="z", oxygen_bc="dirichlet",
        max_step=max_step,
    )
    sealed_thin = integrate_1d(
        grid_z, kin, opt_thin, diff, I0, t_end, along="z", oxygen_bc="sealed",
        max_step=max_step,
    )
    n_x = int(cfg["geometry_x"]["n"])
    x_line, I_gray = _grayscale_I(cfg, n_x)
    grid_x = Grid1D(x=x_line, dx=float(x_line[1] - x_line[0]))
    gray_on = integrate_1d(
        grid_x, kin, opt_thin, diff, I_gray, t_end, along="x", oxygen_bc="sealed",
        diffuse=True, max_step=max_step,
    )
    gray_off = integrate_1d(
        grid_x, kin, opt_thin, diff, I_gray, t_end, along="x", oxygen_bc="sealed",
        diffuse=False, max_step=max_step,
    )
    for name, traj in (
        ("sealed_thick", sealed),
        ("oxygen_wall", wall),
        ("grayscale_D_on", gray_on),
        ("grayscale_D_off", gray_off),
    ):
        if not traj.success:
            raise RuntimeError(f"{name}: {traj.message}")

    A0 = float(absorption(kin.C_I0, 0.0, opt_thick))
    dx_z = grid_z.dx
    summary = {
        "beta_opt_uncured_thick": optical_depth(A0, H),
        "dt_fourier_explicit_s": fourier_explicit_limit(diff.D_R_liquid, dx_z),
        "integrator": "BDF implicit; Fourier limit is documentation only",
        "p_illuminated_face_thick": float(sealed.p[-1, 0]),
        "p_far_face_thick": float(sealed.p[-1, -1]),
        "p_wall_dirichlet": float(wall.p[-1, 0]),
        "p_mid_dirichlet": float(wall.p[-1, n_z // 2]),
        "p_mid_sealed_thin": float(sealed_thin.p[-1, n_z // 2]),
        "flow": "not_used",
    }
    out = output_dir("pegda_stage4")
    (out / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    fig, axes = plt.subplots(2, 2, figsize=(8.6, 6.2), constrained_layout=True)
    ax = axes[0, 0]
    ax.plot(sealed.x * 1e6, sealed.p[-1], color="tab:orange", lw=2, label="sealed, absorber")
    ax.set_xlabel(r"$z$ / $\mu$m")
    ax.set_ylabel("conversion p")
    ax.set_title("Optically thick depth profile")
    ax.legend(frameon=False)

    ax = axes[0, 1]
    ax.plot(wall.x * 1e6, wall.p[-1], color="tab:red", lw=2, label="O2 Dirichlet at z=0")
    ax.plot(
        sealed_thin.x * 1e6,
        sealed_thin.p[-1],
        color="tab:gray",
        lw=2,
        ls="--",
        label="sealed, no absorber",
    )
    ax.set_xlabel(r"$z$ / $\mu$m")
    ax.set_ylabel("conversion p")
    ax.set_title("Oxygen-wall dead zone (thin optics)")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1, 0]
    ax.plot(x_line * 1e6, gray_off.p[-1], color="tab:blue", lw=2, label="D=0")
    ax.plot(x_line * 1e6, gray_on.p[-1], color="tab:green", lw=2, label="D(p) on")
    ax.set_xlabel(r"$x$ / $\mu$m")
    ax.set_ylabel("conversion p")
    ax.set_title(r"G0|G80 grayscale (thin film in $x$)")
    ax.legend(frameon=False)

    ax = axes[1, 1]
    ax.plot(sealed.x * 1e6, sealed.I[-1], color="tab:purple", lw=2)
    ax.set_xlabel(r"$z$ / $\mu$m")
    ax.set_ylabel(r"$I$ / W m$^{-2}$")
    ax.set_title("Beer–Lambert I(z) at t_end")

    fig.suptitle("Stage 4 reaction–diffusion (no flow)", fontsize=11)
    fig.savefig(out / "pegda_stage4_overview.pdf")
    fig.savefig(FIGURES / "pegda_stage4_overview.pdf")
    fig.savefig(FIGURES / "pegda_stage4_overview.png")
    plt.close(fig)
    print(json.dumps(summary, indent=2))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
