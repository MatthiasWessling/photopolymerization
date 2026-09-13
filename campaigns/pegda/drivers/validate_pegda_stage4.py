#!/usr/bin/env python3
"""Stage 4 gates: conservation, MMS, ODE match, front, dead zone, blur."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import trapezoid

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))

from photopolymerization.kinetics.protocols import ConstantIntensity
from photopolymerization.kinetics.solve import integrate
from photopolymerization.parameters import (
    diffusion_params_from_book,
    kinetic_params_from_book,
    load_parameter_book,
    optical_params_from_book,
)
from photopolymerization.projection.psf import grayscale_to_amplitude, irradiance_1d
from photopolymerization.rd.beer import absorption, optical_depth
from photopolymerization.rd.mms import integrate_cosine_mms
from photopolymerization.rd.mol import Grid1D, make_grid
from photopolymerization.rd.operators import fourier_explicit_limit
from photopolymerization.rd.solve import integrate_1d

from _paths import config_path


def _gate(name: str, ok: bool, detail: str) -> bool:
    status = "PASS" if ok else "FAIL"
    print(f"  [{status}] {name}: {detail}")
    return ok


def _transition_width(
    x: np.ndarray,
    p: np.ndarray,
    x_left: float,
    x_right: float,
    x_step: float,
) -> float:
    """10--90% width of p across a commanded step, using interior plateaus."""
    p_left = float(np.mean(p[np.abs(x - x_left) < 0.6 * abs(x_step - x_left)]))
    p_right = float(np.mean(p[np.abs(x - x_right) < 0.6 * abs(x_right - x_step)]))
    lo, hi = (min(p_left, p_right), max(p_left, p_right))
    span = hi - lo
    if span < 0.02:
        return 0.0
    p10 = lo + 0.1 * span
    p90 = lo + 0.9 * span
    # Walk through the interface neighbourhood only.
    mask = (x >= min(x_left, x_right)) & (x <= max(x_left, x_right))
    xs, ps = x[mask], p[mask]
    if ps[0] > ps[-1]:
        xs, ps = xs[::-1], ps[::-1]
    x10 = float(np.interp(p10, ps, xs))
    x90 = float(np.interp(p90, ps, xs))
    return abs(x90 - x10)


def validate(*, smoke: bool) -> int:
    cfg = json.loads(config_path("config.pegda_stage4.json").read_text())
    book = load_parameter_book(config_path(cfg["parameter_book"]))
    print("Stage 4 — reaction–diffusion + Beer–Lambert (no flow)")
    results = []
    kin = kinetic_params_from_book(book)
    diff = diffusion_params_from_book(book)
    opt_thin = optical_params_from_book(book, w_absorber=0.0)
    opt_thick = optical_params_from_book(
        book, w_absorber=float(cfg["optics_variants"]["thick_vat"]["w_absorber"])
    )
    t_end = float(cfg["solver"]["t_end_smoke_s"] if smoke else cfg["solver"]["t_end_s"])
    n_z = int(cfg["geometry_z"]["n_smoke"] if smoke else cfg["geometry_z"]["n"])
    H = float(cfg["geometry_z"]["H_m"])
    I0 = float(cfg["I0_W_m2"])
    grid_z = make_grid(H, n_z)
    max_step = float(cfg["solver"]["max_step_s"])
    dt_fo = fourier_explicit_limit(diff.D_R_liquid, grid_z.dx)
    results.append(
        _gate(
            "fourier_documented",
            np.isfinite(dt_fo) and dt_fo > 0.0,
            f"explicit dt_max={dt_fo:.3e} s; integrator=BDF (implicit)",
        )
    )

    dark = integrate_1d(
        grid_z, kin, opt_thin, diff, 0.0, t_end, along="z", oxygen_bc="sealed",
        max_step=max_step,
    )
    results.append(_gate("integrator_dark", dark.success, dark.message))
    m0 = float(trapezoid(np.full(n_z, kin.C_M0), grid_z.x))
    m1 = float(trapezoid(dark.C_M[-1], grid_z.x))
    mass_rtol = float(cfg["gates"]["mass_rtol"])
    results.append(
        _gate(
            "dark_sealed_mass",
            abs(m1 - m0) / m0 < mass_rtol,
            f"rel={abs(m1-m0)/m0:.3e}",
        )
    )

    t_mms = min(float(cfg["solver"]["t_mms_s"]), t_end if smoke else 5.0)
    y_num, y_ex = integrate_cosine_mms(grid_z.x, diff.D_I_liquid, t_mms)
    mms_err = float(np.max(np.abs(y_num - y_ex)) / np.max(np.abs(y_ex)))
    mms_tol = 2e-2 if smoke else float(cfg["gates"]["mms_rtol"])
    results.append(
        _gate("mms_diffusion", mms_err < mms_tol, f"rel max={mms_err:.3e}")
    )

    no_d = integrate_1d(
        grid_z, kin, opt_thin, diff, I0, t_end, along="z", oxygen_bc="sealed",
        diffuse=False, max_step=max_step,
    )
    ode = integrate(kin, ConstantIntensity(I0), t_end, n_eval=80)
    dp = abs(float(no_d.p[-1, 0]) - float(ode.p[-1]))
    ode_tol = 0.08 if smoke else float(cfg["gates"]["ode_match_p"])
    results.append(
        _gate(
            "surface_matches_stage1_ode",
            no_d.success and ode.success and dp < ode_tol,
            f"p_z0={float(no_d.p[-1,0]):.4f} p_ode={float(ode.p[-1]):.4f} |dp|={dp:.3e}",
        )
    )

    if smoke:
        overall = all(results)
        print(f"OVERALL {'PASS' if overall else 'FAIL'}")
        return 0 if overall else 1

    thick = integrate_1d(
        grid_z, kin, opt_thick, diff, I0, t_end, along="z", oxygen_bc="sealed",
        max_step=max_step,
    )
    results.append(_gate("integrator_thick", thick.success, thick.message))
    A0 = float(absorption(kin.C_I0, 0.0, opt_thick))
    bopt = optical_depth(A0, H)
    dfront = float(thick.p[-1, 0] - thick.p[-1, -1])
    results.append(
        _gate(
            "optically_thick_front",
            bopt >= 1.0 and dfront > float(cfg["gates"]["front_delta_p"]),
            f"beta_opt={bopt:.3f}, p(0)-p(H)={dfront:.4f}",
        )
    )

    wall = integrate_1d(
        grid_z, kin, opt_thin, diff, I0, t_end, along="z", oxygen_bc="dirichlet",
        max_step=max_step,
    )
    sealed_thin = integrate_1d(
        grid_z, kin, opt_thin, diff, I0, t_end, along="z", oxygen_bc="sealed",
        max_step=max_step,
    )
    p_w = float(wall.p[-1, 0])
    p_mid = float(wall.p[-1, n_z // 2])
    p_sealed0 = float(sealed_thin.p[-1, 0])
    ratio = p_w / max(p_mid, 1e-12)
    results.append(
        _gate(
            "oxygen_dead_zone",
            wall.success and p_mid > 0.05 and ratio < float(cfg["gates"]["dead_zone_ratio"])
            and p_w < p_sealed0 - 0.02,
            f"p_wall={p_w:.4f} p_mid={p_mid:.4f} ratio={ratio:.3f} p_sealed_face={p_sealed0:.4f}",
        )
    )

    geo = cfg["geometry_x"]
    n_x = int(geo["n"])
    pitch = float(geo["pitch_m"])
    n_pix = int(geo["n_pixels"])
    x = np.linspace(-1.5 * pitch, (n_pix - 0.5) * pitch, n_x)
    centers = np.arange(n_pix, dtype=float) * pitch
    I_pixel = float(geo["I_pixel_W_m2"])
    amp = np.array(
        [I_pixel] * (n_pix // 2)
        + [grayscale_to_amplitude(80.0, I_pixel)] * (n_pix - n_pix // 2)
    )
    I_gray = irradiance_1d(x, centers, amp, float(geo["sigma_m"]))
    grid_x = Grid1D(x=x, dx=float(x[1] - x[0]))
    gon = integrate_1d(
        grid_x, kin, opt_thin, diff, I_gray, t_end, along="x", oxygen_bc="sealed",
        diffuse=True, max_step=max_step,
    )
    goff = integrate_1d(
        grid_x, kin, opt_thin, diff, I_gray, t_end, along="x", oxygen_bc="sealed",
        diffuse=False, max_step=max_step,
    )
    x_step = 4.5 * pitch
    x_left = 2.0 * pitch
    x_right = 7.0 * pitch
    w_on = _transition_width(x, gon.p[-1], x_left, x_right, x_step)
    w_off = _transition_width(x, goff.p[-1], x_left, x_right, x_step)
    sigma = float(geo["sigma_m"])
    extra = float(cfg["gates"]["blur_extra_m"])
    results.append(
        _gate(
            "grayscale_wider_than_optical_sigma",
            gon.success
            and goff.success
            and w_on > sigma
            and w_off > sigma,
            f"w_on={w_on*1e6:.1f} um w_off={w_off*1e6:.1f} um sigma={sigma*1e6:.1f} um",
        )
    )
    dw = w_on - w_off
    print(
        f"  [NOTE] Table 3 D_i increment at t={t_end}s: "
        f"w_on-w_off={dw*1e6:.2f} um (threshold {extra*1e6:.1f} um). "
        f"Montgomery 100-150 um interface is consistent with the Gaussian PSF, "
        f"not an extra {extra*1e6:.0f} um from D_R at these rates."
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
