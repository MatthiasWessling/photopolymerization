#!/usr/bin/env python3
"""Stage 2 gates: Slutzky Eq. 3 and QSS radicals on the constant-rate fork."""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import replace
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))

from photopolymerization.kinetics.limits import qss_residual, slutzky_oxygen_monomer_residual
from photopolymerization.kinetics.protocols import ConstantIntensity
from photopolymerization.kinetics.solve import integrate
from photopolymerization.parameters import kinetic_params_from_book, load_parameter_book

from _paths import config_path


def _gate(name: str, ok: bool, detail: str) -> bool:
    status = "PASS" if ok else "FAIL"
    print(f"  [{status}] {name}: {detail}")
    return ok


def validate(*, smoke: bool) -> int:
    cfg = json.loads(config_path("config.pegda_stage2.json").read_text())
    book = load_parameter_book(config_path(cfg["parameter_book"]))
    print("Stage 2 — analytic limits")
    results = []
    results.append(
        _gate(
            "chemistry_key",
            book["chemistry_key"] == cfg["chemistry_key_required"],
            book["chemistry_key"],
        )
    )
    params = kinetic_params_from_book(book)
    I0 = float(cfg["solver"]["I_W_m2"])
    t_end = 4.0 if smoke else float(cfg["solver"]["t_end_s"])
    t_skip = 0.2 if smoke else float(cfg["solver"]["t_skip_qss_s"])
    slutzky_tol = 1e-3 if smoke else float(cfg["gates"]["slutzky_rtol"])
    qss_tol = 2e-2 if smoke else float(cfg["gates"]["qss_rtol"])
    intensity = ConstantIntensity(I0)

    const = integrate(params, intensity, t_end, constant_rate=True)
    results.append(_gate("integrator_constant", const.success, const.message))
    r_const = slutzky_oxygen_monomer_residual(const, params)
    results.append(
        _gate(
            "slutzky_eq3_constant_rate",
            r_const < slutzky_tol,
            f"max relative residual={r_const:.3e} (tol={slutzky_tol:.1e})",
        )
    )

    qss_traj = integrate(
        replace(params, C_O0=0.0),
        intensity,
        t_end,
        constant_rate=True,
    )
    results.append(_gate("integrator_qss", qss_traj.success, qss_traj.message))
    r_qss = qss_residual(qss_traj, params, t_skip=t_skip)
    results.append(
        _gate(
            "qss_radicals_oxygen_free",
            r_qss < qss_tol,
            f"max relative residual after t={t_skip}s is {r_qss:.3e} (tol={qss_tol:.1e})",
        )
    )

    if not smoke:
        from photopolymerization.kinetics.closures import kp as kp_fn

        mont = integrate(params, intensity, t_end, constant_rate=False)
        r_mont = slutzky_oxygen_monomer_residual(mont, params)
        print(
            f"  [INFO] slutzky_eq3_montgomery_oxygen_window: residual={r_mont:.3e} "
            "(expected similar while p≈0 during induction; not a fail)"
        )
        ratio = float(kp_fn(0.8, params) / kp_fn(0.0, params))
        results.append(
            _gate(
                "conversion_dependent_kp_is_not_constant",
                ratio < 0.99,
                f"kp(0.8)/kp(0)={ratio:.4f} (fork exists; Eq. 3 is not imposed here)",
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
