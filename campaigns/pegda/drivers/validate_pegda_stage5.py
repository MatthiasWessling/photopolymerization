#!/usr/bin/env python3
"""Stage 5 gates: fast/slow flow, Eq. 3, Uc scaling, 2% predicate labelled."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))

from photopolymerization.flow.params import load_slutzky_book, slutzky_params_from_book
from photopolymerization.flow.plug import (
    PI,
    critical_speed,
    gelled,
    integrate_plug,
    kinematic_fiber_length,
    tau_exp,
)
from photopolymerization.kinetics.limits import slutzky_eq3_residual

from _paths import config_path


def _gate(name: str, ok: bool, detail: str) -> bool:
    status = "PASS" if ok else "FAIL"
    print(f"  [{status}] {name}: {detail}")
    return ok


def _p(book, case: dict):
    return slutzky_params_from_book(
        book, kd=float(case["kd"]), Pi0=float(case["Pi0"]), kt=float(case["kt"])
    )


def validate(*, smoke: bool) -> int:
    cfg = json.loads(config_path("config.pegda_stage5.json").read_text())
    book = load_slutzky_book(config_path(cfg["parameter_book"]))
    print("Stage 5 — Slutzky plug flow (no wall O2, not Zhu P_gel)")
    results = []
    results.append(
        _gate(
            "chemistry_key",
            book["chemistry_key"] == cfg["chemistry_key_required"],
            book["chemistry_key"],
        )
    )
    results.append(
        _gate(
            "not_montgomery_book",
            "slutzky_2019" in book["chemistry_key"],
            "separate PEGDA-575 jet book",
        )
    )
    fig2 = _p(book, cfg["cases"]["fig2_high_oxygen"])
    fig3 = _p(book, cfg["cases"]["fig3_low_oxygen"])
    t2 = integrate_plug(fig2, float(cfg["cases"]["fig2_high_oxygen"]["U"]))
    t3 = integrate_plug(fig3, float(cfg["cases"]["fig3_low_oxygen"]["U"]))
    results.append(_gate("integrator_fig2", t2.success, t2.message))
    results.append(_gate("integrator_fig3", t3.success, t3.message))
    results.append(
        _gate(
            "fig2_no_gel",
            not gelled(t2, fig2),
            f"min M/M0={float(np.min(t2.M_over_M0)):.4f} (predicate 0.98, not P_gel)",
        )
    )
    results.append(
        _gate(
            "fig3_gel_two_percent",
            gelled(t3, fig3),
            f"min M/M0={float(np.min(t3.M_over_M0)):.4f}",
        )
    )

    r3 = slutzky_eq3_residual(
        t3.y[2], t3.y[3], fig3.y0, fig3.M0, fig3.ko, fig3.kp
    )
    eq_tol = 1e-3 if smoke else float(cfg["gates"]["eq3_rtol"])
    results.append(
        _gate("eq3_along_x", r3 < eq_tol, f"max rel residual={r3:.3e}")
    )

    if smoke:
        overall = all(results)
        print(f"OVERALL {'PASS' if overall else 'FAIL'}")
        return 0 if overall else 1

    t_fast = integrate_plug(fig3, float(cfg["cases"]["fast_flow"]["U"]))
    results.append(
        _gate(
            "fast_flow_negligible_conversion",
            t_fast.success
            and float(np.min(t_fast.M_over_M0)) > float(cfg["gates"]["fast_M_min"]),
            f"min M/M0={float(np.min(t_fast.M_over_M0)):.6f}",
        )
    )
    # Photoinitiator still consumed only in the window: downstream kd=0.
    iL = int(np.searchsorted(t3.x, fig3.L, side="left"))
    iL = min(iL, t3.x.size - 1)
    results.append(
        _gate(
            "Pi_frozen_downstream",
            abs(float(t3.y[PI, -1] - t3.y[PI, iL])) / fig3.Pi0 < 1e-6,
            "kd=0 for x>L",
        )
    )

    p_kd = slutzky_params_from_book(book, kd=0.002, Pi0=400.0, kt=1000.0)
    p_kd2 = slutzky_params_from_book(book, kd=0.008, Pi0=400.0, kt=1000.0)
    p_pi = slutzky_params_from_book(book, kd=0.005, Pi0=200.0, kt=1000.0)
    p_pi2 = slutzky_params_from_book(book, kd=0.005, Pi0=400.0, kt=1000.0)
    Uc_kd = critical_speed(p_kd)
    Uc_kd2 = critical_speed(p_kd2)
    Uc_pi = critical_speed(p_pi)
    Uc_pi2 = critical_speed(p_pi2)
    ratio = float(cfg["gates"]["uc_ratio_min"])
    results.append(
        _gate(
            "Uc_increases_with_kd",
            Uc_kd2 > ratio * Uc_kd,
            f"Uc(kd=0.002)={Uc_kd:.4f} Uc(kd=0.008)={Uc_kd2:.4f} m/s",
        )
    )
    results.append(
        _gate(
            "Uc_increases_with_Pi0",
            Uc_pi2 > ratio * Uc_pi,
            f"Uc(Pi0=200)={Uc_pi:.4f} Uc(Pi0=400)={Uc_pi2:.4f} m/s",
        )
    )

    U = float(cfg["pulses"]["U"])
    t_uv = cfg["pulses"]["t_uv_s"]
    L1 = kinematic_fiber_length(U, float(t_uv[0]))
    L2 = kinematic_fiber_length(U, float(t_uv[-1]))
    results.append(
        _gate(
            "pulse_length_increases_with_tuv",
            L2 > L1,
            f"U t_uv = {L1:.4f} -> {L2:.4f} m; tau_exp={tau_exp(float(t_uv[-1]), fig3, U):.3f}",
        )
    )
    print(
        "  [NOTE] Gel predicate is M/M0=0.98 (Slutzky). Zhu P_gel is Stage 6 and is off."
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
