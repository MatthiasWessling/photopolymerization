#!/usr/bin/env python3
"""Stage 1 gates for the local Montgomery four-species ODE."""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))

from photopolymerization.kinetics.protocols import ConstantIntensity
from photopolymerization.kinetics.rhs import CI, CM, CO
from photopolymerization.kinetics.solve import Trajectory, integrate
from photopolymerization.parameters import kinetic_params_from_book, load_parameter_book

from _paths import config_path


def _gate(name: str, ok: bool, detail: str) -> bool:
    status = "PASS" if ok else "FAIL"
    print(f"  [{status}] {name}: {detail}")
    return ok


def _acrylate_balance(traj: Trajectory, C_M0: float) -> float:
    consumed = np.trapezoid(traj.Rp, traj.t)
    delta = C_M0 - traj.y[CM, -1]
    if abs(delta) < 1e-18:
        return abs(consumed)
    return abs(consumed - delta) / abs(delta)


def validate(*, smoke: bool) -> int:
    cfg = json.loads(config_path("config.pegda_stage1.json").read_text())
    book = load_parameter_book(config_path(cfg["parameter_book"]))
    if book["chemistry_key"] != cfg["chemistry_key_required"]:
        print("OVERALL FAIL")
        print("chemistry_key mismatch")
        return 1
    params = kinetic_params_from_book(book)
    proto = cfg["protocols"]
    t_end = 8.0 if smoke else float(cfg["solver"]["t_end_s"])
    I_high = float(proto["I_high_W_m2"])
    I_low = float(proto["I_low_W_m2"])
    t_high = min(float(proto["t_high_s"]), t_end)
    t_low = t_high * (I_high / I_low)

    print("Stage 1 — local constitutive response")
    results = []

    dark = integrate(params, ConstantIntensity(0.0), t_end)
    results.append(_gate("integrator_dark", dark.success, dark.message))
    results.append(
        _gate(
            "dark",
            float(np.max(np.abs(dark.p))) < 1e-10,
            f"max |p|={float(np.max(np.abs(dark.p))):.3e}",
        )
    )

    high = integrate(params, ConstantIntensity(I_high), t_end)
    results.append(_gate("integrator_high", high.success, high.message))
    ymin = float(np.min(high.y))
    results.append(
        _gate(
            "positivity",
            ymin >= -1e-12,
            f"min y={ymin:.3e} mol/m^3",
        )
    )
    results.append(
        _gate(
            "bounds",
            float(np.min(high.p)) >= -1e-10 and float(np.max(high.p)) <= 1.0 + 1e-10,
            f"p in [{float(np.min(high.p)):.4f}, {float(np.max(high.p)):.4f}]",
        )
    )
    dCI = np.diff(high.y[CI])
    results.append(
        _gate(
            "PI_monotone",
            bool(np.all(dCI <= 1e-12)),
            f"max dC_I={float(np.max(dCI)):.3e}",
        )
    )
    bal = _acrylate_balance(high, params.C_M0)
    results.append(
        _gate(
            "acrylate_balance",
            bal < (5e-2 if smoke else 1e-2),
            f"relative residual={bal:.3e}",
        )
    )

    if smoke:
        results.append(
            _gate(
                "conversion_on",
                float(high.p[-1]) > 0.0 or float(np.max(high.y[CO]) - high.y[CO, -1]) > 0.0,
                f"p_end={float(high.p[-1]):.4f}",
            )
        )
    else:
        high_eq = integrate(params, ConstantIntensity(I_high), t_high)
        low_eq = integrate(params, ConstantIntensity(I_low), t_low)
        p1 = float(high_eq.p[-1])
        p2 = float(low_eq.p[-1])
        results.append(
            _gate(
                "reciprocity",
                max(p1, p2) > 0.02 and abs(p1 - p2) > 0.01,
                f"p(I={I_high},{t_high}s)={p1:.4f}; p(I={I_low},{t_low}s)={p2:.4f}",
            )
        )
        ox_free = integrate(replace(params, C_O0=0.0), ConstantIntensity(I_high), t_end)

        def t_at(traj, p_star: float) -> float:
            hit = np.where(traj.p >= p_star)[0]
            if hit.size == 0:
                return float("inf")
            return float(traj.t[hit[0]])

        t_ox = t_at(high, 0.02)
        t_free = t_at(ox_free, 0.02)
        results.append(
            _gate(
                "induction",
                t_ox > t_free + 0.05,
                f"t(p=0.02) with O2={t_ox:.3f}s, without={t_free:.3f}s",
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
