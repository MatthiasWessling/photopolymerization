#!/usr/bin/env python3
"""Stage 6 gates: P_gel, no gel at C0, eta=p^2, G<G_ideal, no eosin ODEs."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))

from photopolymerization.network.macosko import (
    G_ideal,
    G_shear,
    can_gel,
    eta,
    eta_noloop,
    network_params_from_book,
    p_gel,
)
from photopolymerization.parameters import load_parameter_book

from _paths import config_path


def _gate(name: str, ok: bool, detail: str) -> bool:
    status = "PASS" if ok else "FAIL"
    print(f"  [{status}] {name}: {detail}")
    return ok


def validate(*, smoke: bool) -> int:
    cfg = json.loads(config_path("config.pegda_stage6.json").read_text())
    book = load_parameter_book(config_path(cfg["parameter_book"]))
    print("Stage 6 — Zhu network post-processor (eosin Y ODEs off)")
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
            "eosin_odes_off",
            str(book.get("eosin_Y_bleaching_odes", "")).lower() == "off",
            book.get("eosin_Y_bleaching_odes"),
        )
    )
    net = network_params_from_book(book)
    C0 = net.C0_volpct
    results.append(
        _gate(
            "theta_to_1_pgel_to_0",
            p_gel(1e6, net) < 1e-5,
            f"P_gel(huge c)={p_gel(1e6, net):.3e}",
        )
    )
    results.append(
        _gate(
            "pgel_equals_C0_over_c",
            abs(p_gel(20.0, net) - C0 / 20.0) < 1e-12,
            f"P_gel(20%)={p_gel(20.0, net):.4f} C0/c={C0/20.0:.4f}",
        )
    )
    results.append(
        _gate(
            "no_gel_at_or_below_C0",
            (not can_gel(C0, net))
            and (not can_gel(0.5 * C0, net))
            and float(eta(1.0, C0, net, loops=True)) == 0.0,
            f"C0={C0} vol%, eta(p=1)={float(eta(1.0, C0, net)):.3e}",
        )
    )
    p = np.linspace(0.0, 1.0, 21 if smoke else 101)
    et0 = eta_noloop(p)
    results.append(
        _gate(
            "noloop_eta_is_p_squared",
            np.max(np.abs(et0 - p * p)) < float(cfg["gates"]["noloop_rtol"]),
            "Eq. 36",
        )
    )
    et_loop = eta(p, 20.0, net, loops=True)
    results.append(
        _gate(
            "eta_le_1",
            float(np.max(et_loop)) <= float(cfg["gates"]["eta_max"]) + 1e-12,
            f"max eta={float(np.max(et_loop)):.4f}",
        )
    )
    if smoke:
        overall = all(results)
        print(f"OVERALL {'PASS' if overall else 'FAIL'}")
        return 0 if overall else 1

    results.append(
        _gate(
            "loops_reduce_eta_at_full_conversion",
            float(et_loop[-1]) < float(et0[-1]) - 0.02,
            f"eta_loops(p=1)={float(et_loop[-1]):.4f} eta_noloop={float(et0[-1]):.4f}",
        )
    )
    G = float(G_shear(et_loop[-1], 20.0, net))
    Gi = G_ideal(20.0, net)
    results.append(
        _gate(
            "G_below_ideal",
            G < Gi,
            f"G={G:.3e} Pa G_ideal={Gi:.3e} Pa G/G_ideal={G/Gi:.3f}",
        )
    )
    results.append(
        _gate(
            "pregel_eta_zero",
            float(eta(0.5 * p_gel(20.0, net), 20.0, net)) == 0.0,
            "eta=0 for p < P_gel",
        )
    )
    results.append(
        _gate(
            "slutzky_2pct_is_not_pgel",
            abs(p_gel(20.0, net) - 0.02) > 0.05,
            f"Zhu P_gel(20%)={p_gel(20.0, net):.3f} vs Slutzky 0.02",
        )
    )
    print("  [NOTE] Mapping Montgomery p(t) through this eta is hypothetical, not matched.")
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
