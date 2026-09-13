#!/usr/bin/env python3
"""Stage 0 gates: provenance schema and no silent recipe mixing."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))

from photopolymerization.parameters import (
    kinetic_params_from_book,
    load_parameter_book,
    still_required_names,
)

from _paths import config_path


def _gate(name: str, ok: bool, detail: str) -> bool:
    status = "PASS" if ok else "FAIL"
    print(f"  [{status}] {name}: {detail}")
    return ok


def validate(*, smoke: bool) -> int:
    cfg = json.loads(config_path("config.pegda_stage0.json").read_text())
    book = load_parameter_book(config_path(cfg["parameter_book"]))
    required_ode = (
        "m",
        "beta",
        "k_oxygen",
        "kp0",
        "kp_D0",
        "c_p",
        "kt_SD",
        "kt_TD0",
        "C_RD",
        "C_I0",
        "C_R0",
        "C_O0",
        "C_M0",
    )
    print("Stage 0 — parameter book")
    results = []
    results.append(
        _gate(
            "schema",
            True,
            f"loaded {book['chemistry_key']} with {len(book['parameters'])} entries",
        )
    )
    missing = [k for k in required_ode if k not in book["parameters"]]
    results.append(
        _gate("ode_keys", not missing, "ok" if not missing else f"missing {missing}")
    )
    results.append(
        _gate(
            "type_I_flag",
            book["initiator_class"] == "type_I",
            book["initiator_class"],
        )
    )
    mix = book.get("do_not_mix_with", [])
    results.append(
        _gate(
            "do_not_mix",
            len(mix) >= 2,
            f"{len(mix)} foreign chemistries listed",
        )
    )
    params = kinetic_params_from_book(book)
    results.append(
        _gate(
            "positive_inventory",
            params.C_M0 > 0 and params.C_I0 > 0 and params.C_O0 >= 0,
            f"C_M0={params.C_M0}, C_I0={params.C_I0}, C_O0={params.C_O0}",
        )
    )
    still = still_required_names(book)
    if smoke:
        results.append(
            _gate(
                "still_required_logged",
                True,
                f"{still} (smoke: allowed)",
            )
        )
    else:
        results.append(
            _gate(
                "oxygen_flagged",
                "C_O0" in still,
                "C_O0 must remain still_required until measured",
            )
        )
        results.append(
            _gate(
                "beta_units_note",
                "s^2/kg" in book["parameters"]["beta"]["units"],
                book["parameters"]["beta"]["units"],
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
