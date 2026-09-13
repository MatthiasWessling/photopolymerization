#!/usr/bin/env python3
"""Stage 5 plug-flow figure: high/low oxygen profiles, U sweep, pulse lengths."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))

from photopolymerization.flow.params import load_slutzky_book, scaling_Uc, slutzky_params_from_book
from photopolymerization.flow.plug import (
    critical_speed,
    gelled,
    integrate_plug,
    kinematic_fiber_length,
    tau_exp,
)

from _paths import FIGURES, config_path, output_dir


def _params(book, case: dict):
    return slutzky_params_from_book(
        book, kd=float(case["kd"]), Pi0=float(case["Pi0"]), kt=float(case["kt"])
    )


def main() -> None:
    cfg = json.loads(config_path("config.pegda_stage5.json").read_text())
    book = load_slutzky_book(config_path(cfg["parameter_book"]))
    fig2 = _params(book, cfg["cases"]["fig2_high_oxygen"])
    fig3 = _params(book, cfg["cases"]["fig3_low_oxygen"])
    t2 = integrate_plug(fig2, float(cfg["cases"]["fig2_high_oxygen"]["U"]))
    t3 = integrate_plug(fig3, float(cfg["cases"]["fig3_low_oxygen"]["U"]))
    t_fast = integrate_plug(fig3, float(cfg["cases"]["fast_flow"]["U"]))
    for name, traj in (("fig2", t2), ("fig3", t3), ("fast", t_fast)):
        if not traj.success:
            raise RuntimeError(f"{name}: {traj.message}")

    U_grid = np.geomspace(3e-4, 3e-2, 16)
    gel_flag = []
    Mmin = []
    for U in U_grid:
        tr = integrate_plug(fig3, float(U), n_eval=200)
        gel_flag.append(gelled(tr, fig3))
        Mmin.append(float(np.min(tr.M_over_M0)))

    kd_lo, kd_hi = 0.002, 0.008
    p_lo = slutzky_params_from_book(book, kd=kd_lo, Pi0=400.0, kt=1000.0)
    p_hi = slutzky_params_from_book(book, kd=kd_hi, Pi0=400.0, kt=1000.0)
    Uc_lo = critical_speed(p_lo)
    Uc_hi = critical_speed(p_hi)

    pulses = cfg["pulses"]
    U_pulse = float(pulses["U"])
    lengths = [kinematic_fiber_length(U_pulse, t) for t in pulses["t_uv_s"]]
    taus = [tau_exp(t, fig3, U_pulse) for t in pulses["t_uv_s"]]

    summary = {
        "chemistry_key": book["chemistry_key"],
        "gel_predicate": cfg["gel_predicate"],
        "fig2_min_M_over_M0": float(np.min(t2.M_over_M0)),
        "fig3_min_M_over_M0": float(np.min(t3.M_over_M0)),
        "fast_min_M_over_M0": float(np.min(t_fast.M_over_M0)),
        "Uc_kd_lo": Uc_lo,
        "Uc_kd_hi": Uc_hi,
        "scaling_Uc_fig3": scaling_Uc(fig3),
        "pulse_lengths_m": lengths,
        "tau_exp": taus,
        "montgomery_rates": "not_used",
        "zhu_Pgel": "not_used",
    }
    out = output_dir("pegda_stage5")
    (out / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    fig, axes = plt.subplots(2, 2, figsize=(8.6, 6.2), constrained_layout=True)
    ax = axes[0, 0]
    X2 = t2.x / fig2.L
    X3 = t3.x / fig3.L
    ax.plot(X2, t2.y[2] / fig2.y0, color="tab:blue", lw=2, label="Y, Fig. 2")
    ax.plot(X3, t3.y[2] / fig3.y0, color="tab:red", lw=2, label="Y, Fig. 3")
    ax.axvline(1.0, color="0.5", lw=0.8, ls="--")
    ax.set_xlabel(r"$X=x/L$")
    ax.set_ylabel(r"$Y=C_{O_2}/C_{O,0}$")
    ax.set_title("Oxygen along the jet")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[0, 1]
    ax.plot(X2, t2.M_over_M0, color="tab:blue", lw=2, label="Fig. 2 high O2")
    ax.plot(X3, t3.M_over_M0, color="tab:red", lw=2, label="Fig. 3 low O2")
    ax.axhline(0.98, color="tab:green", lw=1.0, ls=":", label=r"$M/M_0=0.98$ (not $P_{\mathrm{gel}}$)")
    ax.set_xlabel(r"$X=x/L$")
    ax.set_ylabel(r"$M/M_0$")
    ax.set_title("Double bonds (2% gel predicate)")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1, 0]
    ax.semilogx(U_grid, Mmin, color="tab:purple", lw=2)
    ax.axhline(0.98, color="tab:green", lw=1.0, ls=":")
    ax.set_xlabel(r"$U$ / m s$^{-1}$")
    ax.set_ylabel(r"min $M/M_0$")
    ax.set_title("Fig. 3 chemistry vs speed")

    ax = axes[1, 1]
    ax.plot(pulses["t_uv_s"], np.array(lengths) * 1e3, "o-", color="tab:orange", lw=2)
    ax.set_xlabel(r"$t_{\mathrm{uv}}$ / s")
    ax.set_ylabel(r"$U t_{\mathrm{uv}}$ / mm")
    ax.set_title(r"Kinematic fiber length (uniform regime)")

    fig.suptitle("Stage 5 Slutzky plug flow (no wall O2 diffusion)", fontsize=11)
    fig.savefig(out / "pegda_stage5_overview.pdf")
    fig.savefig(FIGURES / "pegda_stage5_overview.pdf")
    fig.savefig(FIGURES / "pegda_stage5_overview.png")
    plt.close(fig)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
