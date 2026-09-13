#!/usr/bin/env python3
"""Stage 2: constant-rate identities vs Montgomery conversion-dependent RHS."""

from __future__ import annotations

import json
import sys
from dataclasses import replace
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))

from photopolymerization.kinetics.limits import qss_radical, slutzky_oxygen_monomer_residual
from photopolymerization.kinetics.protocols import ConstantIntensity
from photopolymerization.kinetics.solve import integrate
from photopolymerization.parameters import kinetic_params_from_book, load_parameter_book

from _paths import FIGURES, config_path, output_dir


def main() -> None:
    cfg = json.loads(config_path("config.pegda_stage2.json").read_text())
    params = kinetic_params_from_book(
        load_parameter_book(config_path(cfg["parameter_book"]))
    )
    I0 = float(cfg["solver"]["I_W_m2"])
    t_end = float(cfg["solver"]["t_end_s"])
    intensity = ConstantIntensity(I0)

    const = integrate(params, intensity, t_end, constant_rate=True)
    mont = integrate(params, intensity, t_end, constant_rate=False)
    qss = integrate(
        replace(params, C_O0=0.0),
        intensity,
        t_end,
        constant_rate=True,
    )
    out = output_dir("pegda_stage2")
    np.savez(
        out / "trajectories.npz",
        t=const.t,
        C_O_const=const.y[2],
        C_M_const=const.y[3],
        C_R_qss=qss.y[1],
        C_I_qss=qss.y[0],
        I_qss=qss.I,
        C_O_mont=mont.y[2],
        C_M_mont=mont.y[3],
    )

    frac_M = const.y[3] / params.C_M0
    frac_O = const.y[2] / params.C_O0
    predicted = np.power(np.maximum(frac_M, 1e-30), params.k_oxygen / params.kp0)

    fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.7), constrained_layout=True)
    ax = axes[0]
    ax.plot(const.t, frac_O, color="tab:red", lw=2, label=r"$C_O/C_{O,0}$")
    ax.plot(
        const.t,
        predicted,
        color="tab:blue",
        lw=2,
        ls="--",
        label=r"$(C_M/C_{M,0})^{k_O/k_p}$",
    )
    ax.set_xlabel("t / s")
    ax.set_ylabel("oxygen fraction")
    ax.legend(frameon=False)
    ax.set_title("Slutzky Eq. 3 (constant-rate fork)")

    ax = axes[1]
    pred_R = qss_radical(qss.y[0], qss.I, params)
    ax.plot(qss.t, qss.y[1], color="tab:orange", lw=2, label=r"$C_R$")
    ax.plot(qss.t, pred_R, color="tab:green", lw=2, ls="--", label=r"QSS $\sqrt{R_i/2k_t}$")
    ax.set_xlabel("t / s")
    ax.set_ylabel(r"$C_R$ / mol m$^{-3}$")
    ax.legend(frameon=False)
    ax.set_title("Oxygen-free QSS radicals")

    fig.suptitle("Stage 2 analytic limits (not a Montgomery FTIR match)", fontsize=11)
    fig.savefig(out / "pegda_stage2_overview.pdf")
    fig.savefig(FIGURES / "pegda_stage2_overview.pdf")
    fig.savefig(FIGURES / "pegda_stage2_overview.png")
    plt.close(fig)

    summary = {
        "slutzky_residual_constant": slutzky_oxygen_monomer_residual(const, params),
        "slutzky_residual_montgomery": slutzky_oxygen_monomer_residual(mont, params),
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
