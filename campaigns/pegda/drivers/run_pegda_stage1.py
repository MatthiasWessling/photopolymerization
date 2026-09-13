#!/usr/bin/env python3
"""Integrate Stage 1 protocols and write timeseries + overview figure."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import cumulative_trapezoid

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))

from dataclasses import replace

from photopolymerization.kinetics.protocols import ConstantIntensity
from photopolymerization.kinetics.rhs import CI, CM, CO, CR
from photopolymerization.kinetics.solve import integrate
from photopolymerization.parameters import kinetic_params_from_book, load_parameter_book

from _paths import FIGURES, config_path, output_dir


def _save_csv(path: Path, traj) -> None:
    header = "t_s,I_W_m2,C_I,C_R,C_O,C_M,p,Rp"
    data = np.column_stack(
        [
            traj.t,
            traj.I,
            traj.y[CI],
            traj.y[CR],
            traj.y[CO],
            traj.y[CM],
            traj.p,
            traj.Rp,
        ]
    )
    np.savetxt(path, data, delimiter=",", header=header, comments="")


def main() -> None:
    cfg = json.loads(config_path("config.pegda_stage1.json").read_text())
    book = load_parameter_book(config_path(cfg["parameter_book"]))
    params = kinetic_params_from_book(book)
    proto = cfg["protocols"]
    t_end = float(cfg["solver"]["t_end_s"])
    I_high = float(proto["I_high_W_m2"])
    I_low = float(proto["I_low_W_m2"])

    cases = {
        "dark": integrate(params, ConstantIntensity(0.0), t_end),
        "high_I": integrate(params, ConstantIntensity(I_high), t_end),
        "low_I": integrate(params, ConstantIntensity(I_low), t_end),
        "oxygen_free": integrate(
            replace(params, C_O0=0.0),
            ConstantIntensity(I_high),
            t_end,
        ),
    }
    out = output_dir("pegda_stage1")
    for name, traj in cases.items():
        if not traj.success:
            raise RuntimeError(f"{name} integrator failed: {traj.message}")
        _save_csv(out / f"{name}.csv", traj)

    fig, axes = plt.subplots(2, 2, figsize=(8.4, 6.2), constrained_layout=True)
    ax = axes[0, 0]
    for name, color in (("high_I", "tab:orange"), ("low_I", "tab:blue"), ("dark", "tab:gray")):
        ax.plot(cases[name].t, cases[name].p, color=color, lw=2, label=name)
    ax.set_xlabel("t / s")
    ax.set_ylabel("conversion p")
    ax.legend(frameon=False)
    ax.set_title("Conversion vs time")

    ax = axes[0, 1]
    dose_h = cumulative_trapezoid(cases["high_I"].I, cases["high_I"].t, initial=0.0)
    dose_l = cumulative_trapezoid(cases["low_I"].I, cases["low_I"].t, initial=0.0)
    ax.plot(dose_h, cases["high_I"].p, color="tab:orange", lw=2, label=f"{I_high} W/m$^2$")
    ax.plot(dose_l, cases["low_I"].p, color="tab:blue", lw=2, label=f"{I_low} W/m$^2$")
    ax.set_xlabel("nominal dose $\\int I\\,dt$ / J m$^{-2}$")
    ax.set_ylabel("conversion p")
    ax.legend(frameon=False)
    ax.set_title("Equal-dose comparison (schematic of history dependence)")

    ax = axes[1, 0]
    ax.plot(cases["high_I"].t, cases["high_I"].y[CO], color="tab:red", lw=2, label="with O2")
    ax.plot(
        cases["oxygen_free"].t,
        cases["oxygen_free"].y[CO],
        color="tab:green",
        lw=2,
        ls="--",
        label="C_O0=0",
    )
    ax.set_xlabel("t / s")
    ax.set_ylabel("$C_{O_2}$ / mol m$^{-3}$")
    ax.legend(frameon=False)
    ax.set_title("Oxygen inventory")

    ax = axes[1, 1]
    ax.plot(cases["high_I"].t, cases["high_I"].p, color="tab:orange", lw=2, label="with O2")
    ax.plot(
        cases["oxygen_free"].t,
        cases["oxygen_free"].p,
        color="tab:green",
        lw=2,
        ls="--",
        label="C_O0=0",
    )
    ax.set_xlabel("t / s")
    ax.set_ylabel("conversion p")
    ax.legend(frameon=False)
    ax.set_title("Oxygen induction")

    fig.suptitle(params.chemistry_key, fontsize=11)
    fig.savefig(out / "pegda_stage1_overview.pdf")
    fig.savefig(FIGURES / "pegda_stage1_overview.pdf")
    fig.savefig(FIGURES / "pegda_stage1_overview.png")
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
