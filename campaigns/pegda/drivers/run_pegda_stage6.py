#!/usr/bin/env python3
"""Stage 6: Zhu network maps of conversion. No bleaching ODEs, no kp refit."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))

from photopolymerization.kinetics.protocols import ConstantIntensity
from photopolymerization.kinetics.solve import integrate
from photopolymerization.network.macosko import (
    G_ideal,
    G_shear,
    eta,
    eta_noloop,
    map_conversion,
    network_params_from_book,
    p_gel,
)
from photopolymerization.parameters import kinetic_params_from_book, load_parameter_book

from _paths import FIGURES, config_path, output_dir


def main() -> None:
    cfg = json.loads(config_path("config.pegda_stage6.json").read_text())
    book = load_parameter_book(config_path(cfg["parameter_book"]))
    net = network_params_from_book(book)
    c_demo = float(cfg["c_demo_volpct"])
    p = np.linspace(0.0, 1.0, 201)
    eta_loop = eta(p, c_demo, net, loops=True)
    eta_free = eta_noloop(p)
    G_loop = G_shear(eta_loop, c_demo, net)
    G_id = G_ideal(c_demo, net)

    cs = np.array(cfg["c_sweep_volpct"], dtype=float)
    pg = np.array([p_gel(c, net) for c in cs])

    mont_book = load_parameter_book(config_path("parameter_book.montgomery_2022.json"))
    kin = kinetic_params_from_book(mont_book)
    traj = integrate(kin, ConstantIntensity(64.0), 20.0, n_eval=80)
    hypo = map_conversion(traj.p, c_demo, net, loops=True, hypothetical=True)

    summary = {
        "chemistry_key": net.chemistry_key,
        "eosin_Y_bleaching_odes": cfg["eosin_Y_bleaching_odes"],
        "C0_volpct": net.C0_volpct,
        "p_gel_20pct": p_gel(c_demo, net),
        "eta_p1_loops": float(eta_loop[-1]),
        "eta_p1_noloop": float(eta_free[-1]),
        "G_over_Gideal_p1": float(G_loop[-1] / G_id),
        "hypothetical_montgomery": True,
        "kp_refit": False,
    }
    out = output_dir("pegda_stage6")
    (out / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    fig, axes = plt.subplots(2, 2, figsize=(8.6, 6.2), constrained_layout=True)
    ax = axes[0, 0]
    ax.plot(cs, np.clip(pg, 0, 1.05), "o-", color="tab:red", lw=2)
    ax.axhline(1.0, color="0.5", ls="--", lw=0.8)
    ax.axvline(net.C0_volpct, color="tab:green", ls=":", lw=1.2, label=r"$C_0$")
    ax.set_xlabel(r"$[\mathrm{PEGDA}]_0$ / vol%")
    ax.set_ylabel(r"$P_{\mathrm{gel}}$")
    ax.set_title("Gel conversion vs concentration")
    ax.legend(frameon=False)

    ax = axes[0, 1]
    ax.plot(p, eta_free, color="tab:gray", lw=2, ls="--", label=r"no loops $\eta=p^2$")
    ax.plot(p, eta_loop, color="tab:blue", lw=2, label=rf"loops, {c_demo:.0f} vol%")
    ax.axvline(p_gel(c_demo, net), color="tab:red", ls=":", lw=1.0, label=r"$P_{\mathrm{gel}}$")
    ax.set_xlabel("acrylate conversion p")
    ax.set_ylabel(r"$\eta$")
    ax.set_title("Elastically active fraction")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1, 0]
    ax.plot(p, G_loop / 1000.0, color="tab:purple", lw=2, label="loops")
    ax.axhline(G_id / 1000.0, color="tab:gray", ls="--", lw=1.2, label=r"$G_{\mathrm{ideal}}$")
    ax.set_xlabel("p")
    ax.set_ylabel(r"$G$ / kPa")
    ax.set_title("Phantom modulus (Eq. 38)")
    ax.legend(frameon=False)

    ax = axes[1, 1]
    ax.plot(traj.t, traj.p, color="tab:orange", lw=2, label="Montgomery p(t) [hyp.]")
    ax.plot(traj.t, hypo.eta, color="tab:blue", lw=2, label=r"Zhu $\eta(p)$ [hyp.]")
    ax.set_xlabel("t / s")
    ax.set_ylabel("p or η")
    ax.set_title("Hypothetical Type-I p → Type-II network map")
    ax.legend(frameon=False, fontsize=8)

    fig.suptitle("Stage 6 Zhu network post-processor (eosin Y ODEs off)", fontsize=11)
    fig.savefig(out / "pegda_stage6_overview.pdf")
    fig.savefig(FIGURES / "pegda_stage6_overview.pdf")
    fig.savefig(FIGURES / "pegda_stage6_overview.png")
    plt.close(fig)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
