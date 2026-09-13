"""Schematic figures for the photopolymerization background note.

These panels are ontology cartoons, not experimental data.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle
import numpy as np

OUT = Path(__file__).resolve().parent
plt.rcParams.update(
    {
        "font.size": 10,
        "axes.labelsize": 11,
        "figure.dpi": 150,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)


def _box(ax, xy, w, h, text, fc="#e8f1fb", ec="#1f4e79"):
    patch = FancyBboxPatch(
        xy,
        w,
        h,
        boxstyle="round,pad=0.02,rounding_size=0.08",
        facecolor=fc,
        edgecolor=ec,
        linewidth=1.2,
    )
    ax.add_patch(patch)
    ax.text(
        xy[0] + w / 2,
        xy[1] + h / 2,
        text,
        ha="center",
        va="center",
        wrap=True,
    )


def fig_species():
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 6.2)
    ax.axis("off")
    ax.set_title("Chemical inventory (schematic)")

    _box(ax, (0.3, 4.4), 2.6, 1.4, "Photoinitiator\n(ground state)", fc="#fff4d6")
    _box(ax, (3.7, 4.4), 2.6, 1.4, "Primary radicals\n$R^\\bullet$")
    _box(ax, (7.1, 4.4), 2.6, 1.4, "Growing chains\n$P^\\bullet$")
    _box(ax, (0.3, 2.2), 2.6, 1.4, "Acrylate groups\n$M$ (PEGDA)", fc="#e7f6e7")
    _box(ax, (3.7, 2.2), 2.6, 1.4, "Dissolved $\\mathrm{O}_2$\n(inhibitor)", fc="#fde8e8")
    _box(ax, (7.1, 2.2), 2.6, 1.4, "Dead polymer\n(no radical)", fc="#f0f0f0")

    ax.annotate("", xy=(3.65, 5.1), xytext=(2.95, 5.1), arrowprops=dict(arrowstyle="->", lw=1.4))
    ax.text(3.3, 5.45, "light", ha="center", fontsize=9)
    ax.annotate("", xy=(7.05, 5.1), xytext=(6.35, 5.1), arrowprops=dict(arrowstyle="->", lw=1.4))
    ax.text(6.7, 5.45, "initiation", ha="center", fontsize=9)
    ax.annotate("", xy=(8.4, 3.65), xytext=(8.4, 4.35), arrowprops=dict(arrowstyle="->", lw=1.4))
    ax.text(8.7, 4.0, "termination", ha="left", fontsize=9)
    ax.annotate("", xy=(2.95, 2.9), xytext=(3.65, 2.9), arrowprops=dict(arrowstyle="->", lw=1.4, color="#a33"))
    ax.text(3.3, 1.85, "quenches radicals", ha="center", fontsize=9, color="#a33")
    ax.annotate("", xy=(7.05, 2.9), xytext=(2.95, 4.55), arrowprops=dict(arrowstyle="->", lw=1.2))
    ax.text(4.6, 3.55, "propagation\n(conversion $p$)", ha="center", fontsize=9)

    ax.text(
        5.0,
        0.45,
        "Light does not polymerize PEGDA directly. It generates radicals that compete\n"
        "for acrylate addition versus oxygen scavenging and radical–radical termination.",
        ha="center",
        va="center",
        fontsize=9,
    )
    fig.savefig(OUT / "fig_background_species.pdf", bbox_inches="tight")
    fig.savefig(OUT / "fig_background_species.png", bbox_inches="tight")
    plt.close(fig)


def fig_reciprocity():
    """Qualitative cartoon of failed exposure reciprocity (not measured data)."""
    dose = np.linspace(0, 1, 200)
    # Schematic: higher intensity reaches a given conversion at different dose.
    p_low = 1 - np.exp(-2.2 * dose)
    p_high = 1 - np.exp(-3.4 * dose)
    p_high *= 0.92 / p_high[-1] * 0.95
    p_low *= 0.78 / p_low[-1]

    fig, ax = plt.subplots(figsize=(6.2, 4.0))
    ax.plot(dose, p_low, color="tab:blue", lw=2.2, label="Low intensity, long time")
    ax.plot(dose, p_high, color="tab:orange", lw=2.2, label="High intensity, short time")
    ax.set_xlabel("Nominal dose $\\propto I\\,t$ (normalized, schematic)")
    ax.set_ylabel("Acrylate conversion $p$")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.legend(frameon=False, loc="lower right")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_title("Equal dose need not mean equal conversion (schematic)")
    ax.text(
        0.02,
        0.98,
        "Cartoon only: shape inspired by bimolecular termination\n"
        "and oxygen inhibition, not a fit to a specific resin.",
        transform=ax.transAxes,
        va="top",
        fontsize=8,
        color="#444",
    )
    fig.savefig(OUT / "fig_background_reciprocity.pdf", bbox_inches="tight")
    fig.savefig(OUT / "fig_background_reciprocity.png", bbox_inches="tight")
    plt.close(fig)


def fig_two_velocities():
    fig, axes = plt.subplots(1, 2, figsize=(7.6, 3.6), constrained_layout=True)

    ax = axes[0]
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 6)
    ax.axis("off")
    ax.set_title("Pattern velocity $\\mathbf{v}_\\mathrm{pattern}$")
    ax.add_patch(Rectangle((0.5, 1.2), 9, 3.2, fill=False, lw=1.2, ec="#333"))
    ax.text(5, 0.6, "resin (stationary in this panel)", ha="center")
    ax.add_patch(Rectangle((2.2, 1.6), 2.4, 2.4, facecolor="#ffe08a", edgecolor="#b8860b", alpha=0.9))
    ax.annotate("", xy=(7.2, 2.8), xytext=(4.8, 2.8), arrowprops=dict(arrowstyle="->", lw=2, color="#b8860b"))
    ax.text(6.2, 3.3, "scanned mask\nor moving image", ha="center", fontsize=9)

    ax = axes[1]
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 6)
    ax.axis("off")
    ax.set_title("Material velocity $\\mathbf{u}$")
    ax.add_patch(Rectangle((0.5, 1.2), 9, 3.2, fill=False, lw=1.2, ec="#333"))
    ax.add_patch(Rectangle((3.6, 1.5), 2.8, 2.6, facecolor="#cfe8ff", edgecolor="#1f4e79", alpha=0.9))
    ax.text(5, 2.8, "fixed UV\nwindow", ha="center", fontsize=9)
    ax.annotate("", xy=(8.6, 2.8), xytext=(1.2, 2.8), arrowprops=dict(arrowstyle="->", lw=2, color="#1f4e79"))
    ax.text(5, 0.6, "resin flows through a stationary illumination zone", ha="center")

    fig.savefig(OUT / "fig_background_two_velocities.pdf", bbox_inches="tight")
    fig.savefig(OUT / "fig_background_two_velocities.png", bbox_inches="tight")
    plt.close(fig)


def fig_network_stack():
    fig, ax = plt.subplots(figsize=(6.8, 4.4))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 7)
    ax.axis("off")
    layers = [
        (5.2, "#d6eaf8", "4. Mechanics: $G \\propto \\eta\\,[\\mathrm{PEGDA}]_0\\,k_BT$\n   (elastically active chains, after gelation)"),
        (3.7, "#d5f5e3", "3. Network: gelation, loops, $\\eta$\n   (not every reacted acrylate is a load-bearing crosslink)"),
        (2.2, "#fdebd0", "2. Conversion $p = 1 - [M]/[M]_0$\n   (fraction of acrylates that have reacted)"),
        (0.7, "#f5b7b1", "1. Light history $I(\\mathbf{x},t)$ and species fields\n   (initiator, radicals, oxygen, monomer)"),
    ]
    for y, fc, text in layers:
        _box(ax, (0.6, y), 8.8, 1.25, text, fc=fc, ec="#333")
    ax.set_title("Do not collapse these layers into a single “dose”")
    fig.savefig(OUT / "fig_background_network_stack.pdf", bbox_inches="tight")
    fig.savefig(OUT / "fig_background_network_stack.png", bbox_inches="tight")
    plt.close(fig)


def main():
    fig_species()
    fig_reciprocity()
    fig_two_velocities()
    fig_network_stack()
    print(f"Wrote schematics to {OUT}")


if __name__ == "__main__":
    main()
