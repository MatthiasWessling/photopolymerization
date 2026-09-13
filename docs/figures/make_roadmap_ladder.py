"""Stage-ladder schematic for the campaign roadmap note."""

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

OUT = Path(__file__).resolve().parent
plt.rcParams.update({"font.size": 9, "pdf.fonttype": 42, "ps.fonttype": 42})

STAGES = [
    ("0", "Parameter\nbook"),
    ("1", "Local\nkinetics"),
    ("2", "Analytic\nlimits"),
    ("3", "Projection\nfield"),
    ("4", "Space &\noxygen"),
    ("5", "Plug\nflow"),
    ("6", "Network\n$\\eta,G$"),
]


def main():
    fig, ax = plt.subplots(figsize=(7.8, 2.8))
    ax.set_xlim(0, 14.2)
    ax.set_ylim(0, 4.2)
    ax.axis("off")
    ax.set_title("Campaign ladder (stop at each gate; do not skip)")

    for i, (n, label) in enumerate(STAGES):
        x = 0.35 + i * 2.0
        box = FancyBboxPatch(
            (x, 1.55),
            1.7,
            1.55,
            boxstyle="round,pad=0.03,rounding_size=0.08",
            facecolor="#d6eaf8" if i < 3 else "#d5f5e3" if i < 6 else "#fdebd0",
            edgecolor="#1f4e79",
            linewidth=1.2,
        )
        ax.add_patch(box)
        ax.text(x + 0.85, 2.7, f"Stage {n}", ha="center", va="center", fontweight="bold")
        ax.text(x + 0.85, 2.05, label, ha="center", va="center")
        if i < len(STAGES) - 1:
            ax.annotate(
                "",
                xy=(x + 1.95, 2.32),
                xytext=(x + 1.72, 2.32),
                arrowprops=dict(arrowstyle="->", lw=1.4, color="#333"),
            )
            ax.text(x + 1.85, 2.85, "gate", ha="center", fontsize=7, color="#a33")

    ax.text(
        7.1,
        0.55,
        "Parked after Stage 6: moving pattern, free-volume Dobson chemistry, heat, inverse design.",
        ha="center",
        fontsize=8,
        color="#444",
    )
    fig.savefig(OUT / "fig_roadmap_stage_ladder.pdf", bbox_inches="tight")
    fig.savefig(OUT / "fig_roadmap_stage_ladder.png", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
