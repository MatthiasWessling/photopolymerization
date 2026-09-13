#!/usr/bin/env python3
"""Write the Stage 0 provenance table and status figure."""

from __future__ import annotations

import csv
import json
import sys
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))

from photopolymerization.parameters import load_parameter_book, still_required_names

from _paths import FIGURES, config_path, output_dir


def main() -> None:
    cfg = json.loads(config_path("config.pegda_stage0.json").read_text())
    book_path = config_path(cfg["parameter_book"])
    book = load_parameter_book(book_path)
    out = output_dir("pegda_stage0")

    rows = []
    for name, entry in book["parameters"].items():
        rows.append(
            {
                "name": name,
                "value": entry["value"],
                "units": entry["units"],
                "status": entry["status"],
                "source": entry["source"],
            }
        )
    csv_path = out / "provenance_table.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    counts = Counter(r["status"] for r in rows)
    fig, ax = plt.subplots(figsize=(6.2, 3.6), constrained_layout=True)
    labels = ["transcribed", "calculated", "assumed", "still_required"]
    values = [counts.get(k, 0) for k in labels]
    colors = ["#1f4e79", "#2e86ab", "#f4a261", "#c0392b"]
    ax.bar(labels, values, color=colors)
    ax.set_ylabel("Number of parameters")
    ax.set_title(book["chemistry_key"])
    fig.savefig(out / "pegda_stage0_overview.pdf")
    fig.savefig(FIGURES / "pegda_stage0_overview.pdf")
    fig.savefig(FIGURES / "pegda_stage0_overview.png")
    plt.close(fig)

    summary = {
        "chemistry_key": book["chemistry_key"],
        "n_parameters": len(rows),
        "counts": dict(counts),
        "still_required": still_required_names(book),
        "csv": str(csv_path),
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
