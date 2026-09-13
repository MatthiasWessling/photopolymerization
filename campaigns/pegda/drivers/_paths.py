"""Resolve campaign config and output paths."""

from __future__ import annotations

from pathlib import Path

CAMPAIGN_ROOT = Path(__file__).resolve().parents[1]
CONFIGS = CAMPAIGN_ROOT / "configs"
REPO_ROOT = CAMPAIGN_ROOT.parents[1]
OUTPUTS = REPO_ROOT / "outputs"
FIGURES = REPO_ROOT / "docs" / "figures"


def config_path(name: str) -> Path:
    return CONFIGS / name


def output_dir(stem: str) -> Path:
    path = OUTPUTS / f"{stem}_out"
    path.mkdir(parents=True, exist_ok=True)
    return path
