"""Stage 0 provenance tests."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from photopolymerization.parameters import (
    kinetic_params_from_book,
    load_parameter_book,
)

CONFIGS = Path(__file__).resolve().parents[1] / "configs"
BOOK = CONFIGS / "parameter_book.montgomery_2022.json"


def test_book_loads():
    book = load_parameter_book(BOOK)
    assert book["chemistry_key"] == "montgomery_2022_pegda250_irgacure819"
    assert book["initiator_class"] == "type_I"


def test_rejects_bad_status(tmp_path: Path):
    book = json.loads(BOOK.read_text())
    book["parameters"]["beta"]["status"] = "guessed"
    bad = tmp_path / "bad.json"
    bad.write_text(json.dumps(book))
    with pytest.raises(ValueError, match="illegal status"):
        load_parameter_book(bad)


def test_kinetic_extract():
    params = kinetic_params_from_book(load_parameter_book(BOOK))
    assert params.C_M0 == pytest.approx(8880.0)
    assert params.beta == pytest.approx(2.70e-3)
    assert params.m == 2.0
