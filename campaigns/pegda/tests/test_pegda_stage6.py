"""Stage 6 Zhu network tests (no Type-II bleaching)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from numpy.testing import assert_allclose

from photopolymerization.network.macosko import (
    G_ideal,
    G_shear,
    can_gel,
    eta,
    eta_noloop,
    network_params_from_book,
    p_gel,
    theta,
)
from photopolymerization.parameters import load_parameter_book

BOOK = (
    Path(__file__).resolve().parents[1] / "configs" / "parameter_book.zhu_2020.json"
)


def _net():
    return network_params_from_book(load_parameter_book(BOOK))


def test_pgel_identity():
    net = _net()
    assert_allclose(p_gel(20.0, net), net.C0_volpct / 20.0)


def test_never_gels_at_C0():
    net = _net()
    assert not can_gel(net.C0_volpct, net)
    assert float(eta(1.0, net.C0_volpct, net)) == 0.0


def test_noloop_p_squared():
    p = np.array([0.0, 0.3, 1.0])
    assert_allclose(eta_noloop(p), p * p)


def test_theta_unity_pgel_vanishes():
    net = _net()
    assert p_gel(1.0e9, net) < 1e-8
    assert_allclose(theta(1.0e9, net), 1.0, rtol=1e-8)


def test_loops_below_ideal_modulus():
    net = _net()
    et = float(eta(1.0, 20.0, net, loops=True))
    assert et < 1.0
    assert G_shear(et, 20.0, net) < G_ideal(20.0, net)
