"""NFRC processing and fresh-fuel fabrication: mass balance (T4.3) and isotopics."""
import math

import pytest

from fuel_cycle.nfrc import (
    atomic_mass, decay, electrorefine, fabricate, mass_number, compute_tru_metrics,
)

RG_PU = {"Pu239": 0.55, "Pu240": 0.25, "Pu241": 0.12, "Pu242": 0.08}
DU = {"U235": 0.002, "U238": 0.998}


def _mol(comp):
    return sum(g / atomic_mass(n) for n, g in comp.items())


def _tru_frac(hm_vec):
    return sum(f for n, f in hm_vec.items() if not n.startswith("U"))


def test_atomic_mass_is_physical():
    assert atomic_mass("Pu239") == pytest.approx(239.052, abs=1e-3)
    assert atomic_mass("Am242_m1") == pytest.approx(242.06, abs=0.01)


def test_mass_number():
    assert mass_number("Pu239") == 239
    assert mass_number("Am242_m1") == 242
    assert mass_number("U235") == 235


def test_decay_pu241_half_life():
    out = decay({"Pu241": 100.0}, 14.29)
    assert out["Pu241"] == pytest.approx(50.0)
    assert out["Am241"] == pytest.approx(50.0, rel=1e-5)   # beta decay: same A, ~no mass loss


def test_decay_conserves_atoms_and_loses_alpha_mass():
    hm = {"Pu238": 100.0, "Pu241": 50.0, "Cm244": 10.0, "U238": 1000.0}
    out = decay(hm, 2.0)
    assert _mol(out) == pytest.approx(_mol(hm), rel=1e-12)
    assert sum(out.values()) < sum(hm.values())   # alpha decays shed He


def test_decay_zero_time_is_identity():
    hm = {"Pu241": 10.0, "Np239": 1.0}
    out = {n: g for n, g in decay(hm, 0.0).items() if g > 0.0}
    assert out == pytest.approx(hm)


def test_electrorefine_conserves_mass():
    hm = {"U238": 900.0, "U235": 2.0, "Pu239": 80.0, "Am241": 3.0, "Np237": 1.0}
    tru, u, waste = electrorefine(hm, tru_recovery=0.999, u_recovery=0.95)
    total = sum(tru.values()) + sum(u.values()) + sum(waste.values())
    assert total == pytest.approx(sum(hm.values()), rel=1e-12)
    assert set(tru) == {"Pu239", "Am241", "Np237"}
    assert sum(u.values()) == pytest.approx(0.95 * 902.0)


def test_startup_reload_is_all_external_rg_pu():
    """Before any discharge reaches the NFRC, TRU comes entirely from the external feed."""
    fab = fabricate(10, 100.0, 0.21, {}, {}, DU, RG_PU)
    assert fab.from_tru_stock == {}
    assert _mol(fab.from_external) == pytest.approx(0.21 * 1000.0)
    assert _mol(fab.from_makeup_u) == pytest.approx(0.79 * 1000.0)
    assert _tru_frac(fab.hm_vec) == pytest.approx(0.21)
    assert fab.hm_vec["Pu239"] / fab.hm_vec["Pu240"] == pytest.approx(0.55 / 0.25)


def test_recovered_tru_used_before_external():
    stock = {"Pu239": atomic_mass("Pu239") * 300.0, "Am241": atomic_mass("Am241") * 10.0}   # 310 mol
    fab = fabricate(10, 100.0, 0.21, stock, {}, DU, RG_PU)    # needs 210 mol TRU
    assert fab.from_external == {}
    assert _mol(fab.from_tru_stock) == pytest.approx(210.0)
    assert _mol(stock) == pytest.approx(100.0)                # breeder surplus stays in stock
    # pro-rata draw keeps the stock's isotopic ratio
    ratio = (fab.from_tru_stock["Am241"] / atomic_mass("Am241")) / (fab.from_tru_stock["Pu239"] / atomic_mass("Pu239"))
    assert ratio == pytest.approx(10 / 300)


def test_shortfall_topped_up_from_external():
    stock = {"Pu239": atomic_mass("Pu239") * 50.0}
    fab = fabricate(10, 100.0, 0.21, stock, {}, DU, RG_PU)
    assert _mol(fab.from_tru_stock) == pytest.approx(50.0)
    assert _mol(fab.from_external) == pytest.approx(160.0)
    assert _mol(stock) == pytest.approx(0.0, abs=1e-9)
    assert _tru_frac(fab.hm_vec) == pytest.approx(0.21)


def test_recycled_u_replaces_makeup():
    u_stock = {"U238": atomic_mass("U238") * 500.0, "U235": atomic_mass("U235") * 1.0}
    fab = fabricate(10, 100.0, 0.21, {}, u_stock, DU, RG_PU, recycle_u=True)
    assert _mol(fab.from_u_stock) == pytest.approx(501.0)
    assert _mol(fab.from_makeup_u) == pytest.approx(289.0)


def test_recycle_u_off_leaves_u_stock_untouched():
    u_stock = {"U238": 238.0 * 500.0}
    fab = fabricate(10, 100.0, 0.21, {}, u_stock, DU, RG_PU, recycle_u=False)
    assert fab.from_u_stock == {}
    assert u_stock == {"U238": 238.0 * 500.0}


def test_nfrc_to_fabrication_mass_balance():
    """Stocks before + NFRC products + inputs == stocks after + fabricated (T4.3)."""
    discharged = {"U238": 80_000.0, "U235": 150.0, "Pu239": 12_000.0,
                  "Pu240": 6_000.0, "Pu241": 1_000.0, "Am241": 300.0, "Cm244": 20.0}
    tru_stock = {"Pu239": 2_000.0}
    u_stock: dict = {}
    before = sum(tru_stock.values()) + sum(u_stock.values())

    hm = decay(discharged, 1.0)
    decay_loss = sum(discharged.values()) - sum(hm.values())
    tru, u, waste = electrorefine(hm)
    for stock, prod in ((tru_stock, tru), (u_stock, u)):
        for n, g in prod.items():
            stock[n] = stock.get(n, 0.0) + g

    fab = fabricate(30, 350.0, 0.21, tru_stock, u_stock, DU, RG_PU, recycle_u=True)
    after = sum(tru_stock.values()) + sum(u_stock.values())

    lhs = before + sum(discharged.values()) - decay_loss - sum(waste.values()) \
        + sum(fab.from_external.values()) + sum(fab.from_makeup_u.values())
    rhs = after + fab.total_g
    assert lhs == pytest.approx(rhs, rel=1e-12)


def test_compute_tru_metrics():
    tru_g, ma = compute_tru_metrics({"U238": 100.0, "Pu239": 30.0, "Am241": 10.0})
    assert tru_g == 40.0
    assert ma == pytest.approx(0.25)
    assert compute_tru_metrics({"U238": 1.0}) == (0.0, 0.0)
