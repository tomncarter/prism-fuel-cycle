"""Fissile-loading search: root finder on synthetic k(x) and the loadout knobs."""
import random

import pytest

from fuel_cycle.search import (
    apply_loading_scale, loading_knob, scale_bounds, secant_search,
)


def _linear_core(k1=1.27, worth=0.55, sigma=0.0, seed=1):
    """k(x) with dk/k = worth dx/x around x=1; optional Gaussian noise."""
    rng = random.Random(seed)
    calls = []

    def evaluate(x):
        calls.append(x)
        k = k1 * (1.0 + worth * (x - 1.0))
        return k + rng.gauss(0.0, sigma), sigma
    return evaluate, calls


def test_converges_on_noiseless_core():
    evaluate, calls = _linear_core()
    res = secant_search(evaluate, target=1.005, tol_pcm=10)
    assert res.converged
    assert res.k == pytest.approx(1.005, abs=1e-4)
    assert len(calls) <= 4
    assert len(res.history) == len(calls)


def test_converges_with_monte_carlo_noise():
    evaluate, _ = _linear_core(sigma=0.0012, seed=7)
    res = secant_search(evaluate, target=1.005, tol_pcm=100)
    assert res.converged
    assert abs(res.k - 1.005) <= max(100e-5, 2 * 0.0012)


def test_accepts_first_point_inside_tolerance():
    evaluate, calls = _linear_core(k1=1.0055)
    res = secant_search(evaluate, target=1.005, tol_pcm=100)
    assert res.converged and calls == [1.0]


def test_reports_unconverged_when_pinned_at_bound():
    """A core that cannot reach the target inside the bounds is flagged, not faked."""
    evaluate, _ = _linear_core(k1=0.80, worth=0.1)
    res = secant_search(evaluate, target=1.005, bounds=(0.2, 1.5))
    assert not res.converged
    assert res.x == pytest.approx(1.5)


BREAKEVEN = {
    "inner_fuel": {"universe": "u_inner_fuel", "fuel_params": {"pu_hm_frac": 0.21, "zr_frac": 0.1}},
    "outer_fuel": {"universe": "u_outer_fuel", "fuel_params": {"pu_hm_frac": 0.21, "zr_frac": 0.1}},
    "radial_blanket": {"universe": "u_radial_blanket", "fuel_params": {"u235_frac": 0.002}},
    "reflector": {"universe": "u_reflector"},
}
UNF = {
    "inner_fuel": {"universe": "u_inner_fuel", "fuel_params": {"pu_hm_frac": 0.189}},
    "outer_fuel": {"universe": "u_outer_fuel", "fuel_params": {"pu_hm_frac": 0.227}},
}
WG_U = {
    "inner_fuel": {"universe": "u_inner_fuel",
                   "fuel_params": {"u_vector": {"weapons_grade": 0.224, "depleted": 0.776}}},
    "outer_fuel": {"universe": "u_outer_fuel",
                   "fuel_params": {"u_vector": {"weapons_grade": 0.224, "depleted": 0.776}}},
}
HALEU = {
    "inner_fuel": {"universe": "u_inner_fuel", "fuel_params": {"u_vector": "haleu"}},
}


def test_knob_selection():
    assert loading_knob(BREAKEVEN) == "tru"
    assert loading_knob(WG_U) == "u_fissile"
    assert loading_knob(HALEU) == "u_fissile"


def test_tru_scale_keeps_zone_ratio_and_leaves_blankets():
    out = apply_loading_scale(UNF, 0.9)
    assert out["inner_fuel"]["fuel_params"]["pu_hm_frac"] == pytest.approx(0.1701)
    assert out["outer_fuel"]["fuel_params"]["pu_hm_frac"] == pytest.approx(0.2043)
    blk = apply_loading_scale(BREAKEVEN, 0.9)
    assert blk["radial_blanket"] == BREAKEVEN["radial_blanket"]
    assert BREAKEVEN["inner_fuel"]["fuel_params"]["pu_hm_frac"] == 0.21   # input untouched


def test_u_blend_scale_backfills_with_depleted():
    out = apply_loading_scale(WG_U, 0.8)["inner_fuel"]["fuel_params"]["u_vector"]
    assert out["weapons_grade"] == pytest.approx(0.1792)
    assert out["depleted"] == pytest.approx(0.8208)


def test_named_u_vector_becomes_blend():
    out = apply_loading_scale(HALEU, 0.85)["inner_fuel"]["fuel_params"]["u_vector"]
    assert out == pytest.approx({"haleu": 0.85, "depleted": 0.15})
    with pytest.raises(ValueError):
        apply_loading_scale(HALEU, 1.2)


def test_bounds():
    assert scale_bounds(BREAKEVEN) == pytest.approx((0.2, 0.35 / 0.21))
    assert scale_bounds(UNF)[1] == pytest.approx(0.35 / 0.227)
    assert scale_bounds(HALEU)[1] == pytest.approx(1.0)
    assert scale_bounds(WG_U)[1] == pytest.approx(2.0)
