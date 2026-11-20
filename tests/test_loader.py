"""Tests for core_maps.loader axial geometry."""
import json

from core_maps.loader import load_axial
from geometry.params import Prism3D


def _write(tmp_path, doc):
    path = tmp_path / "loadout.json"
    path.write_text(json.dumps(doc))
    return str(path)


def test_loadout_height_overrides_default(tmp_path):
    assert load_axial(_write(tmp_path, {"active_height_cm": 66.0})).active_height_cm == 66.0


def test_missing_height_falls_back_to_default(tmp_path):
    assert load_axial(_write(tmp_path, {})) == Prism3D()


# ---- ring + override layouts -------------------------------------------------
import glob

import pytest

from core_maps.loader import load_core, resolve_layout, rotate60, type_counts


def test_rotate60_cycles_ring_corners():
    seen = [rotate60(0, 3, k) for k in range(6)]
    assert seen == [(0, 3), (-3, 3), (-3, 0), (0, -3), (3, -3), (3, 0)]
    assert rotate60(0, 3, 6) == (0, 3)


def test_sym_expands_overrides():
    layout = {"rings": [{"range": [0, 3], "type": "fuel"}],
              "overrides": [{"q": 0, "r": 3, "sym": 3, "type": "rod"},
                            {"q": -3, "r": 3, "sym": 3, "type": "uss"}]}
    t = resolve_layout(layout)
    assert {p for p, v in t.items() if v == "rod"} == {(0, 3), (-3, 0), (3, -3)}
    assert {p for p, v in t.items() if v == "uss"} == {(-3, 3), (0, -3), (3, 0)}


@pytest.mark.parametrize("overrides, msg", [
    ([{"q": 0, "r": 5, "type": "rod"}], "outside"),
    ([{"q": 0, "r": 3, "type": "rod"}, {"q": 0, "r": 3, "type": "gem"}], "twice"),
    ([{"q": 0, "r": 3, "sym": 4, "type": "rod"}], "sym"),
])
def test_bad_overrides_raise(overrides, msg):
    layout = {"rings": [{"range": [0, 3], "type": "fuel"}], "overrides": overrides}
    with pytest.raises(ValueError, match=msg):
        resolve_layout(layout)


def test_ring_gap_raises():
    with pytest.raises(ValueError, match="no type"):
        resolve_layout({"rings": [{"range": [0, 1], "type": "a"}, {"range": [3, 3], "type": "b"}]})


def test_expected_counts_enforced(tmp_path):
    doc = {"layout": {"mode": "rings", "rings": [{"range": [0, 1], "type": "fuel"}]},
           "assembly_types": {"fuel": {"universe": "u_inner_fuel"}},
           "expected_counts": {"fuel": 6}}
    with pytest.raises(ValueError, match="expected_counts"):
        load_core(_write(tmp_path, doc))
    doc["expected_counts"] = {"fuel": 7}
    assert type_counts(load_core(_write(tmp_path, doc))[0]) == {"fuel": 7}


@pytest.mark.parametrize("path", sorted(glob.glob("core_loadouts/*.json")))
def test_every_loadout_loads(path):
    load_core(path)   # also checks expected_counts where a loadout declares them
