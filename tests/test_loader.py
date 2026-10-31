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
