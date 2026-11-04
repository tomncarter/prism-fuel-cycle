"""Operating temperatures reach every material, including fuel-cycle copies."""
import openmc
import pytest

from build_prism_core import _require_temperatures, build_model
from fuel_cycle.cycle_runner import _clone_material, _fresh_material
from materials.tier1_homog import (
    T_COOLANT_K, T_INLET_K, T_INTERIOR_K, T_OUTLET_K,
    make_material_from_composition, make_smeared_driver, sodium_density,
)

LOADOUT = "core_loadouts/prism_breakeven_k103.json"


def test_sodium_density_at_core_conditions():
    assert sodium_density(T_INLET_K) == pytest.approx(0.867, abs=1e-3)
    assert sodium_density(T_COOLANT_K) == pytest.approx(0.851, abs=1e-3)
    assert sodium_density(T_OUTLET_K) == pytest.approx(0.835, abs=1e-3)


def test_every_model_material_is_hot():
    model = build_model(LOADOUT, "smoke")
    assert model.settings.temperature["method"] == "interpolation"
    temps = {m.name: m.temperature for m in model.materials}
    assert None not in temps.values()
    assert min(temps.values()) >= T_INLET_K
    assert temps["fuel_driver_inner_homog"] == T_INTERIOR_K
    assert temps["upper_plenum"] == T_OUTLET_K
    assert temps["lower_plenum"] == T_INLET_K


def test_builder_guard_refuses_missing_temperature():
    with pytest.raises(ValueError, match="cold"):
        _require_temperatures([openmc.Material(name="cold")])


def test_fuel_cycle_copies_keep_temperature():
    ref = make_smeared_driver("ref", 0.15, u_vector="depleted")
    assert _clone_material(ref, "clone", 1.0).temperature == T_INTERIOR_K
    fresh = _fresh_material("fresh", ref, {"U238": 0.8, "Pu239": 0.2}, 1.0)
    assert fresh.temperature == T_INTERIOR_K
    assert make_material_from_composition("c", {"U238": 1.0}).temperature == T_INTERIOR_K


def test_temperature_survives_xml_round_trip(tmp_path):
    """Assemblies carried between cycles are rebuilt from materials.xml."""
    path = tmp_path / "materials.xml"
    openmc.Materials([make_smeared_driver("m", 0.15, u_vector="depleted")]).export_to_xml(path)
    assert openmc.Materials.from_xml(path)[0].temperature == T_INTERIOR_K
