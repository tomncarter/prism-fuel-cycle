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


# ---- pin-design volume fractions ---------------------------------------------
from geometry.params import AssemblyGeom
from materials.tier1_homog import (
    LEGACY_RECIPE, PIN_DESIGNS, recipe_from_cell_fractions, resolve_recipe,
)


def test_cell_to_interior_round_trip():
    """Interior fractions x interior area + explicit duct/gap == the book's cell fractions."""
    g = AssemblyGeom()
    a = g.assembly_pitch_cm ** 2
    f_in = g.duct_inner_flat_to_flat_cm ** 2 / a
    f_duct = (g.duct_flat_to_flat_cm ** 2 - g.duct_inner_flat_to_flat_cm ** 2) / a
    f_gap = 1 - g.duct_flat_to_flat_cm ** 2 / a
    for d in PIN_DESIGNS.values():
        rc = recipe_from_cell_fractions(d["fuel"], d["coolant"], d["structure"], bond=d["bond"])
        assert rc.fuel * f_in == pytest.approx(d["fuel"])
        assert rc.sodium * f_in + f_gap == pytest.approx(d["coolant"] + d["bond"])
        assert rc.steel * f_in + f_duct == pytest.approx(d["structure"])


def test_resolve_recipe_options():
    assert resolve_recipe({"pin_design": "legacy"}) == LEGACY_RECIPE
    assert resolve_recipe({}) == resolve_recipe({"pin_design": "anl_sprism_driver"})
    rc = resolve_recipe({"volume_fractions": {"basis": "interior", "fuel": 0.3, "coolant": 0.4, "structure": 0.2}})
    assert (rc.fuel, rc.sodium, rc.steel) == (0.3, 0.4, 0.2)
    with pytest.raises(ValueError):
        resolve_recipe({"pin_design": "nope"})


def test_fuel_fraction_sets_hm_density():
    hm = lambda m: sum(d for n, d in m.get_nuclide_atom_densities().items() if n.startswith(("U", "Pu")))
    a = make_smeared_driver("a", 0.2, u_vector="depleted", recipe=resolve_recipe({"pin_design": "tc_d090"}))
    b = make_smeared_driver("b", 0.2, u_vector="depleted", recipe=resolve_recipe({"pin_design": "tc_d120"}))
    assert hm(b) / hm(a) == pytest.approx(0.376 / 0.279, rel=1e-6)


def test_anl_designs_fill_the_cell():
    """ANL fractions sum to 100% of the cell, so the interior has no void."""
    for k in ("anl_abr_driver", "anl_sprism_driver", "anl_hcr_driver", "anl_blanket"):
        rc = resolve_recipe({"pin_design": k})
        assert rc.fuel + rc.sodium + rc.steel == pytest.approx(1.0, abs=2e-3)


def test_lwr_sf_tru_vector():
    from materials.tier1_homog import PU_VECTORS, TRU_WT_LWR_SF
    assert sum(TRU_WT_LWR_SF.values()) == pytest.approx(100.0, abs=0.01)
    v = PU_VECTORS["lwr_sf_tru"]
    assert sum(v.values()) == pytest.approx(1.0)
    m = make_smeared_driver("unf", 0.2, pu_vector="lwr_sf_tru", u_vector="depleted")
    dens = m.get_nuclide_atom_densities()
    assert {"Np237", "Am241", "Cm244"} <= set(dens)
