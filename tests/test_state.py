"""CycleState JSON persistence."""
import json

from fuel_cycle.state import CycleState, InCoreEntry, CoolingEntry, CycleRecord


def _state():
    return CycleState(
        loadout_name="prism_breakeven",
        current_cycle=4,
        in_core={
            "inner_fuel_q+0_r+0_c004": InCoreEntry(
                material_id="inner_fuel_q+0_r+0_c004", q=0, r=0,
                assembly_type="inner_fuel", cycle_loaded=4,
                position_id="inner_fuel_q+0_r+0",
                fresh_hm_vec={"U238": 0.79, "Pu239": 0.21},
            )
        },
        cooling_queue=[CoolingEntry(
            material_id="outer_fuel_q+5_r+0_c001", assembly_type="outer_fuel",
            cycles_remaining=1, source_cycle=3, hm_g={"Pu239": 1234.5},
        )],
        cycle_records=[CycleRecord(
            cycle=3, boec_keff=1.004, boec_keff_unc=1e-4, eoec_keff=1.0,
            eoec_keff_unc=1e-4, conversion_ratio=1.05, tru_mass_kg=3000.0,
            ma_tru_frac=0.01, discharged_tru_kg=900.0, fresh_tru_kg=950.0,
            external_tru_kg=50.0,
        )],
        tru_stock_g={"Pu239": 10.0},
        u_stock_g={"U238": 20.0},
        ledger_g={"external_tru": 50_000.0},
    )


def test_round_trip(tmp_path):
    path = tmp_path / "cycle_state.json"
    s = _state()
    s.to_json(str(path))
    assert CycleState.from_json(str(path)) == s


def test_loads_state_written_before_stocks_and_compositions(tmp_path):
    """cycle_state.json files from the earlier schema still load."""
    path = tmp_path / "old.json"
    path.write_text(json.dumps({
        "loadout_name": "x", "current_cycle": 1,
        "in_core": {"a": {"material_id": "a", "q": 0, "r": 0,
                          "assembly_type": "inner_fuel", "cycle_loaded": 0,
                          "cycles_in_core": 1}},
        "cooling_queue": [{"material_id": "b", "assembly_type": "outer_fuel",
                           "cycles_remaining": 1, "source_cycle": 0}],
        "cycle_records": [{"cycle": 0, "boec_keff": 1.0, "boec_keff_unc": 0.0,
                           "eoec_keff": 1.0, "eoec_keff_unc": 0.0,
                           "conversion_ratio": None, "tru_mass_kg": 0.0,
                           "ma_tru_frac": 0.0, "discharged_tru_kg": 0.0,
                           "fresh_tru_kg": 0.0}],
    }))
    s = CycleState.from_json(str(path))
    assert s.tru_stock_g == {}
    assert s.cooling_queue[0].hm_g == {}
    assert s.in_core["a"].position_id == ""


def test_ledger_accumulates():
    s = _state()
    s.add_to_ledger("external_tru", 5.0)
    s.add_to_ledger("waste", 1.0)
    assert s.ledger_g == {"external_tru": 50_005.0, "waste": 1.0}
