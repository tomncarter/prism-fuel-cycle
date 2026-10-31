"""Batch bookkeeping: staggered startup, discharge cadence, cooling enforcement (T2.4)."""
from collections import Counter
from types import SimpleNamespace

from fuel_cycle.batch_manager import BatchConfig, BatchManager
from fuel_cycle.state import CycleState, InCoreEntry, assembly_id


def _positions(n_inner=12, n_outer=9, n_blanket=8):
    out = []
    for t, n in (("inner_fuel", n_inner), ("outer_fuel", n_outer), ("radial_blanket", n_blanket)):
        for i in range(n):
            out.append(SimpleNamespace(q=i, r=-i, type=t, material_id=f"{t}_q{i:+d}_r{-i:+d}"))
    out.append(SimpleNamespace(q=99, r=0, type="reflector", material_id="reflector_q+99_r+0"))
    return out


FUEL_TYPES = frozenset({"inner_fuel", "outer_fuel", "radial_blanket"})


def _start(config=None):
    positions = _positions()
    mgr = BatchManager(config or BatchConfig.prism_breeder(blanket_cycles=4))
    state = CycleState.initialize(positions, "test", fuel_types=FUEL_TYPES)
    mgr.stagger_initial_core(state)
    return positions, mgr, state


def _reload(mgr, state, positions, cycle):
    fresh = [
        InCoreEntry(
            material_id=assembly_id(p.material_id, cycle + 1), q=p.q, r=p.r,
            assembly_type=p.type, cycle_loaded=cycle + 1, position_id=p.material_id,
        )
        for p in mgr.vacant_fuel_positions(state, positions, fuel_types=FUEL_TYPES)
    ]
    mgr.load_fresh(state, fresh)
    return fresh


def test_initialize_tracks_only_fuel_and_names_are_unique():
    positions, _, state = _start()
    assert len(state.in_core) == 29
    assert all(mid.endswith("_c000") for mid in state.in_core)
    assert all(e.position_id and e.position_id in mid for mid, e in state.in_core.items())


def test_stagger_spreads_ages_evenly():
    _, _, state = _start()
    ages = Counter((e.assembly_type, e.cycles_in_core) for e in state.in_core.values())
    assert [ages[("inner_fuel", a)] for a in range(3)] == [4, 4, 4]
    assert [ages[("outer_fuel", a)] for a in range(3)] == [3, 3, 3]
    assert [ages[("radial_blanket", a)] for a in range(4)] == [2, 2, 2, 2]


def test_one_batch_discharges_each_cycle():
    positions, mgr, state = _start()
    for cycle in range(6):
        _, discharged = mgr.step(state, cycle)
        by_type = Counter(e.assembly_type for e in discharged)
        assert by_type["inner_fuel"] == 4
        assert by_type["outer_fuel"] == 3
        assert by_type["radial_blanket"] == 2
        _reload(mgr, state, positions, cycle)
        assert len(state.in_core) == 29


def test_reload_fills_every_vacant_slot_with_a_new_name():
    positions, mgr, state = _start()
    mgr.step(state, 0)
    fresh = _reload(mgr, state, positions, 0)
    assert all(e.material_id.endswith("_c001") for e in fresh)
    assert mgr.vacant_fuel_positions(state, positions, fuel_types=FUEL_TYPES) == []
    slots = [e.position_id for e in state.in_core.values()]
    assert len(slots) == len(set(slots))


def test_cooling_enforced_for_one_full_cycle():
    """T2.4: nothing discharged after cycle c reaches the NFRC before the refueling after c+1."""
    positions, mgr, state = _start()
    discharged_at = {}
    for cycle in range(8):
        ready, discharged = mgr.step(state, cycle)
        for e in discharged:
            discharged_at[e.material_id] = cycle
        for e in ready:
            assert cycle - discharged_at[e.material_id] >= 1
            assert e.source_cycle == discharged_at[e.material_id]
        _reload(mgr, state, positions, cycle)
    # first NFRC feed arrives at the refueling after cycle 1
    _, mgr2, state2 = _start()
    ready0, _ = mgr2.step(state2, 0)
    assert ready0 == []


def test_driver_residence_is_three_cycles():
    positions, mgr, state = _start()
    loaded_at = {mid: 0 for mid in state.in_core}
    for cycle in range(9):
        _, discharged = mgr.step(state, cycle)
        for e in discharged:
            if e.material_id.endswith("_c000"):
                continue   # startup assemblies carry nominal ages
            if e.assembly_type != "radial_blanket":
                assert cycle + 1 - loaded_at[e.material_id] == 3
        for e in _reload(mgr, state, positions, cycle):
            loaded_at[e.material_id] = cycle + 1


def test_equilibrium_check_needs_a_full_window():
    _, mgr, state = _start()
    assert not mgr.is_at_equilibrium(state)
