# fuel_cycle/mass_balance.py
"""
Whole-system heavy-metal balance for a cycle-runner run directory (T4.3).

Every gram of HM that entered the system (initial core, external TRU feed,
makeup U) must be somewhere now:

    initial + external + makeup
        = core burnup (sum of BOC - EOC over cycles)
        + in core now (carried-over EOC + fresh assemblies awaiting their first cycle)
        + cooling queue + recovered U/TRU stocks
        + NFRC waste + out-of-core decay loss

Also checks that each assembly carried into the next cycle starts with the HM
it ended on, and that the cumulative fabricated mass matches the fresh
assemblies actually built.

Usage:
    python -m fuel_cycle.mass_balance runs/<run_dir> [loadout.json]
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import openmc.deplete

from fuel_cycle.cycle_runner import (
    CycleRunner, CycleRunnerConfig, _fresh_material, _hm_grams,
)
from fuel_cycle.state import CycleState

GATE_REL = 1e-3   # T4.3: 0.1% of HM


def _total(comp: dict) -> float:
    return sum(comp.values())


def check(run_dir: str, loadout_path: str = "core_loadouts/prism_breakeven.json") -> bool:
    run = Path(run_dir)
    state = CycleState.from_json(str(run / CycleRunner.STATE_FILE))
    runner = CycleRunner(loadout_path, run_dir, CycleRunnerConfig(cycle_length_days=1.0))
    runner._load_loadout()
    vol = runner._assembly_vol_cm3

    burnup_g = 0.0
    fresh_built_g = 0.0
    max_carry_err = 0.0
    prev_eoc: dict = {}
    eoc_hm: dict = {}
    print(f"{'cycle':>5} {'BOC HM kg':>11} {'EOC HM kg':>11} {'burned kg':>10} "
          f"{'fresh kg':>9} {'carry err':>10}")
    for c in range(state.current_cycle):
        cdir = run / f"cycle_{c:03d}"
        results = openmc.deplete.Results(str(cdir / "depletion_results.h5"))
        boc = {n: _total(_hm_grams(m)) for n, m in runner._export(str(cdir), results, 0).items()}
        eoc_mats = runner._export(str(cdir), results, -1)
        eoc_hm = {n: _hm_grams(m) for n, m in eoc_mats.items()}
        eoc = {n: _total(h) for n, h in eoc_hm.items()}

        fresh = [n for n in boc if n.endswith(f"_c{c:03d}") and c > 0]
        fresh_built_g += sum(boc[n] for n in fresh)
        for n in boc:
            if n in prev_eoc:
                max_carry_err = max(max_carry_err, abs(boc[n] - prev_eoc[n]) / prev_eoc[n])
        burned = sum(boc.values()) - sum(eoc.values())
        burnup_g += burned
        print(f"{c:5d} {sum(boc.values())/1e3:11.1f} {sum(eoc.values())/1e3:11.1f} "
              f"{burned/1e3:10.2f} {sum(boc[n] for n in fresh)/1e3:9.1f} {max_carry_err:10.2e}")
        prev_eoc = eoc

    # In core now: carried-over assemblies at last EOC + fresh ones not yet irradiated
    in_core_g = 0.0
    pending_g = 0.0
    for mid, e in state.in_core.items():
        if e.cycle_loaded == state.current_cycle:
            ref = runner._ref_mats[e.assembly_type]
            g = _total(_hm_grams(_fresh_material(mid, ref, e.fresh_hm_vec, vol)))
            pending_g += g
            in_core_g += g
        else:
            in_core_g += _total(eoc_hm[mid])
    fresh_built_g += pending_g

    L = state.ledger_g
    cooling_g = sum(_total(e.hm_g) for e in state.cooling_queue)
    stocks_g = _total(state.tru_stock_g) + _total(state.u_stock_g)
    inputs = L["initial_u"] + L["initial_tru"] + L.get("external_tru", 0.0) + L.get("makeup_u", 0.0)
    accounted = (burnup_g + in_core_g + cooling_g + stocks_g
                 + L.get("waste", 0.0) + L.get("decay_loss", 0.0))
    rel = abs(inputs - accounted) / inputs
    fab_rel = abs(fresh_built_g - L.get("fabricated", 0.0)) / max(L.get("fabricated", 1.0), 1.0)

    kg = lambda g: f"{g/1e3:12.2f} kg"
    print("\nInputs:   initial U", kg(L["initial_u"]), " initial TRU", kg(L["initial_tru"]))
    print("          external TRU", kg(L.get("external_tru", 0.0)),
          " makeup U", kg(L.get("makeup_u", 0.0)))
    print("Where:    burned", kg(burnup_g), " in core", kg(in_core_g),
          " cooling", kg(cooling_g))
    print("          stocks", kg(stocks_g), " waste", kg(L.get("waste", 0.0)),
          " decay loss", kg(L.get("decay_loss", 0.0)))
    print(f"\nSystem HM balance error   : {rel:.2e}  (gate {GATE_REL:.0e})")
    print(f"Fabricated vs built fresh : {fab_rel:.2e}")
    print(f"Max carry-over HM error   : {max_carry_err:.2e}")
    ok = rel < GATE_REL and fab_rel < GATE_REL and max_carry_err < GATE_REL
    print("T4.3", "PASS" if ok else "FAIL")
    return ok


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    sys.exit(0 if check(*sys.argv[1:3]) else 1)
