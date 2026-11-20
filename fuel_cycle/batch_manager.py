# fuel_cycle/batch_manager.py
"""
Batch management: discharge, in-vessel cooling, NFRC feed, fresh reload.

No OpenMC dependency — pure bookkeeping on CycleState.

Batch scheme
------------
Each fuel zone (inner_fuel, outer_fuel) operates independently with a
configurable maximum dwell time (max_cycles_in_core).  On each refueling:

  1. Age all in-core assemblies by one cycle.
  2. Discharge assemblies whose cycles_in_core == max_cycles_in_core.
  3. Add discharged assemblies to the cooling queue.
  4. Age the cooling queue (decrement cycles_remaining).
  5. Release assemblies with cycles_remaining == 0 → ready for NFRC.

The caller (cycle_runner) is responsible for:
  - Extracting EOC compositions for discharged assemblies.
  - Processing ready_for_nfrc assemblies through nfrc.electrorefine().
  - Calling load_fresh() with newly fabricated InCoreEntry objects.

Spatial shuffling note
----------------------
For Tier-1 (homogenised, per-zone shared materials) the physical position of
an assembly within a zone does not affect the neutronics.  The batch manager
therefore uses a simple FIFO discharge policy: the oldest assemblies in each
zone are discharged first, and fresh assemblies are loaded at the positions
vacated by the discharged batch.

For Tier-2 (per-pin, spatially resolved), a true in-to-out shuffle would move
assemblies from outer positions to inner positions across zones over successive
cycles.  That extension is implemented by subclassing or replacing this module.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List

from fuel_cycle.state import CycleState, InCoreEntry, CoolingEntry


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

@dataclass
class ZoneConfig:
    """Irradiation and cooling parameters for one assembly type."""
    max_cycles_in_core: int = 3  # irradiation dwell time in cycles
    cooling_cycles: int = 1      # in-vessel cooling hold before NFRC


@dataclass
class BatchConfig:
    """Batch configuration keyed by assembly type name."""
    zones: Dict[str, ZoneConfig] = field(default_factory=dict)
    default: ZoneConfig = field(default_factory=ZoneConfig)

    def for_type(self, asm_type: str) -> ZoneConfig:
        return self.zones.get(asm_type, self.default)

    @classmethod
    def prism_unf_recycle(cls) -> BatchConfig:
        """Standard PRISM UNF Recycle batch config (3-batch, 1-cycle cooling)."""
        return cls(zones={
            "inner_fuel": ZoneConfig(max_cycles_in_core=3, cooling_cycles=1),
            "outer_fuel": ZoneConfig(max_cycles_in_core=3, cooling_cycles=1),
        })

    @classmethod
    def prism_breeder(cls, blanket_cycles: int = 4) -> BatchConfig:
        """3-batch drivers plus internal/radial blankets with their own residence."""
        cfg = cls.prism_unf_recycle()
        for t in ("internal_blanket", "radial_blanket"):
            cfg.zones[t] = ZoneConfig(max_cycles_in_core=blanket_cycles, cooling_cycles=1)
        return cfg


# ---------------------------------------------------------------------------
# Manager
# ---------------------------------------------------------------------------

class BatchManager:
    """
    Manages the multi-batch fuel cycle bookkeeping for one PRISM core.

    All methods mutate the CycleState passed to them in-place and also
    return the relevant subsets for the caller to act on.
    """

    def __init__(self, config: BatchConfig):
        self.config = config

    # --------------------------------------------------------------- startup

    def stagger_initial_core(self, state: CycleState) -> None:
        """
        Give the all-fresh startup core a spread of nominal ages so that
        1/max_cycles_in_core of each zone discharges at every refueling,
        instead of the whole core leaving together.  Ages are dealt round-robin
        in slot order within each zone, which spreads each batch across the zone.
        """
        counters: dict[str, int] = {}
        for entry in sorted(state.in_core.values(), key=lambda e: e.position_id):
            n = self.config.for_type(entry.assembly_type).max_cycles_in_core
            i = counters.get(entry.assembly_type, 0)
            entry.cycles_in_core = i % n
            counters[entry.assembly_type] = i + 1

    # ------------------------------------------------------------------ step

    def step(
        self,
        state: CycleState,
        completed_cycle: int,
    ) -> tuple[list[CoolingEntry], list[CoolingEntry]]:
        """
        Advance state by one refueling event after cycle `completed_cycle`.

        Mutates state in-place.

        Returns:
            ready_for_nfrc:    assemblies that have finished cooling and are
                               ready for electrorefining this cycle.
            newly_discharged:  assemblies just removed from the core and
                               entered into the cooling queue.
        """
        # 1. Age all in-core assemblies
        for entry in state.in_core.values():
            entry.cycles_in_core += 1

        # 2. Discharge assemblies that have reached their dwell limit
        newly_discharged: list[CoolingEntry] = []
        remaining: dict[str, InCoreEntry] = {}
        for mid, entry in state.in_core.items():
            zone = self.config.for_type(entry.assembly_type)
            if entry.cycles_in_core >= zone.max_cycles_in_core:
                newly_discharged.append(CoolingEntry(
                    material_id=mid,
                    assembly_type=entry.assembly_type,
                    cycles_remaining=zone.cooling_cycles,
                    source_cycle=completed_cycle,
                ))
            else:
                remaining[mid] = entry

        state.in_core = remaining

        # 3 & 4. Advance the existing cooling queue BEFORE adding the newly
        # discharged assemblies — assemblies discharged this cycle do not start
        # their cooling countdown until the next refueling.
        ready_for_nfrc: list[CoolingEntry] = []
        still_cooling: list[CoolingEntry] = []
        for entry in state.cooling_queue:
            entry.cycles_remaining -= 1
            if entry.cycles_remaining <= 0:
                ready_for_nfrc.append(entry)
            else:
                still_cooling.append(entry)

        state.cooling_queue = still_cooling
        state.cooling_queue.extend(newly_discharged)

        return ready_for_nfrc, newly_discharged

    # --------------------------------------------------------------- reload

    def load_fresh(
        self,
        state: CycleState,
        fresh_entries: list[InCoreEntry],
    ) -> None:
        """
        Register freshly fabricated assemblies as in-core.
        Called by cycle_runner after NFRC processing and material construction.
        """
        for entry in fresh_entries:
            state.in_core[entry.material_id] = entry

    def vacant_fuel_positions(
        self,
        state: CycleState,
        all_positions,  # list[AssemblyPosition]
        fuel_types: frozenset = frozenset({"inner_fuel", "outer_fuel"}),
    ) -> list:
        """
        Return AssemblyPositions that are fuel slots not currently occupied
        in state.in_core.  These are the slots to fill with fresh assemblies.
        """
        occupied_ids = {e.position_id for e in state.in_core.values()}
        return [
            p for p in all_positions
            if p.type in fuel_types and p.material_id not in occupied_ids
        ]

    # ------------------------------------------------------ equilibrium check

    def is_at_equilibrium(
        self,
        state: CycleState,
        window: int = 3,
        keff_tolerance_pcm: float = 100.0,
        cr_tolerance: float = 0.01,
    ) -> bool:
        """
        Return True if the last `window` cycles show stable BOEC k-eff and CR.

        Stability criteria (T7.1 / T7.3):
          - BOEC k-eff variation < keff_tolerance_pcm over the window
          - CR variation < cr_tolerance (relative) over the window
        """
        records = state.cycle_records
        if len(records) < window:
            return False

        recent = records[-window:]
        keffs = [r.boec_keff for r in recent]
        keff_swing = (max(keffs) - min(keffs)) * 1e5  # pcm
        if keff_swing >= keff_tolerance_pcm:
            return False

        crs = [r.conversion_ratio for r in recent if r.conversion_ratio is not None]
        if len(crs) == window:
            cr_range = (max(crs) - min(crs)) / (sum(crs) / len(crs))
            if cr_range >= cr_tolerance:
                return False

        return True
