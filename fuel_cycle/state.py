# fuel_cycle/state.py
"""
Persistent cycle state — dataclasses and JSON serialisation.

No OpenMC dependency.  The state is written to cycle_state.json in the run
directory at the end of every cycle so that interrupted runs can restart cleanly.

In-core compositions are NOT stored here; they live in per-cycle
depletion_results.h5 files.  This file records the bookkeeping needed to
drive the next cycle: which assembly is where, how old it is, the heavy metal
of each assembly in the cooling queue, and the recovered U/TRU stocks.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field, asdict
from typing import Dict, List, Optional


@dataclass
class InCoreEntry:
    """Tracks a single assembly currently in the active core."""
    material_id: str      # unique per assembly (slot id + load cycle); OpenMC material name
    q: int                # axial hex coordinate
    r: int
    assembly_type: str    # "inner_fuel", "outer_fuel", etc.
    cycle_loaded: int     # cycle index when this assembly was first loaded
    cycles_in_core: int = 0  # incremented each refueling
    position_id: str = ""    # lattice slot (AssemblyPosition.material_id) it occupies
    fresh_hm_vec: Dict[str, float] = field(default_factory=dict)
    # HM atom fractions as fabricated; used to build the material on its first cycle


@dataclass
class CoolingEntry:
    """Tracks a discharged assembly in its in-vessel cooling hold."""
    material_id: str
    assembly_type: str
    cycles_remaining: int  # cooling cycles left before NFRC processing
    source_cycle: int      # cycle index whose depletion_results.h5 holds the EOC composition
    hm_g: Dict[str, float] = field(default_factory=dict)  # EOC heavy metal at discharge


@dataclass
class CycleRecord:
    """Summary metrics for one completed irradiation cycle."""
    cycle: int
    boec_keff: float
    boec_keff_unc: float
    eoec_keff: float
    eoec_keff_unc: float
    conversion_ratio: Optional[float]    # time-averaged tally CR; None if tallies missing
    tru_mass_kg: float                   # total in-core TRU at BOC
    ma_tru_frac: float                   # MA / TRU weight fraction at BOC
    discharged_tru_kg: float             # TRU sent to NFRC this cycle
    fresh_tru_kg: float                  # TRU loaded in fresh assemblies this cycle
    fissile_boc_kg: float = 0.0          # in-core U235 + Pu239 + Pu241 at BOC
    fissile_eoc_kg: float = 0.0
    external_tru_kg: float = 0.0         # external feed used to fabricate this reload
    makeup_u_kg: float = 0.0             # fresh DU/NU used this reload (drivers + blankets)
    tru_stock_kg: float = 0.0            # recovered TRU left after fabrication (surplus)
    reload_scale: float = 1.0            # reload TRU fraction / startup TRU fraction
    reload_tru_frac: Dict[str, float] = field(default_factory=dict)  # driver type -> TRU/HM atom frac
    predicted_boec_keff: Optional[float] = None  # reload search's k for the next cycle


@dataclass
class CycleState:
    """
    Full mutable state of the fuel cycle simulation.

    Serialised to / from cycle_state.json.
    """
    loadout_name: str
    current_cycle: int
    in_core: Dict[str, InCoreEntry]      # material_id -> InCoreEntry
    cooling_queue: List[CoolingEntry]
    cycle_records: List[CycleRecord]
    tru_stock_g: Dict[str, float] = field(default_factory=dict)  # recovered TRU awaiting fab
    u_stock_g: Dict[str, float] = field(default_factory=dict)    # recovered U awaiting fab
    ledger_g: Dict[str, float] = field(default_factory=dict)
    # cumulative HM flows: external_tru, makeup_u, waste, decay_loss, nfrc_in, fabricated

    # ------------------------------------------------------------------ I/O --

    def to_json(self, path: str) -> None:
        with open(path, "w") as f:
            json.dump(asdict(self), f, indent=2)

    @classmethod
    def from_json(cls, path: str) -> CycleState:
        with open(path) as f:
            d = json.load(f)
        return cls(
            loadout_name=d["loadout_name"],
            current_cycle=d["current_cycle"],
            in_core={
                mid: InCoreEntry(**v)
                for mid, v in d["in_core"].items()
            },
            cooling_queue=[CoolingEntry(**v) for v in d["cooling_queue"]],
            cycle_records=[CycleRecord(**v) for v in d["cycle_records"]],
            tru_stock_g=d.get("tru_stock_g", {}),
            u_stock_g=d.get("u_stock_g", {}),
            ledger_g=d.get("ledger_g", {}),
        )

    def add_to_ledger(self, key: str, grams: float) -> None:
        self.ledger_g[key] = self.ledger_g.get(key, 0.0) + grams

    # ----------------------------------------------------------- factory -----

    @classmethod
    def initialize(
        cls,
        positions,          # list[AssemblyPosition] from core_maps.loader
        loadout_name: str,
        fuel_types: frozenset = frozenset({"inner_fuel", "outer_fuel"}),
    ) -> CycleState:
        """
        Build cycle 0 state with all fuel assemblies freshly loaded (age 0).
        Non-fuel positions (reflector, shield, etc.) are not tracked.
        """
        in_core = {}
        for p in positions:
            if p.type not in fuel_types:
                continue
            mid = assembly_id(p.material_id, 0)
            in_core[mid] = InCoreEntry(
                material_id=mid,
                q=p.q,
                r=p.r,
                assembly_type=p.type,
                cycle_loaded=0,
                cycles_in_core=0,
                position_id=p.material_id,
            )
        return cls(
            loadout_name=loadout_name,
            current_cycle=0,
            in_core=in_core,
            cooling_queue=[],
            cycle_records=[],
        )


def assembly_id(position_id: str, cycle_loaded: int) -> str:
    """Unique assembly name: the slot it is loaded into plus its load cycle."""
    return f"{position_id}_c{cycle_loaded:03d}"
