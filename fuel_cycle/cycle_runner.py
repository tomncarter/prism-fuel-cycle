# fuel_cycle/cycle_runner.py
"""
Cycle-managed fuel cycle loop for PRISM.

Orchestrates:
  - Per-cycle OpenMC model construction (one material per assembly)
  - Transport + depletion execution
  - BOC/EOC composition extraction from depletion_results.h5
  - Batch management (discharge, cooling, NFRC, fresh reload)
  - Recovered U/TRU stocks and external makeup (reactor-grade Pu by default)
  - CycleState persistence between cycles
  - Equilibrium detection

Assembly naming: each assembly gets a unique material name, the lattice slot it
occupies plus its load cycle (e.g. "inner_fuel_q+1_r-2_c004").  Slots are
reused; names are not, so a fresh assembly never inherits the composition of
the one it replaced.

Directory layout produced:
  {run_dir}/
    cycle_state.json         ← mutable state ledger (updated each cycle)
    cycle_000/
      materials.xml
      geometry.xml
      settings.xml
      tallies.xml
      openmc_simulation_n{i}.h5
      depletion_results.h5
    cycle_001/
      ...
"""
from __future__ import annotations

import json
import math
import os
import shutil
import warnings
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Optional

import openmc
import openmc.deplete

from core_maps.loader import load_axial, load_core, AssemblyPosition
from build_prism_core import build_model_per_assembly, RUN_PROFILES, CHAIN_FILE, POWER_MWth
from geometry.params import AssemblyGeom
from materials.tier1_homog import make_materials_from_loadout, PU_VECTORS
from fuel_cycle.state import CycleState, InCoreEntry, CycleRecord, assembly_id
from fuel_cycle.nfrc import (
    NFRCConfig, HEAVY_METAL_ELEMENTS, TRU_ELEMENTS,
    nuclide_element, electrorefine, decay, fabricate, compute_tru_metrics,
)
from fuel_cycle.batch_manager import BatchConfig, BatchManager
from fuel_cycle.search import DRIVER_UNIVERSES, secant_search


# Assembly universes that carry depletable heavy metal -> zone-default material key
BLANKET_UNIVERSES = {
    "u_internal_blanket": "fuel_blanket_internal_homog",
    "u_radial_blanket":   "fuel_blanket_radial_homog",
}

FISSILE_NUCLIDES = ("U235", "Pu239", "Pu241")
AVOGADRO = 6.02214076e23


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

@dataclass
class CycleRunnerConfig:
    """All tunable parameters for a cycle runner session."""
    # Depletion
    chain_file: str = CHAIN_FILE
    power_MWth: float = POWER_MWth
    cycle_length_days: Optional[float] = None
    # irradiation days per cycle; None = the loadout's "cycle_length_days"
    # (Triplett 2012: 23 months breakeven/breeder, 12 months UNF recycle), else 365
    first_steps_days: tuple = (5.0, 25.0) # short leading steps let Xe-135 / Sm-149 saturate
    max_step_days: float = 60.0           # remaining days split into steps no longer than this

    # OpenMC transport statistics: a build_prism_core.RUN_PROFILES name
    profile: str = "standard"

    # Batch management (blanket zones only apply when the loadout has blankets)
    batch: BatchConfig = field(default_factory=BatchConfig.prism_breeder)

    # NFRC (external feed defaults to reactor-grade Pu)
    nfrc: NFRCConfig = field(default_factory=NFRCConfig)

    # Reload fissile search: each reload's driver TRU fraction is set so the
    # next cycle's BOEC k-eff hits k_target (T1.1); above k_max the run halts (T7.4)
    search_reload: bool = True
    k_target: float = 1.005
    k_max: float = 1.010
    search_profile: str = "quick"
    search_tol_pcm: float = 100.0
    search_max_iter: int = 8

    # Convergence (T7.1 / T7.3)
    equilibrium_window: int = 3            # cycles to check for stability
    keff_tolerance_pcm: float = 100.0
    cr_tolerance: float = 0.01

    @property
    def timesteps(self) -> list[float]:
        """Depletion steps (days) covering one cycle."""
        lead = [d for d in self.first_steps_days]
        rest = self.cycle_length_days - sum(lead)
        if rest <= 0.0:
            return [self.cycle_length_days]
        n = math.ceil(rest / self.max_step_days)
        return lead + [rest / n] * n

    @property
    def power_W(self) -> float:
        return self.power_MWth * 1e6


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _clone_material(src: openmc.Material, name: str, volume: float) -> openmc.Material:
    """Copy of src under a new name and id, with the per-assembly volume set."""
    m = openmc.Material(name=name)
    for nuc in src.nuclides:
        m.add_nuclide(nuc.name, nuc.percent, nuc.percent_type)
    m.set_density(src.density_units, src.density)
    m.volume = volume
    m.depletable = True
    return m


def _fresh_material(
    name: str, ref: openmc.Material, hm_vec: dict, volume: float
) -> openmc.Material:
    """
    Fresh assembly with the zone reference's Zr/Na/steel atom densities and the
    same total HM atom density, filled with the fabricated HM vector.  Keeping
    HM atoms fixed makes the fabrication mass ledger exact.
    """
    dens = ref.get_nuclide_atom_densities()   # atom/b-cm
    hm_total = sum(d for n, d in dens.items() if _is_hm(n))
    m = openmc.Material(name=name)
    for nuc, d in dens.items():
        if not _is_hm(nuc):
            m.add_nuclide(nuc, d)
    for nuc, frac in hm_vec.items():
        if frac > 0.0:
            m.add_nuclide(nuc, frac * hm_total)
    m.set_density("sum")
    m.volume = volume
    m.depletable = True
    return m


def _is_hm(nuc: str) -> bool:
    return nuclide_element(nuc) in HEAVY_METAL_ELEMENTS


def _hm_grams(mat: openmc.Material) -> dict:
    """{nuclide: grams} of heavy metal in a material with its volume set."""
    return {n: mat.get_mass(n) for n in mat.get_nuclides() if _is_hm(n)}


def _hm_vec(mat: openmc.Material) -> dict:
    """HM atom fractions of a material."""
    dens = {n: d for n, d in mat.get_nuclide_atom_densities().items() if _is_hm(n)}
    total = sum(dens.values())
    return {n: d / total for n, d in dens.items()}


def _hm_mol(mat: openmc.Material) -> float:
    """Moles of heavy metal in a material with its volume set."""
    dens = mat.get_nuclide_atom_densities()
    return sum(d for n, d in dens.items() if _is_hm(n)) * 1e24 * mat.volume / AVOGADRO


def _sum_comps(comps) -> dict:
    total: dict = {}
    for comp in comps:
        for n, g in comp.items():
            total[n] = total.get(n, 0.0) + g
    return total


def _tru_g(comp: dict) -> float:
    return sum(g for n, g in comp.items() if nuclide_element(n) in TRU_ELEMENTS)


def _fissile_g(comp: dict) -> float:
    return sum(comp.get(n, 0.0) for n in FISSILE_NUCLIDES)


def _read_keff(results: openmc.deplete.Results) -> tuple[float, float, float, float]:
    """Return (boec_keff, boec_unc, eoec_keff, eoec_unc) from depletion results."""
    _, eigenvalues = results.get_keff(time_units="d")
    boec_k, boec_u = float(eigenvalues[0, 0]),  float(eigenvalues[0, 1])
    eoec_k, eoec_u = float(eigenvalues[-1, 0]), float(eigenvalues[-1, 1])
    return boec_k, boec_u, eoec_k, eoec_u


def _compute_cr(cycle_dir: str, timesteps: list[float]) -> Optional[float]:
    """
    Cycle conversion ratio from the cr_fissile_production / cr_fissile_destruction
    tallies, time-weighted over the beginning-of-step transport solves
    (openmc_simulation_n{i}.h5, i = 0 .. n_steps-1).

    Returns None if any statepoint or tally is missing.
    """
    prod = dest = 0.0
    for i, dt in enumerate(timesteps):
        sp_path = os.path.join(cycle_dir, f"openmc_simulation_n{i}.h5")
        if not os.path.exists(sp_path):
            return None
        try:
            with openmc.StatePoint(sp_path) as sp:
                prod += dt * float(sp.get_tally(name="cr_fissile_production").mean.sum())
                dest += dt * float(sp.get_tally(name="cr_fissile_destruction").mean.sum())
        except (KeyError, LookupError):
            return None
    return prod / dest if dest > 0.0 else None


# ---------------------------------------------------------------------------
# Main runner
# ---------------------------------------------------------------------------

class CycleRunner:
    """
    Manages a multi-cycle PRISM fuel cycle simulation.

    Usage:
        runner = CycleRunner(
            loadout_path="core_loadouts/prism_breakeven.json",
            run_dir="runs/prism_breakeven__cycles",
        )
        runner.run(n_cycles=20)
    """

    STATE_FILE = "cycle_state.json"

    def __init__(
        self,
        loadout_path: str,
        run_dir: str,
        config: Optional[CycleRunnerConfig] = None,
    ):
        self.loadout_path = loadout_path
        self.run_dir = Path(run_dir)
        self.cfg = config or CycleRunnerConfig()
        self.batch_mgr = BatchManager(self.cfg.batch)

        # Loaded lazily in run()
        self._positions: list[AssemblyPosition] = []
        self._asm_types: dict = {}
        self._base_mats_dict: dict = {}
        self._assembly_vol_cm3: float = 0.0
        self._ref_mats: dict = {}        # assembly type -> zone-default material (volume set)
        self._eoc_cache: dict = {}       # cycle -> {material name: EOC openmc.Material}
        self._driver_types: frozenset = frozenset()
        self._blanket_types: frozenset = frozenset()

    # ------------------------------------------------------------------ public

    def run(self, n_cycles: int, stop_at_equilibrium: bool = True) -> CycleState:
        """
        Run up to n_cycles irradiation cycles, saving state after each.

        If stop_at_equilibrium is True, stops early when the batch manager
        reports that the equilibrium criteria are satisfied.

        Returns the final CycleState.
        """
        self.run_dir.mkdir(parents=True, exist_ok=True)
        self._load_loadout()

        state = self._load_or_initialize_state()

        for _ in range(n_cycles):
            cycle = state.current_cycle
            print(f"\n{'='*60}")
            print(f"  Cycle {cycle:03d}  ({self.cfg.cycle_length_days:.0f} d  "
                  f"@ {self.cfg.power_MWth} MWth)")
            print(f"{'='*60}")

            record = self._run_one_cycle(state, cycle)
            state.cycle_records.append(record)
            state.current_cycle += 1
            state.to_json(self._state_path)

            print(f"  BOC k-eff = {record.boec_keff:.5f} ± {record.boec_keff_unc:.5f}")
            print(f"  EOC k-eff = {record.eoec_keff:.5f} ± {record.eoec_keff_unc:.5f}")
            print(f"  TRU in-core (BOC) = {record.tru_mass_kg:.1f} kg   "
                  f"MA/TRU = {record.ma_tru_frac:.4f}")
            print(f"  Fissile BOC -> EOC = {record.fissile_boc_kg:.1f} -> "
                  f"{record.fissile_eoc_kg:.1f} kg")
            if record.conversion_ratio is not None:
                print(f"  CR = {record.conversion_ratio:.3f}")
            if record.predicted_boec_keff is not None:
                print(f"  Reload scale {record.reload_scale:.4f} -> next BOEC k "
                      f"{record.predicted_boec_keff:.5f} (search)")
            print(f"  Reload: {record.fresh_tru_kg:.1f} kg TRU "
                  f"({record.external_tru_kg:.1f} kg external), "
                  f"TRU stock left {record.tru_stock_kg:.1f} kg")

            if stop_at_equilibrium and self.batch_mgr.is_at_equilibrium(
                state,
                window=self.cfg.equilibrium_window,
                keff_tolerance_pcm=self.cfg.keff_tolerance_pcm,
                cr_tolerance=self.cfg.cr_tolerance,
            ):
                print(f"\n  Equilibrium reached after cycle {cycle}.")
                break

        return state

    # --------------------------------------------------------------- private

    def _load_loadout(self) -> None:
        self._positions, self._asm_types = load_core(self.loadout_path)
        if self.cfg.cycle_length_days is None:
            with open(self.loadout_path) as f:
                self.cfg.cycle_length_days = float(json.load(f).get("cycle_length_days", 365.0))
        _, self._base_mats_dict = make_materials_from_loadout(self._asm_types)
        ag = AssemblyGeom()
        self._axial = load_axial(self.loadout_path)
        f2f = ag.duct_inner_flat_to_flat_cm
        self._assembly_vol_cm3 = (math.sqrt(3) / 2.0) * f2f ** 2 * self._axial.active_height_cm

        drivers, blankets = set(), set()
        for type_name, at in self._asm_types.items():
            universe = at.get("universe", "")
            mat_key = DRIVER_UNIVERSES.get(universe) or BLANKET_UNIVERSES.get(universe)
            if mat_key is None:
                continue
            (drivers if universe in DRIVER_UNIVERSES else blankets).add(type_name)
            ref = self._base_mats_dict[mat_key]
            self._ref_mats[type_name] = _clone_material(
                ref, f"_ref_{type_name}", self._assembly_vol_cm3
            )
        self._driver_types = frozenset(drivers)
        self._blanket_types = frozenset(blankets)

    @property
    def _fuel_types(self) -> frozenset:
        return self._driver_types | self._blanket_types

    @property
    def _state_path(self) -> str:
        return str(self.run_dir / self.STATE_FILE)

    def _cycle_dir(self, cycle: int) -> str:
        d = self.run_dir / f"cycle_{cycle:03d}"
        d.mkdir(exist_ok=True)
        return str(d)

    def _load_or_initialize_state(self) -> CycleState:
        if os.path.exists(self._state_path):
            state = CycleState.from_json(self._state_path)
            print(f"  Resuming from cycle {state.current_cycle} "
                  f"(loaded {self._state_path})")
        else:
            state = CycleState.initialize(
                self._positions,
                loadout_name=Path(self.loadout_path).stem,
                fuel_types=self._fuel_types,
            )
            self.batch_mgr.stagger_initial_core(state)
            print(f"  Initialising fresh cycle state for '{state.loadout_name}' "
                  f"({len(state.in_core)} assemblies)")
        return state

    def _export(self, cycle_dir: str, results, index: int) -> dict:
        """{material name: openmc.Material} at depletion step `index`."""
        # Re-reading materials.xml reuses ids already taken in this process;
        # the copies are only read from, so the duplicate-id warnings are noise.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", openmc.IDWarning)
            mats = results.export_to_materials(
                index, path=os.path.join(cycle_dir, "materials.xml")
            )
        out = {}
        for m in mats:
            if m.depletable:
                m.volume = self._assembly_vol_cm3
                out[m.name] = m
        return out

    def _eoc(self, cycle: int) -> dict:
        """EOC materials of a completed cycle, cached."""
        if cycle not in self._eoc_cache:
            cdir = self._cycle_dir(cycle)
            dep_file = os.path.join(cdir, "depletion_results.h5")
            if not os.path.exists(dep_file):
                raise FileNotFoundError(
                    f"Missing depletion results for cycle {cycle}: {dep_file}"
                )
            self._eoc_cache = {cycle: self._export(cdir, openmc.deplete.Results(dep_file), -1)}
        return self._eoc_cache[cycle]

    def _build_per_assembly_mats(self, in_core: dict, cycle: int) -> dict:
        """
        {material_id: openmc.Material} for every assembly in core this cycle:
          - cycle-0 assemblies: copy of the zone-default material
          - assemblies fabricated for this cycle: built from fresh_hm_vec
          - everything else: EOC composition from the previous cycle,
            fission products included
        """
        per_asm: dict = {}
        for mid, entry in in_core.items():
            ref = self._ref_mats[entry.assembly_type]
            if entry.cycle_loaded == cycle and entry.fresh_hm_vec:
                per_asm[mid] = _fresh_material(
                    mid, ref, entry.fresh_hm_vec, self._assembly_vol_cm3
                )
            elif cycle == 0:
                per_asm[mid] = _clone_material(ref, mid, self._assembly_vol_cm3)
            else:
                prev_eoc = self._eoc(cycle - 1)
                per_asm[mid] = _clone_material(prev_eoc[mid], mid, self._assembly_vol_cm3)
        return per_asm

    def _positions_for(self, in_core: dict) -> list[AssemblyPosition]:
        """Lattice positions with each fuel slot renamed to the assembly occupying it."""
        slot_to_asm = {e.position_id: mid for mid, e in in_core.items()}
        out = []
        for p in self._positions:
            if p.type in self._fuel_types:
                if p.material_id not in slot_to_asm:
                    raise RuntimeError(f"Fuel slot {p.material_id} is empty")
                p = replace(p, material_id=slot_to_asm[p.material_id])
            out.append(p)
        return out

    def _deplete(self, model: openmc.Model, cycle_dir: str) -> str:
        """Run (or resume) this cycle's depletion; return the results path."""
        timesteps = self.cfg.timesteps
        dep_h5 = os.path.join(cycle_dir, "depletion_results.h5")
        prev_results = None
        remaining = timesteps
        if os.path.exists(dep_h5):
            prev_results = openmc.deplete.Results(dep_h5)
            done = len(prev_results) - 1   # len includes the initial state
            if done >= len(timesteps):
                print(f"  Depletion already complete ({done} steps); reusing results")
                return dep_h5
            remaining = timesteps[done:]
            print(f"  Resuming depletion at step {done}/{len(timesteps)}")

        orig_dir = os.getcwd()
        os.chdir(cycle_dir)
        try:
            operator = openmc.deplete.CoupledOperator(
                model, self.cfg.chain_file, prev_results=prev_results
            )
            integrator = openmc.deplete.CECMIntegrator(
                operator,
                timesteps=remaining,
                power=self.cfg.power_W,
                timestep_units="d",
            )
            integrator.integrate()
        finally:
            os.chdir(orig_dir)
        return dep_h5

    def _run_one_cycle(self, state: CycleState, cycle: int) -> CycleRecord:
        """Build model, run depletion, update state, return CycleRecord."""
        cycle_dir = self._cycle_dir(cycle)
        prof = RUN_PROFILES[self.cfg.profile]

        # ---- Build per-assembly materials and model -------------------------
        per_asm_mats = self._build_per_assembly_mats(state.in_core, cycle)
        model = build_model_per_assembly(
            positions=self._positions_for(state.in_core),
            per_asm_fuel_mats=per_asm_mats,
            base_mats_dict=self._base_mats_dict,
            settings_override={
                "particles": prof["particles"],
                "batches":   prof["batches"],
                "inactive":  prof["inactive"],
            },
            axial=self._axial,
        )
        model.export_to_xml(directory=cycle_dir)

        if cycle == 0:
            boc_hm = _sum_comps(_hm_grams(m) for m in per_asm_mats.values())
            state.add_to_ledger("initial_tru", _tru_g(boc_hm))
            state.add_to_ledger("initial_u", sum(boc_hm.values()) - _tru_g(boc_hm))

        # ---- Run depletion and extract results ------------------------------
        results = openmc.deplete.Results(self._deplete(model, cycle_dir))
        boec_k, boec_u, eoec_k, eoec_u = _read_keff(results)

        boc_mats = self._export(cycle_dir, results, 0)
        eoc_mats = self._export(cycle_dir, results, -1)
        self._eoc_cache = {cycle: eoc_mats}
        in_core_ids = list(state.in_core)
        boc_hm = _sum_comps(_hm_grams(boc_mats[mid]) for mid in in_core_ids)
        eoc_hm = _sum_comps(_hm_grams(eoc_mats[mid]) for mid in in_core_ids)
        tru_g, ma_frac = compute_tru_metrics(boc_hm)
        cr = _compute_cr(cycle_dir, self.cfg.timesteps)

        # ---- Discharge and cooling ------------------------------------------
        ready_for_nfrc, newly_discharged = self.batch_mgr.step(state, cycle)
        for entry in newly_discharged:
            entry.hm_g = _hm_grams(eoc_mats[entry.material_id])

        # ---- NFRC: decay, electrorefine, add to U/TRU stocks ----------------
        nfrc = self.cfg.nfrc
        discharged_tru_g = 0.0
        for entry in ready_for_nfrc:
            cooling = self.cfg.batch.for_type(entry.assembly_type).cooling_cycles
            hm = decay(entry.hm_g, cooling * self.cfg.cycle_length_days / 365.25)
            tru, u, waste = electrorefine(hm, nfrc.tru_recovery, nfrc.u_recovery)
            for stock, product in ((state.tru_stock_g, tru), (state.u_stock_g, u)):
                for n, g in product.items():
                    stock[n] = stock.get(n, 0.0) + g
            discharged_tru_g += _tru_g(hm)
            state.add_to_ledger("nfrc_in", sum(entry.hm_g.values()))
            state.add_to_ledger("decay_loss", sum(entry.hm_g.values()) - sum(hm.values()))
            state.add_to_ledger("waste", sum(waste.values()))

        # ---- Fabricate fresh assemblies for the vacated slots ---------------
        scale, predicted = self._search_reload(state, cycle)
        fresh_entries, fab_totals = self._reload(state, cycle, scale)
        self.batch_mgr.load_fresh(state, fresh_entries)

        return CycleRecord(
            cycle=cycle,
            boec_keff=boec_k,
            boec_keff_unc=boec_u,
            eoec_keff=eoec_k,
            eoec_keff_unc=eoec_u,
            conversion_ratio=cr,
            tru_mass_kg=tru_g / 1000.0,
            ma_tru_frac=ma_frac,
            discharged_tru_kg=discharged_tru_g / 1000.0,
            fresh_tru_kg=fab_totals["tru"] / 1000.0,
            fissile_boc_kg=_fissile_g(boc_hm) / 1000.0,
            fissile_eoc_kg=_fissile_g(eoc_hm) / 1000.0,
            external_tru_kg=fab_totals["external"] / 1000.0,
            makeup_u_kg=fab_totals["makeup_u"] / 1000.0,
            tru_stock_kg=sum(state.tru_stock_g.values()) / 1000.0,
            reload_scale=scale,
            reload_tru_frac=fab_totals["tru_frac"],
            predicted_boec_keff=predicted,
        )

    def _search_reload(self, state: CycleState, cycle: int) -> tuple[float, Optional[float]]:
        """
        Choose the reload scale (fresh driver TRU fraction / startup TRU fraction)
        so the next cycle's BOEC k-eff hits k_target.  Each trial fabricates on
        copies of the stocks and runs a transport-only solve of the next core.

        Returns (scale, predicted BOEC k); (1.0, None) when the search is off or
        no driver slot is being reloaded.  Halts the run if the best loading
        still exceeds k_max (T7.4) or the search does not converge.
        """
        vacant_drivers = [
            p for p in self.batch_mgr.vacant_fuel_positions(
                state, self._positions, fuel_types=self._fuel_types)
            if p.type in self._driver_types
        ]
        if not self.cfg.search_reload or not vacant_drivers:
            return 1.0, None

        cfg = self.cfg
        search_root = Path(self._cycle_dir(cycle)) / "reload_search"
        n_iter = 0

        def evaluate(scale: float) -> tuple[float, float]:
            nonlocal n_iter
            fresh, _ = self._reload(state, cycle, scale, dry_run=True)
            in_core = dict(state.in_core)
            in_core.update({e.material_id: e for e in fresh})
            workdir = search_root / f"iter_{n_iter:02d}"
            n_iter += 1
            k, k_unc = self._transport_keff(in_core, cycle + 1, workdir)
            print(f"  reload search: scale {scale:.5f} -> k = {k:.5f} ± {k_unc:.5f}")
            return k, k_unc

        res = secant_search(
            evaluate, cfg.k_target, x0=self._last_scale(state),
            tol_pcm=cfg.search_tol_pcm, max_iter=cfg.search_max_iter,
            bounds=(0.2, self._max_scale()),
        )
        with open(Path(self._cycle_dir(cycle)) / "reload_search.json", "w") as f:
            json.dump({"k_target": cfg.k_target, "profile": cfg.search_profile,
                       **res.to_dict()}, f, indent=2)
        shutil.rmtree(search_root, ignore_errors=True)

        if not res.converged:
            raise RuntimeError(
                f"Reload search after cycle {cycle} did not converge "
                f"(best scale {res.x:.4f}, k = {res.k:.5f} ± {res.k_unc:.5f}); "
                f"see reload_search.json"
            )
        if res.k - 2.0 * res.k_unc > cfg.k_max:
            raise RuntimeError(
                f"T7.4 halt: reload for cycle {cycle + 1} gives BOEC k = {res.k:.5f} "
                f"> {cfg.k_max}"
            )
        return res.x, res.k

    def _last_scale(self, state: CycleState) -> float:
        """Start each search from the previous reload's answer."""
        for r in reversed(state.cycle_records):
            if r.predicted_boec_keff is not None:
                return r.reload_scale
        return 1.0

    def _max_scale(self) -> float:
        """Keep the reload TRU fraction under the material clamp (0.35 TRU/HM)."""
        fracs = [self._startup_tru_frac(t) for t in self._driver_types]
        return min([2.0] + [0.35 / f for f in fracs if f > 0.0])

    def _startup_tru_frac(self, asm_type: str) -> float:
        return sum(f for n, f in _hm_vec(self._ref_mats[asm_type]).items()
                   if nuclide_element(n) in TRU_ELEMENTS)

    def _transport_keff(self, in_core: dict, cycle: int, workdir: Path) -> tuple[float, float]:
        """BOC k-eff of the core `in_core` would form at the start of `cycle`."""
        prof = RUN_PROFILES[self.cfg.search_profile]
        workdir.mkdir(parents=True, exist_ok=True)
        model = build_model_per_assembly(
            positions=self._positions_for(in_core),
            per_asm_fuel_mats=self._build_per_assembly_mats(in_core, cycle),
            base_mats_dict=self._base_mats_dict,
            settings_override={
                "particles": prof["particles"],
                "batches":   prof["batches"],
                "inactive":  prof["inactive"],
            },
            axial=self._axial,
        )
        with openmc.StatePoint(model.run(cwd=str(workdir))) as sp:
            k = sp.keff
        return k.n, k.s

    def _reload(
        self, state: CycleState, cycle: int, scale: float = 1.0, dry_run: bool = False,
    ) -> tuple[list[InCoreEntry], dict]:
        """
        Fabricate one fresh assembly per vacant slot.  Drivers get their zone's
        startup TRU/HM atom fraction × scale: recovered TRU first, the external
        feed (reactor-grade Pu by default) for any shortfall.  Blankets are fresh
        makeup U at their startup isotopics.

        dry_run: fabricate from copies of the stocks and leave the ledger alone
        (used by the reload search to try a loading without committing it).
        """
        nfrc = self.cfg.nfrc
        makeup_u = {"U235": nfrc.du_u235_frac, "U238": 1.0 - nfrc.du_u235_frac}
        external = PU_VECTORS[nfrc.external_feed]
        tru_stock = dict(state.tru_stock_g) if dry_run else state.tru_stock_g
        u_stock = dict(state.u_stock_g) if dry_run else state.u_stock_g

        vacant = self.batch_mgr.vacant_fuel_positions(
            state, self._positions, fuel_types=self._fuel_types
        )
        by_type: dict[str, list[AssemblyPosition]] = {}
        for p in vacant:
            by_type.setdefault(p.type, []).append(p)

        totals = {"tru": 0.0, "external": 0.0, "makeup_u": 0.0, "tru_frac": {}}
        fresh: list[InCoreEntry] = []
        for asm_type, slots in sorted(by_type.items()):
            ref = self._ref_mats[asm_type]
            if asm_type in self._driver_types:
                tru_frac = self._startup_tru_frac(asm_type) * scale
                if tru_frac <= 0.0:
                    raise NotImplementedError(
                        f"'{asm_type}' starts on U-only fuel; its U/TRU transition "
                        f"reload (recovered TRU plus external enriched U) is not "
                        f"implemented yet."
                    )
                fab = fabricate(
                    len(slots), _hm_mol(ref), tru_frac, tru_stock, u_stock,
                    makeup_u_vec=makeup_u, external_vec=external,
                    recycle_u=nfrc.recycle_u,
                )
                totals["tru_frac"][asm_type] = tru_frac
            else:
                fab = fabricate(
                    len(slots), _hm_mol(ref), 0.0, {}, {},
                    makeup_u_vec=_hm_vec(ref), external_vec={},
                )

            ext_g = sum(fab.from_external.values())
            makeup_g = sum(fab.from_makeup_u.values())
            totals["tru"] += sum(fab.from_tru_stock.values()) + ext_g
            totals["external"] += ext_g
            totals["makeup_u"] += makeup_g
            if not dry_run:
                state.add_to_ledger("external_tru", ext_g)
                state.add_to_ledger("makeup_u", makeup_g)
                state.add_to_ledger("fabricated", fab.total_g)

            for p in slots:
                mid = assembly_id(p.material_id, cycle + 1)
                fresh.append(InCoreEntry(
                    material_id=mid,
                    q=p.q, r=p.r,
                    assembly_type=p.type,
                    cycle_loaded=cycle + 1,
                    cycles_in_core=0,
                    position_id=p.material_id,
                    fresh_hm_vec=fab.hm_vec,
                ))

        # Stocks are mutated in place by fabricate(); drop drained isotopes
        for stock in (tru_stock, u_stock):
            for n in [n for n, g in stock.items() if g <= 1e-9]:
                del stock[n]
        return fresh, totals


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(
        description="Multi-cycle PRISM recycle run (reactor-grade Pu by default)."
    )
    ap.add_argument("loadout", nargs="?", default="core_loadouts/prism_breakeven.json",
                    help="core loadout JSON (default: RG-Pu breakeven core)")
    ap.add_argument("--cycles", type=int, default=20, help="maximum cycles to run")
    ap.add_argument("--run-dir", help="default: runs/<loadout>__cycles[__<profile>]")
    ap.add_argument("--profile", choices=RUN_PROFILES, default="standard",
                    help="transport statistics (default: standard)")
    ap.add_argument("--cycle-days", type=float,
                    help="irradiation days per cycle (default: the loadout's Triplett value)")
    ap.add_argument("--external-feed", choices=PU_VECTORS, default="reactor_grade",
                    help="Pu vector for startup / shortfall TRU (default: reactor_grade)")
    ap.add_argument("--no-recycle-u", action="store_true",
                    help="fresh DU matrix only; discard recovered U")
    ap.add_argument("--blanket-cycles", type=int, default=4,
                    help="blanket residence in cycles (default: 4)")
    ap.add_argument("--no-stop", action="store_true",
                    help="run all cycles even after equilibrium is detected")
    ap.add_argument("--no-search", action="store_true",
                    help="reload at the startup TRU fraction instead of searching")
    ap.add_argument("--k-target", type=float, default=1.005,
                    help="BOEC k-eff target for the reload search (default: 1.005)")
    ap.add_argument("--search-profile", choices=RUN_PROFILES, default="quick",
                    help="transport statistics for reload-search solves (default: quick)")
    ap.add_argument("--k-max", type=float, default=1.010,
                    help="BOEC k-eff above which a reload halts the run (default: 1.010, T7.4)")
    ap.add_argument("--max-step-days", type=float, default=60.0,
                    help="longest depletion step after the short leading steps (default: 60)")
    opts = ap.parse_args()

    name = Path(opts.loadout).stem
    run_dir = opts.run_dir or (
        f"runs/{name}__cycles" if opts.profile == "standard"
        else f"runs/{name}__cycles__{opts.profile}"
    )
    cfg = CycleRunnerConfig(
        cycle_length_days=opts.cycle_days,
        profile=opts.profile,
        batch=BatchConfig.prism_breeder(blanket_cycles=opts.blanket_cycles),
        nfrc=NFRCConfig(external_feed=opts.external_feed, recycle_u=not opts.no_recycle_u),
        search_reload=not opts.no_search,
        k_target=opts.k_target,
        k_max=opts.k_max,
        max_step_days=opts.max_step_days,
        search_profile=opts.search_profile,
    )
    runner = CycleRunner(loadout_path=opts.loadout, run_dir=run_dir, config=cfg)
    final_state = runner.run(n_cycles=opts.cycles, stop_at_equilibrium=not opts.no_stop)
    print(f"\nCompleted {len(final_state.cycle_records)} cycles.")
    print(f"State saved to: {runner._state_path}")
