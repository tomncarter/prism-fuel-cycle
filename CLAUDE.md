# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

All simulation commands must be run through the Docker wrapper:

```bash
# Transport-only (quick test, ~minutes)
./run_openmc.sh python build_prism_core.py core_loadouts/prism_unf_recycle.json

# Full depletion (~90 min)
./run_openmc.sh python build_prism_core.py core_loadouts/prism_unf_recycle.json --deplete

# Extract results to Parquet
./run_openmc.sh python postprocess.py
```

Run profiles (`--profile smoke|quick|cycle|standard|production`, default `standard`) set
particles/batches and the depletion schedule — see `RUN_PROFILES` in
`build_prism_core.py`. Non-standard profiles write to `runs/<loadout>__<profile>/`.
```bash
# Fuel-test matrix (core_loadouts/fueltest_*.json: HALEU, WG-U, RG-Pu, WG-Pu, WG-U+Pu)
./run_fuel_tests.sh quick             # BOC k-eff screen, then postprocess
./run_fuel_tests.sh quick --deplete   # ~1-year depletion (6 steps)
```

Geometry figures (2-D slices + 3-D ray-traced views, plot mode only, no transport):
```bash
./run_openmc.sh python visualize_geometry.py [loadout]   # default prism_breeder → figures/geometry/<loadout>/
```

Fuel-volume-fraction sweep (loading search per Table 14-3 pin design, then TRU and CR):
```bash
./run_openmc.sh python -m fuel_cycle.vf_sweep run core_loadouts/prism_breeder.json   # → runs/vf_sweep/
./run_openmc.sh python -m fuel_cycle.vf_sweep summarize                              # summary.csv + vf_sweep.png
```

Per-assembly burnup for a cycle run (written automatically when `fuel_cycle.cycle_runner` finishes):
```bash
./run_openmc.sh python -m fuel_cycle.burnup_map runs/<run>__cycles__cycle [--cycle N | --all]   # → <run>/burnup/
```

Files under `runs/` are written by the container as root — delete them with
`./run_openmc.sh rm -rf runs/<dir>`.

Rebuild the Docker image after Dockerfile changes:
```bash
docker build -t prism-openmc .
```

`./run_openmc.sh` pins each container to one NUMA node (24 cores / 48 threads):
`NUMA_NODE=0|1` picks the node (default 0), `THREADS=n` overrides the thread count
(default: all of the node's threads). Run two jobs at once on nodes 0 and 1.

## Architecture

### Model tiers
- **Tier 1 (current)**: Homogenized zones, shared material per assembly type. Fast but no spatial burnup detail.
- **Tier 2 (working; benchmarking against Triplett Table II)**: Unique material per assembly, GEMs/control rods, criticality search, equilibrium cycle via `fuel_cycle/`.
- **Tier 3 (future)**: Pin-level geometry.

### Core data flow

```
JSON loadout → core_maps/loader.py → build_prism_core.py → OpenMC XML → statepoint.*.h5
                                                                       ↓
                                                              postprocess.py → results/*.parquet
```

1. **`build_prism_core.py`** — Main entry point. Loads a JSON loadout, calls `geometry/` and `materials/` modules to build the OpenMC `Model`, sets up tallies, exports XML, and runs transport or depletion.
2. **`geometry/`** — Three modules: `params.py` (frozen dataclasses for dimensions), `universe_factory.py` (builds the 3-cell hex assembly universes), `hexmap.py` (hex coordinate math + ring-ordered lattice builder).
3. **`materials/tier1_homog.py`** — All Tier-1 material definitions. Tunable constants (`FISSILE_SCALE`, `PU_HM_INNER/OUTER`, smear fractions) live at the top. `make_materials_from_loadout()` is the primary factory.
4. **`core_maps/`** — `loader.py` parses JSON loadouts into `AssemblyPosition` objects. `abr_like_assignment.py` provides the default ring-based zoning (397 assemblies, rings 0–11).
5. **`postprocess.py`** — Reads all `runs/*/` directories, extracts depletion results into 5 Parquet files under `results/`.
6. **`fuel_cycle/`** — Tier-2 multi-cycle loop, run via `python -m fuel_cycle.cycle_runner`; builds per-assembly models with `build_prism_core.build_model_per_assembly`: `state.py` (persistent cycle JSON), `nfrc.py` (electrorefining + fuel fab), `batch_manager.py` (discharge/cooling/reload scheduling), `cycle_runner.py` (per-cycle orchestration, reload k search), plus `loading_search.py`, `search.py`, `mass_balance.py`, `burnup_map.py`, `vf_sweep.py`. Drivers with no Pu (HALEU, WG-U) can't be reloaded yet (`NotImplementedError` in `CycleRunner._reload`).

### JSON loadout format

Core configurations live in `core_loadouts/*.json`. Layout modes:
- `"rings_with_overrides"` — ring ranges, then per-position overrides; `"sym": 3|6` expands an
  override to its rotational images. Used by the four Triplett cores.
- `"rings"` — ring ranges only (e.g., `{"range": [0, 3], "type": "inner_fuel"}`)
- `"positions"` — explicit per-slot `(q, r)` coordinates (for shuffled Tier-2 states)

`"expected_counts"` ({type: n}) makes the loader fail if the core doesn't match (Triplett
Table III for the PRISM cores). `tests/test_loader.py` loads every loadout.

A loadout may set `active_height_cm` (default 101.6, the `Prism3D` value).

Each `assembly_type` entry specifies the OpenMC universe name and `fuel_params` (e.g., `pu_hm_frac`, `zr_frac`).

### Hex lattice convention

The OpenMC `HexLattice` (default orientation) has flat-top cells in vertical columns; each
ring starts directly above the centre and runs **clockwise**. `hexmap.axial_to_xy` must use
the same flat-top frame (`(0, n)` straight up), and `ring_angle_sort` sorts clockwise from
90°. Before 2026-11-15 `axial_to_xy` was pointy-top, which shifted each ring by about n/2
cells, so overrides (rods, GEMs, notches) landed in the wrong slots. `tests/test_hexmap.py`
checks every slot against the real lattice.

### Root universe structure

9 cells: active core (hex lattice), IHX annulus, vessel wall, upper/lower plenums, upper/lower end plates, outer void. Defined inline in `build_prism_core.py`. `VesselGeom` dataclass in `geometry/params.py` holds all radial/axial bounds.

### Tallies

- `fission_inner` / `fission_outer` — fission rate per driver zone (CellFilter on that zone's interior cells)
- `flux_spectrum_*` — energy-binned flux per zone
- `blanket_breeding_internal` / `blanket_breeding_radial` — absorption + fission scores (separate tallies, not combined — ordering was non-deterministic with a shared CellFilter)
- `leakage_radial_total`, `leakage_axial_top/bot` + spectra — on vacuum boundary surfaces

### Key constants (all in `materials/tier1_homog.py`; assembly geometry in `geometry/params.py` is ANL-AFCI-202's 16.142 cm pitch / 0.394 cm duct)

| Constant | Value | Notes |
|---|---|---|
| `FISSILE_SCALE` | 1.0 | Criticality tuning knob — scales Pu/fissile content |
| `PU_HM_INNER` | 0.20 | Inner driver Pu/HM atom fraction |
| `PU_HM_OUTER` | 0.17 | Outer driver Pu/HM atom fraction |
| `PIN_DESIGNS` / `DEFAULT_PIN_DESIGN` | `anl_sprism_driver` | Whole-cell volume fractions per pin design, converted to interior fractions around the explicit duct/gap (`recipe_from_cell_fractions`; bond Na → sodium). Sources: ANL-AFCI-202 (`sources/134264.pdf`: `anl_abr_driver`, `anl_sprism_driver`, `anl_hcr_driver` 271 pins; `anl_blanket` 127 pins) and Till & Chang Table 14-3 (`tc_d090`–`tc_d120`). Loadouts set top-level `"pin_design"`, per-type `fuel_params.pin_design`, or `"volume_fractions"`; `"legacy"` = pre-2026-11-15 interior 0.27/0.53/0.20 |
| `ZR_METAL_FRAC` | 0.10 | Zr **weight** fraction of fuel alloy |
| `T_INTERIOR_K` | 900 | Smeared interior temperature (assumed, no paper value) |
| `T_COOLANT_K` | 702.65 | Mean of 360/499 °C sodium (Triplett Table I); sets Na density |
| Thermal power | 840 MWth | Used for depletion normalization |
| Depletion step | 30 days | Default timestep |

### Temperature

`run_openmc.sh` uses the official ENDF/B-VIII.0 library (250–2500 K) with temperature
interpolation.  Runs before this change used 293.6 K data and 0.92 g/cm³ sodium, so their
k-eff, reload loadings and matched-k enrichments are cold-state values.

### OpenMC API notes (v0.15.3)

- `HexagonalPrism`: use `-prism` for interior region (not `~prism`)
- `Model.export_to_xml(directory=...)` not `path=`
- `IndependentSource(constraints={"fissionable": True})` not `only_fissionable=`
- `settings.statepoint["filename"]` is ignored by OpenMC; actual output is `statepoint.{batch}.h5`; use the path returned by `model.run()`

### Results schema

All Parquet files in `results/` join on `run_tag`. Key files:
- `metadata.parquet` — provenance (run_tag, statepoint path, OpenMC version, particle counts)
- `keff.parquet`, `inventories.parquet`, `tallies.parquet`, `spectra.parquet`

Analysis notebooks in `notebooks/` (`01_criticality.ipynb` through `05_inventories.ipynb`) demonstrate standard queries against these files.
