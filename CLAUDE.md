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

Files under `runs/` are written by the container as root — delete them with
`./run_openmc.sh rm -rf runs/<dir>`.

Rebuild the Docker image after Dockerfile changes:
```bash
docker build -t prism-openmc .
```

Exit code 127 from `./run_openmc.sh` is harmless (shell wrapper artifact).

## Architecture

### Model tiers
- **Tier 1 (current)**: Homogenized zones, shared material per assembly type. Fast but no spatial burnup detail.
- **Tier 2 (planned)**: Unique material per assembly, GEMs/control rods, criticality search, equilibrium cycle via `fuel_cycle/`.
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
6. **`fuel_cycle/`** — Tier-2 infrastructure (not yet wired to main build): `state.py` (persistent cycle JSON), `nfrc.py` (electrorefining + fuel fab), `batch_manager.py` (discharge/cooling/reload scheduling), `cycle_runner.py` (per-cycle orchestration).

### JSON loadout format

Core configurations live in `core_loadouts/*.json`. Two layout modes:
- `"rings"` — compact range notation (e.g., `{"range": [0, 3], "type": "inner_fuel"}`)
- `"positions"` — explicit per-slot `(q, r)` coordinates (for shuffled Tier-2 states)

A loadout may set `active_height_cm` (default 101.6, the `Prism3D` value).

Each `assembly_type` entry specifies the OpenMC universe name and `fuel_params` (e.g., `pu_hm_frac`, `zr_frac`).

### Hex lattice convention

OpenMC `HexagonalLattice` uses **clockwise** ring ordering starting from the +x direction. `hexmap.ring_angle_sort` sorts with `(start - a) % 360` ascending to match this convention. Getting this wrong silently mis-assigns assembly types.

### Root universe structure

9 cells: active core (hex lattice), IHX annulus, vessel wall, upper/lower plenums, upper/lower end plates, outer void. Defined inline in `build_prism_core.py`. `VesselGeom` dataclass in `geometry/params.py` holds all radial/axial bounds.

### Tallies

- `fission_split` — fission rate per zone (CellFilter on all fuel cells)
- `flux_spectrum_*` — energy-binned flux per zone
- `blanket_breeding_internal` / `blanket_breeding_radial` — absorption + fission scores (separate tallies, not combined — ordering was non-deterministic with a shared CellFilter)
- `leakage_radial_total`, `leakage_axial_top/bot` + spectra — on vacuum boundary surfaces

### Key constants (all in `materials/tier1_homog.py`)

| Constant | Value | Notes |
|---|---|---|
| `FISSILE_SCALE` | 1.0 | Criticality tuning knob — scales Pu/fissile content |
| `PU_HM_INNER` | 0.20 | Inner driver Pu/HM atom fraction |
| `PU_HM_OUTER` | 0.17 | Outer driver Pu/HM atom fraction |
| `SMEAR_FUEL_METAL` | 0.27 | Interior volume fraction (fuel slug, 75% smear) |
| `SMEAR_SODIUM` | 0.53 | Interior volume fraction (coolant + bond Na) |
| `SMEAR_STEEL` | 0.20 | Interior volume fraction (cladding + wire) |
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
