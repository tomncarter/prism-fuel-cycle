# PRISM Fuel Cycle — OpenMC Simulation

Tier-1 (homogenized) neutronics model of a PRISM/ABR-like sodium-cooled fast reactor,
built to study blanket breeding performance under different uranium feed strategies.

## Architecture

```
build_prism_core.py      ← main entry point
    │
    ├── core_maps/       ← radial zone assignment (which ring = which assembly type)
    ├── geometry/        ← assembly geometry, hex lattice construction
    ├── materials/       ← homogenized material definitions
    │
    └── runs/            ← one subdirectory per permutation (git-ignored)
         └── NUint_NUrad/
              ├── statepoint.80.h5
              └── depletion_results.h5

postprocess.py           ← extracts results → results/*.parquet (git-ignored)
notebooks/               ← analysis notebooks reading from results/
```

## Permutation system

Four blanket feed configurations are run automatically:

| Directory | Internal blanket | Radial blanket |
|-----------|-----------------|----------------|
| `NUint_NUrad` | Natural U | Natural U |
| `NUint_DUrad` | Natural U | Depleted U |
| `DUint_NUrad` | Depleted U | Natural U |
| `DUint_DUrad` | Depleted U | Depleted U |

## Running

```bash
# Build custom Docker image (first time, or after Dockerfile changes)
docker build -t prism-openmc .

# Transport only (quick, ~5 min per permutation)
./run_openmc.sh python build_prism_core.py

# Full depletion (~90 min per permutation)
./run_openmc.sh python build_prism_core.py --deplete

# Extract results to Parquet
./run_openmc.sh python postprocess.py
```

## Results

`postprocess.py` writes five Parquet files to `results/`:

| File | Contents |
|------|----------|
| `keff.parquet` | k-eff vs time for each depleted run |
| `inventories.parquet` | Key nuclide masses vs time per material |
| `tallies.parquet` | Scalar tallies at final step (fission split, breeding, leakage) |
| `spectra.parquet` | Energy-binned flux and leakage spectra at final step |
| `metadata.parquet` | Run provenance (OpenMC version, particle counts, timestamps) |

## Tier roadmap

| Tier | Description | Status |
|------|-------------|--------|
| 1 | Homogenized geometry, no reactivity control | **Complete (POC)** |
| 2 | Homogenized + GEMs + control rods, criticality search | Planned |
| 3 | Pin-level geometry | Future |

## Key design parameters

| Parameter | Value |
|-----------|-------|
| Assembly pitch | 16.142 cm |
| Active height | 101.60 cm |
| Max ring | 11 |
| Pu/HM (inner) | 0.20 |
| Pu/HM (outer) | 0.17 |
| Thermal power | 840 MWth |
| Depletion step | 30 days |
| Pu vector | Pu239=55%, Pu240=25%, Pu241=12%, Pu242=8% |
