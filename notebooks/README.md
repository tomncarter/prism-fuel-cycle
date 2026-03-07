# notebooks/

Analysis notebooks. All read from `../results/*.parquet` — run `postprocess.py` first.

| Notebook | Contents |
|----------|----------|
| `01_criticality.ipynb` | k-eff vs time for each depleted run; reactivity swing |
| `02_breeding.ipynb` | Blanket capture rates, breeding ratio proxy, NU vs DU comparison |
| `03_flux_spectra.ipynb` | Energy-binned flux in driver and blanket zones |
| `04_leakage.ipynb` | Radial and axial leakage fractions, leakage spectra |
| `05_inventories.ipynb` | Nuclide mass evolution (U-235/238, Pu vector, Am, Np) over depletion |

## Running

Notebooks run inside the standard Jupyter environment (not inside the OpenMC
container — they only need `pandas`, `pyarrow`, and `matplotlib`):

```bash
jupyter lab notebooks/
```

Or from the container if preferred:

```bash
./run_openmc.sh bash -c "pip install jupyter -q && jupyter nbconvert --to notebook \
  --execute notebooks/01_criticality.ipynb --output notebooks/01_criticality.ipynb"
```

## Data joins

All parquet files share the `run_tag` column (e.g. `"NUint_NUrad"`). Join on
`run_tag` to correlate k-eff with inventory or tally data. The `metadata.parquet`
table is the authoritative record of what each `run_tag` represents.

```python
import pandas as pd
keff = pd.read_parquet("results/keff.parquet")
meta = pd.read_parquet("results/metadata.parquet")
df   = keff.merge(meta[["run_tag","internal_feed","radial_feed"]], on="run_tag")
```
