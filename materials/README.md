# materials/

## tier1_homog.py

All material definitions for the Tier-1 homogenized model. The key design choice is
that each assembly interior is a single smeared material — fuel metal, sodium coolant,
and steel structure are blended into one OpenMC `Material` with a volume-fraction
weighted composition.

### Tunable constants (top of file)

| Constant | Value | Description |
|----------|-------|-------------|
| `PU_HM_INNER` | 0.20 | Pu/HM atom fraction, inner driver zone |
| `PU_HM_OUTER` | 0.17 | Pu/HM atom fraction, outer driver zone |
| `FISSILE_SCALE` | 0.923 | Global Pu multiplier for criticality tuning (Tier-1 only) |
| `ZR_METAL_FRAC` | 0.10 | Zr fraction within fuel metal (U-Pu-Zr alloy) |
| `SMEAR_FUEL_METAL` | 0.35 | Volume fraction of fuel metal in smeared interior |
| `SMEAR_SODIUM` | 0.55 | Volume fraction of sodium in smeared interior |
| `SMEAR_STEEL` | 0.10 | Volume fraction of steel in smeared interior |
| `NU_U235_FRAC` | 0.00711 | U-235 atom fraction for natural uranium blanket feed |
| `DU_U235_FRAC` | 0.0020 | U-235 atom fraction for depleted uranium blanket feed |

### Material functions

**`make_smeared_driver(name, pu_hm)`** — driver zone interior. Computes the
U-Pu-Zr metal composition from `pu_hm`, applies `FISSILE_SCALE`, then blends with
Na and Fe/Cr at the smear fractions. The Pu isotopic vector is reactor-grade:
Pu239=55%, Pu240=25%, Pu241=12%, Pu242=8%.

**`make_smeared_blanket(name, u235_frac)`** — blanket interior. Same smear structure
but uranium only (no Pu at BOC). The `u235_frac` selects NU vs DU feed.

**`make_blended_na_steel(name, steel_vol_frac)`** — vessel/pool regions. Pure
Na+Fe/Cr blend with volume-averaged density. Used for IHX pool, upper plenum, lower
plenum.

**`make_materials(internal_feed, radial_feed)`** — top-level factory called by
`build_prism_core.py`. `internal_feed` and `radial_feed` are `"NU"` or `"DU"`.
Returns `(openmc.Materials, dict)` where the dict is keyed by material role name
(e.g. `"fuel_driver_inner_homog"`) for lookup by `UniverseFactory`.

### `FISSILE_SCALE` note

This knob exists to bring the Tier-1 model to near-criticality without a formal
criticality search. It will be retired in Tier-2 in favour of an explicit Pu/HM
search. It should not be used in production runs where absolute reactivity matters.

### Depletable flag

Driver and blanket materials are marked `depletable = True` (OpenMC default).
Structural materials (sodium, HT9, reflector, shield, B4C, vessel regions) are
explicitly set to `depletable = False` to exclude them from the depletion solver.
