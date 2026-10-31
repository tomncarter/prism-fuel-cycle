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
| `FISSILE_SCALE` | 1.0 | Global Pu multiplier for criticality tuning (Tier-1 only) |
| `ZR_METAL_FRAC` | 0.10 | Zr **weight** fraction of the U-(Pu)-Zr alloy (default; overridable via `fuel_params.zr_frac`) |
| `SMEAR_FUEL_METAL` | 0.27 | Interior volume fraction of fuel slug (~0.36 clad-inner × 0.75 smear) |
| `SMEAR_SODIUM` | 0.53 | Interior volume fraction of coolant + bond Na |
| `SMEAR_STEEL` | 0.20 | Interior volume fraction of HT9 cladding + wire wrap |
| `NU_U235_FRAC` | 0.00711 | U-235 atom fraction for natural uranium blanket feed |
| `DU_U235_FRAC` | 0.0020 | U-235 atom fraction for depleted uranium blanket feed |

### Material functions

**`make_smeared_driver(name, pu_hm, pu_vector, u_vector, zr_frac)`** — driver zone interior.
Computes the U-Pu-Zr (or U-Zr when `pu_hm == 0`) metal composition, applies
`FISSILE_SCALE` to the Pu fraction only, then blends with Na and Fe/Cr at the smear
fractions.  `make_smeared_fuel` takes each component at its theoretical density
(alloy by specific-volume mixing of 19.05 g/cm³ HM and 6.52 g/cm³ Zr; Na 0.92;
HT9 7.8) weighted by its volume fraction, so amounts are true atom/b-cm
(`set_density("sum")`).

- `pu_vector`: `"reactor_grade"` (Pu239=55%, Pu240=25%, Pu241=12%, Pu242=8%) or
  `"weapons_grade"` (93.4/6.0/0.5/0.1%).
- `u_vector`: a `U_VECTORS` name — `"u238"` (legacy default, pure U-238),
  `"depleted"`, `"natural"`, `"haleu"` (19.75% U-235), `"weapons_grade"` (93.2% HEU)
  — or a blend dict by atom weight, e.g. `{"weapons_grade": 0.224, "depleted": 0.776}`
  for HEU down-blended into DU.

Loadout `fuel_params` map directly: `pu_hm_frac`, `pu_vector`, `u_vector`.

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
