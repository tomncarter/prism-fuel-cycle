# materials/tier1_homog.py
from dataclasses import dataclass

import openmc

from geometry.params import AssemblyGeom

# -----------------------------
# Tunable knobs (Tier-1)
# -----------------------------
# Blanket feed options (U-235 atom fractions)
NU_U235_FRAC = 0.00711   # Natural uranium
DU_U235_FRAC = 0.0020    # Depleted uranium (0.20 wt% tails assay, atom-fraction approx)

BLANKET_FEEDS = {"NU": NU_U235_FRAC, "DU": DU_U235_FRAC}

# Pu isotopic vectors (atom fractions within Pu)
PU_VEC_REACTOR_GRADE = {
    "Pu239": 0.55,
    "Pu240": 0.25,
    "Pu241": 0.12,
    "Pu242": 0.08,
}

PU_VEC_WEAPONS_GRADE = {
    "Pu239": 0.934,
    "Pu240": 0.060,
    "Pu241": 0.005,
    "Pu242": 0.001,
}

PU_VECTORS = {
    "reactor_grade": PU_VEC_REACTOR_GRADE,
    "weapons_grade": PU_VEC_WEAPONS_GRADE,
}

# Default for backward compatibility
PU_VEC = PU_VEC_REACTOR_GRADE

# U isotopic vectors (atom fractions within U).  Tier-1 approximations;
# U-234 follows the enrichment cascade, U-236 only appears in recycled HEU.
U_VECTORS = {
    "u238":          {"U238": 1.0},   # legacy driver matrix (pre-fuel-test loadouts)
    "depleted":      {"U235": DU_U235_FRAC, "U238": 1.0 - DU_U235_FRAC},
    "natural":       {"U234": 0.000055, "U235": NU_U235_FRAC, "U238": 1.0 - NU_U235_FRAC - 0.000055},
    "haleu":         {"U234": 0.0018, "U235": 0.1975, "U238": 0.8007},
    "weapons_grade": {"U234": 0.0100, "U235": 0.9320, "U236": 0.0020, "U238": 0.0560},
}

# Fuel metal composition: Zr weight fraction of the U-(Pu)-Zr alloy
ZR_METAL_FRAC = 0.10  # U-10Zr / U-Pu-10Zr (10 wt%, ~23 at%)

# -----------------------------
# Pin design -> smeared volume fractions
# -----------------------------
# Each design gives whole-lattice-cell volume fractions (pin bundle + duct +
# inter-assembly gap).  "bond" is the bond sodium, modelled as in-core sodium;
# designs that exclude it (Till & Chang) carry bond = 0 and leave that volume void.
#
# ANL-AFCI-202, "Advanced Burner Reactor 1000MWth Reference Concept" (2007),
# sources/134264.pdf.  All at 16.142 cm pitch, 0.432 cm gap, 0.394 cm HT9 duct,
# i.e. this model's assembly geometry.  Volume fractions at fabrication / cold.
#   anl_abr_driver    Table II.1-2   reference ABR metal burner: 271 pins, 0.755 cm OD, 75% smear
#   anl_sprism_driver Table III.6-10 S-PRISM metal driver: 271 pins, 0.744 cm OD, 75% smear
#   anl_hcr_driver    Table III.6-10 high-conversion metal driver: 271 pins, 0.808 cm OD, 75% smear
#   anl_blanket       Table III.6-11 metal blanket (U-Zr): 127 pins, 1.201 cm OD, 85% smear
#                     (its wire wrap prints as 1.300 cm, evidently 0.130)
#
# Till & Chang, "Plentiful Energy" (2011), Table 14-3: 169-pin metal-fuel
# assemblies, 0.05 cm clad, 0.3 cm duct, pitch growing with pin OD (16.2-19.2 cm);
# bond sodium excluded.  Only the fractions are used here.  The 1.0 cm row prints
# fuel 0.415, which breaks the row sum and the Fig. 14-5 trend; 0.313 is used.
PIN_DESIGNS = {
    "anl_abr_driver":    dict(pins=271, pin_od_cm=0.755, fuel=0.292, bond=0.098, coolant=0.353, structure=0.257),
    "anl_sprism_driver": dict(pins=271, pin_od_cm=0.744, fuel=0.283, bond=0.094, coolant=0.366, structure=0.257),
    "anl_hcr_driver":    dict(pins=271, pin_od_cm=0.808, fuel=0.343, bond=0.114, coolant=0.285, structure=0.258),
    "anl_blanket":       dict(pins=127, pin_od_cm=1.201, fuel=0.446, bond=0.078, coolant=0.266, structure=0.210),
    "tc_d090": dict(pins=169, pin_od_cm=0.9, fuel=0.279, bond=0.0, coolant=0.377, structure=0.203),
    "tc_d100": dict(pins=169, pin_od_cm=1.0, fuel=0.313, bond=0.0, coolant=0.345, structure=0.189),
    "tc_d110": dict(pins=169, pin_od_cm=1.1, fuel=0.347, bond=0.0, coolant=0.316, structure=0.178),
    "tc_d120": dict(pins=169, pin_od_cm=1.2, fuel=0.376, bond=0.0, coolant=0.290, structure=0.168),
}
# Triplett's breakeven core (21.1% TRU, CR 1.06, 106/147 MWd/kg, 23 months)
# matches ANL's S-PRISM metal column (21.4%, 1.05, 106/149, 21.7 months).
DEFAULT_PIN_DESIGN = "anl_sprism_driver"


@dataclass(frozen=True)
class SmearRecipe:
    """Volume fractions of the assembly interior (inside the duct wall)."""
    fuel: float     # as-fabricated U-(Pu)-Zr slug
    sodium: float   # coolant (no bond sodium for Table 14-3 designs)
    steel: float    # cladding + wire wrap

    def __post_init__(self):
        for k in ("fuel", "sodium", "steel"):
            if not 0.0 <= getattr(self, k) <= 1.0:
                raise ValueError(f"interior {k} fraction {getattr(self, k):.4f} outside [0, 1]")
        if self.fuel + self.sodium + self.steel > 1.0 + 1e-9:
            raise ValueError(f"interior fractions sum to {self.fuel + self.sodium + self.steel:.4f} > 1")


# Pre-2026-11-15 recipe (interior fractions, bond sodium included); reproduces old runs.
LEGACY_RECIPE = SmearRecipe(fuel=0.27, sodium=0.53, steel=0.20)


def recipe_from_cell_fractions(fuel: float, coolant: float, structure: float,
                               geom: AssemblyGeom | None = None, bond: float = 0.0) -> SmearRecipe:
    """
    Interior fractions from whole-cell fractions.  The duct wall (steel) and the
    inter-assembly gap (sodium) are explicit cells in UniverseFactory, so their
    share of the lattice cell comes off the structure and coolant fractions.
    Bond sodium is added to the interior sodium.
    """
    g = geom or AssemblyGeom()
    a_cell = g.assembly_pitch_cm ** 2           # hex areas scale as flat-to-flat^2
    f_in = g.duct_inner_flat_to_flat_cm ** 2 / a_cell
    f_duct = (g.duct_flat_to_flat_cm ** 2 - g.duct_inner_flat_to_flat_cm ** 2) / a_cell
    f_gap = 1.0 - g.duct_flat_to_flat_cm ** 2 / a_cell
    return SmearRecipe(fuel=fuel / f_in,
                       sodium=(coolant + bond - f_gap) / f_in,
                       steel=(structure - f_duct) / f_in)


def resolve_recipe(fuel_params: dict | None = None) -> SmearRecipe:
    """
    Recipe for one assembly type from its fuel_params:
      "pin_design": a PIN_DESIGNS key, or "legacy"
      "volume_fractions": {"basis": "cell"|"interior", "fuel", "coolant", "structure", "bond"?}
    Neither -> DEFAULT_PIN_DESIGN.
    """
    fp = fuel_params or {}
    vf = fp.get("volume_fractions")
    if vf is not None:
        basis = vf.get("basis", "cell")
        if basis == "cell":
            return recipe_from_cell_fractions(vf["fuel"], vf["coolant"], vf["structure"],
                                              bond=vf.get("bond", 0.0))
        if basis == "interior":
            return SmearRecipe(fuel=vf["fuel"], sodium=vf["coolant"], steel=vf["structure"])
        raise ValueError(f"volume_fractions basis must be 'cell' or 'interior', got {basis!r}")
    name = fp.get("pin_design", DEFAULT_PIN_DESIGN)
    if name == "legacy":
        return LEGACY_RECIPE
    if name not in PIN_DESIGNS:
        raise ValueError(f"unknown pin_design {name!r}; choose from {sorted(PIN_DESIGNS)} or 'legacy'")
    d = PIN_DESIGNS[name]
    return recipe_from_cell_fractions(d["fuel"], d["coolant"], d["structure"], bond=d["bond"])


# Operating temperatures (K).  Sodium inlet/outlet are Triplett 2012 Table I
# (360 / 499 C).  T_INTERIOR_K is an assumption: the paper gives no fuel
# temperature.  A smeared interior holds fuel, Na and steel in one material, so
# it takes the temperature that matters most, the Doppler broadening of the
# actinides; ~600 C is typical for sodium-bonded metal fuel.
T_INLET_K    = 633.15
T_OUTLET_K   = 772.15
T_COOLANT_K  = 0.5 * (T_INLET_K + T_OUTLET_K)
T_INTERIOR_K = 900.0


def sodium_density(temperature_K: float) -> float:
    """Liquid sodium density (g/cm3), Fink & Leibowitz (ANL/RE-95/2), 371-2503.7 K."""
    tau = 1.0 - temperature_K / 2503.7
    return (219.0 + 275.32 * tau + 511.58 * tau ** 0.5) / 1000.0


# Theoretical densities (g/cm3) for the volume-fraction mixing below
_RHO_HM_METAL = 19.05   # alpha-U; delta-stabilised Pu is close enough for Tier-1
_RHO_ZR_METAL = 6.52

_RHO_NA    = sodium_density(T_COOLANT_K)   # at the mean coolant temperature
_RHO_STEEL = 7.80   # HT9/steel density (g/cm3)

# Pu/HM targets for inner/outer zones (start here, then tune)
PU_HM_INNER = 0.20
PU_HM_OUTER = 0.17

# Single scalar multiplier to tune criticality quickly
# Increase -> higher k. Start at 1.0 then try 1.1, 1.2, etc.
FISSILE_SCALE = 1

# -----------------------------
# Vessel / pool region smearing
# -----------------------------
# Volume fraction of steel (Fe/Cr) in each region; remainder is sodium.
# Increase -> more neutron absorption/scattering in that region.
IHX_STEEL_VOL_FRAC          = 0.30   # IHX tube bundles + structure in radial pool
UPPER_PLENUM_STEEL_VOL_FRAC = 0.15   # lighter above core (fission gas plenums, light structure)
LOWER_PLENUM_STEEL_VOL_FRAC = 0.35   # heavier below core (pumps, lower internals)


def _normalize(d):
    s = sum(d.values())
    return {k: v / s for k, v in d.items()}


def resolve_u_vector(u_vector) -> dict:
    """
    Resolve a U isotopic spec into normalized atom fractions.

    u_vector: a U_VECTORS name (e.g. "haleu"), or a blend of named vectors
              by atom weight, e.g. {"weapons_grade": 0.224, "depleted": 0.776}
              for HEU down-blended into a DU matrix.
    """
    if isinstance(u_vector, str):
        return _normalize(U_VECTORS[u_vector])
    blend: dict = {}
    for vec_name, weight in u_vector.items():
        for iso, frac in U_VECTORS[vec_name].items():
            blend[iso] = blend.get(iso, 0.0) + weight * frac
    return _normalize(blend)


def make_smeared_driver(
    name: str,
    pu_hm: float = 0.0,
    pu_vector: str = "reactor_grade",
    u_vector="u238",
    zr_frac: float = ZR_METAL_FRAC,
    recipe: SmearRecipe | None = None,
) -> openmc.Material:
    """
    Smeared interior material approximating: (fuel metal + Na + steel).
    Fuel metal is U-Pu-Zr (or U-Zr when pu_hm == 0) with the specified
    Pu and U isotopic vectors.

    pu_vector: "reactor_grade" (default) or "weapons_grade"
    u_vector:  U_VECTORS name or blend dict (see resolve_u_vector)
    zr_frac:   Zr weight fraction of the fuel alloy

    FISSILE_SCALE applies to the Pu fraction only; U-Zr cores are tuned via
    their U enrichment / blend.
    """
    pu_hm = max(0.0, min(0.35, pu_hm * FISSILE_SCALE))  # clamp for sanity
    u_hm = 1.0 - pu_hm

    pu_vec = _normalize(PU_VECTORS.get(pu_vector, PU_VEC_REACTOR_GRADE))
    u_vec  = resolve_u_vector(u_vector)

    hm_vec = {iso: u_hm * frac for iso, frac in u_vec.items()}
    if pu_hm > 0.0:
        for iso, frac in pu_vec.items():
            hm_vec[iso] = hm_vec.get(iso, 0.0) + pu_hm * frac

    return make_smeared_fuel(name, hm_vec, zr_frac=zr_frac, recipe=recipe)


def make_smeared_blanket(name: str, u235_frac: float,
                         recipe: SmearRecipe | None = None) -> openmc.Material:
    """
    Smeared blanket interior: NU or DU metal + sodium + steel, no Pu.
    u235_frac: atom fraction of U-235 (remainder is U-238).
    """
    return make_smeared_fuel(
        name, {"U235": u235_frac, "U238": 1.0 - u235_frac}, zr_frac=0.0, recipe=recipe
    )


def _avogadro_b_cm() -> float:
    return openmc.data.AVOGADRO * 1e-24


def _fuel_alloy_density(zr_frac: float) -> float:
    """U-(Pu)-Zr alloy density (g/cm3) by specific-volume mixing of HM and Zr."""
    return 1.0 / ((1.0 - zr_frac) / _RHO_HM_METAL + zr_frac / _RHO_ZR_METAL)


def _steel_number_density() -> float:
    """HT9 surrogate (88/12 at% Fe/Cr) atom density in atom/b-cm."""
    m = 0.88 * openmc.data.atomic_weight("Fe") + 0.12 * openmc.data.atomic_weight("Cr")
    return _RHO_STEEL * _avogadro_b_cm() / m


def make_smeared_fuel(name: str, hm_vec: dict, zr_frac: float = ZR_METAL_FRAC,
                      recipe: SmearRecipe | None = None) -> openmc.Material:
    """
    Smeared interior (fuel metal + Na + steel) for an arbitrary heavy-metal
    isotopic vector.  Shared recipe behind drivers, blankets and the fresh
    assemblies fabricated by fuel_cycle from recycled U/TRU.

    Each component is taken at its theoretical density and weighted by its
    interior volume fraction from `recipe` (default: DEFAULT_PIN_DESIGN), so
    nuclide amounts are true atom/b-cm.  Any remainder is void.

    hm_vec:  {nuclide: atom fraction within HM} (normalised here)
    zr_frac: Zr weight fraction of the fuel alloy (0 for blankets)
    """
    hm_vec = _normalize(hm_vec)
    na_av = _avogadro_b_cm()
    rc = recipe or resolve_recipe()

    m = openmc.Material(name=name)

    # Fuel slug: HM and Zr partial densities of the alloy, scaled by volume fraction.
    rho_alloy = _fuel_alloy_density(zr_frac)
    m_hm = sum(frac * openmc.data.atomic_mass(iso) for iso, frac in hm_vec.items())
    n_hm = rc.fuel * rho_alloy * (1.0 - zr_frac) * na_av / m_hm
    for iso, frac in hm_vec.items():
        if frac > 0.0:
            m.add_nuclide(iso, n_hm * frac)

    if zr_frac > 0.0:
        n_zr = rc.fuel * rho_alloy * zr_frac * na_av / openmc.data.atomic_weight("Zr")
        m.add_element("Zr", n_zr)

    # Coolant sodium (plus bond sodium in the legacy recipe):
    m.add_element("Na", rc.sodium * _RHO_NA * na_av / openmc.data.atomic_weight("Na"))

    # Cladding + wire wrap (simple Fe/Cr surrogate):
    n_steel = rc.steel * _steel_number_density()
    m.add_element("Fe", n_steel * 0.88)
    m.add_element("Cr", n_steel * 0.12)

    m.set_density("sum")
    m.temperature = T_INTERIOR_K
    return m


# Mass density of a fresh U-10Zr driver interior (default pin design); the fixed
# density used when rebuilding materials from evolved {nuclide: grams} compositions.
_RC0 = resolve_recipe()
RHO_INTERIOR = (_RC0.fuel * _fuel_alloy_density(ZR_METAL_FRAC)
                + _RC0.sodium * _RHO_NA + _RC0.steel * _RHO_STEEL)


def make_blended_na_steel(
    name: str, steel_vol_frac: float, temperature_K: float = T_COOLANT_K
) -> openmc.Material:
    """
    Homogenized sodium + steel (Fe/Cr) region.
    steel_vol_frac: volume fraction of steel; remainder is sodium.
    temperature_K:  region temperature; also sets the sodium density.
    Composition amounts are volume-fraction weighted (OpenMC normalises internally).
    Density is volume-averaged.
    """
    na_frac = 1.0 - steel_vol_frac
    m = openmc.Material(name=name)
    m.add_element("Na", na_frac)
    m.add_element("Fe", steel_vol_frac * 0.88)
    m.add_element("Cr", steel_vol_frac * 0.12)
    m.set_density("g/cm3", na_frac * sodium_density(temperature_K) + steel_vol_frac * _RHO_STEEL)
    m.temperature = temperature_K
    return m


# Maps universe name -> (mats_dict key, material kind)
# "kind" tells make_materials_from_loadout which factory to call.
_UNIVERSE_TO_MAT = {
    "u_inner_fuel":       ("fuel_driver_inner_homog",    "driver"),
    "u_outer_fuel":       ("fuel_driver_outer_homog",    "driver"),
    "u_internal_blanket": ("fuel_blanket_internal_homog","blanket"),
    "u_radial_blanket":   ("fuel_blanket_radial_homog",  "blanket"),
}


def _make_structural_mats() -> dict:
    """Build the set of non-fuel materials shared by all configurations."""
    mats = {}

    sodium = openmc.Material(name="sodium")
    sodium.add_element("Na", 1.0)
    sodium.set_density("g/cm3", _RHO_NA)
    sodium.depletable = False
    mats["sodium"] = sodium

    ht9 = openmc.Material(name="ht9")
    ht9.add_element("Fe", 0.88)
    ht9.add_element("Cr", 0.12)
    ht9.set_density("g/cm3", 7.8)
    ht9.depletable = False
    mats["ht9"] = ht9

    reflector = openmc.Material(name="reflector")
    reflector.add_element("Fe", 1.0)
    reflector.set_density("g/cm3", 7.8)
    reflector.depletable = False
    mats["reflector"] = reflector

    shield = openmc.Material(name="shield")
    shield.add_element("Fe", 1.0)
    shield.set_density("g/cm3", 7.8)
    shield.depletable = False
    mats["shield"] = shield

    b4c = openmc.Material(name="b4c")
    b4c.add_element("B", 4.0)
    b4c.add_element("C", 1.0)
    b4c.set_density("g/cm3", 2.52)
    b4c.depletable = False
    mats["b4c"] = b4c

    b4c_ult = openmc.Material(name="b4c_ult")
    b4c_ult.add_element("B", 4.0)
    b4c_ult.add_element("C", 1.0)
    b4c_ult.set_density("g/cm3", 2.6)
    b4c_ult.depletable = False
    mats["b4c_ult"] = b4c_ult

    for m in mats.values():
        m.temperature = T_COOLANT_K

    for key, name, frac, temperature_K in [
        ("pool_ihx",     "pool_ihx",     IHX_STEEL_VOL_FRAC,          T_COOLANT_K),
        ("upper_plenum", "upper_plenum", UPPER_PLENUM_STEEL_VOL_FRAC, T_OUTLET_K),
        ("lower_plenum", "lower_plenum", LOWER_PLENUM_STEEL_VOL_FRAC, T_INLET_K),
    ]:
        m = make_blended_na_steel(name, frac, temperature_K)
        m.depletable = False
        mats[key] = m

    return mats


def make_materials_from_loadout(asm_types: dict):
    """
    Build materials from a loadout's assembly_types dict.

    Replaces make_materials() for loadout-driven runs. Produces the same
    mats_dict key names expected by UniverseFactory.build_all().

    Returns (openmc.Materials, mats_dict) — same signature as make_materials().
    """
    mats = _make_structural_mats()

    for type_name, at in asm_types.items():
        universe = at.get("universe", "")
        mapping = _UNIVERSE_TO_MAT.get(universe)
        if mapping is None:
            continue  # reflector, shield, control etc. — structural only

        mat_key, kind = mapping
        fp = at.get("fuel_params", {})

        if kind == "driver":
            mats[mat_key] = make_smeared_driver(
                name=mat_key,
                pu_hm=fp.get("pu_hm_frac", 0.0),
                pu_vector=fp.get("pu_vector", "reactor_grade"),
                u_vector=fp.get("u_vector", "u238"),
                zr_frac=fp.get("zr_frac", ZR_METAL_FRAC),
                recipe=resolve_recipe(fp),
            )
        elif kind == "blanket":
            mats[mat_key] = make_smeared_blanket(
                name=mat_key,
                u235_frac=fp["u235_frac"],
                recipe=resolve_recipe(fp),
            )

    return openmc.Materials(list(mats.values())), mats


def make_material_from_composition(
    name: str,
    composition_g: dict,    # {nuclide: grams} — full smeared assembly composition
    rho: float = RHO_INTERIOR,
    depletable: bool = True,
) -> openmc.Material:
    """
    Build an OpenMC Material from a {nuclide: grams} composition dict.

    Used by cycle_runner to construct per-assembly materials for each cycle
    from the composition evolved in the previous cycle's depletion.

    The density is kept fixed at rho (default RHO_INTERIOR, fresh U-10Zr driver)
    consistent with the Tier-1 smeared material approximation.  Nuclide
    amounts are passed as relative proportions; OpenMC normalises internally.
    """
    m = openmc.Material(name=name)
    for nuc, grams in composition_g.items():
        if grams > 0.0:
            m.add_nuclide(nuc, grams)   # OpenMC normalises to sum = 1
    m.set_density("g/cm3", rho)
    m.depletable = depletable
    m.temperature = T_INTERIOR_K
    return m


def make_materials(internal_feed: str = "NU", radial_feed: str = "NU"):
    """Legacy entry point: build materials from hardcoded constants + NU/DU blanket choice."""
    mats = _make_structural_mats()

    mats["fuel_driver_inner_homog"] = make_smeared_driver("fuel_driver_inner_homog", PU_HM_INNER)
    mats["fuel_driver_outer_homog"] = make_smeared_driver("fuel_driver_outer_homog", PU_HM_OUTER)

    int_u235 = BLANKET_FEEDS[internal_feed]
    rad_u235 = BLANKET_FEEDS[radial_feed]
    mats["fuel_blanket_internal_homog"] = make_smeared_blanket(
        f"fuel_blanket_internal_{internal_feed.lower()}", int_u235
    )
    mats["fuel_blanket_radial_homog"] = make_smeared_blanket(
        f"fuel_blanket_radial_{radial_feed.lower()}", rad_u235
    )

    return openmc.Materials(list(mats.values())), mats