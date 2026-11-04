# materials/tier1_homog.py
import openmc

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

# Volume fractions of the assembly interior (inside the duct; the duct wall
# and inter-assembly gap are modelled explicitly in UniverseFactory).
# Typical metal-fuel SFR driver: fuel column inside the cladding at 75% smear
# density, HT9 cladding + wire wrap, and coolant + bond sodium filling the rest.
SMEAR_FUEL_METAL = 0.27   # as-fabricated fuel slug (~0.36 clad-inner x 0.75 smear)
SMEAR_SODIUM     = 0.53   # coolant + bond sodium
SMEAR_STEEL      = 0.20   # HT9 cladding + wire wrap

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

    return make_smeared_fuel(name, hm_vec, zr_frac=zr_frac)


def make_smeared_blanket(name: str, u235_frac: float) -> openmc.Material:
    """
    Smeared blanket interior: NU or DU metal + sodium + steel, no Pu.
    u235_frac: atom fraction of U-235 (remainder is U-238).
    """
    return make_smeared_fuel(
        name, {"U235": u235_frac, "U238": 1.0 - u235_frac}, zr_frac=0.0
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


def make_smeared_fuel(name: str, hm_vec: dict, zr_frac: float = ZR_METAL_FRAC) -> openmc.Material:
    """
    Smeared interior (fuel metal + Na + steel) for an arbitrary heavy-metal
    isotopic vector.  Shared recipe behind drivers, blankets and the fresh
    assemblies fabricated by fuel_cycle from recycled U/TRU.

    Each component is taken at its theoretical density and weighted by its
    SMEAR_* volume fraction, so nuclide amounts are true atom/b-cm.

    hm_vec:  {nuclide: atom fraction within HM} (normalised here)
    zr_frac: Zr weight fraction of the fuel alloy (0 for blankets)
    """
    hm_vec = _normalize(hm_vec)
    na_av = _avogadro_b_cm()

    m = openmc.Material(name=name)

    # Fuel slug: HM and Zr partial densities of the alloy, scaled by volume fraction.
    rho_alloy = _fuel_alloy_density(zr_frac)
    m_hm = sum(frac * openmc.data.atomic_mass(iso) for iso, frac in hm_vec.items())
    n_hm = SMEAR_FUEL_METAL * rho_alloy * (1.0 - zr_frac) * na_av / m_hm
    for iso, frac in hm_vec.items():
        if frac > 0.0:
            m.add_nuclide(iso, n_hm * frac)

    if zr_frac > 0.0:
        n_zr = SMEAR_FUEL_METAL * rho_alloy * zr_frac * na_av / openmc.data.atomic_weight("Zr")
        m.add_element("Zr", n_zr)

    # Coolant + bond sodium:
    m.add_element("Na", SMEAR_SODIUM * _RHO_NA * na_av / openmc.data.atomic_weight("Na"))

    # Cladding + wire wrap (simple Fe/Cr surrogate):
    n_steel = SMEAR_STEEL * _steel_number_density()
    m.add_element("Fe", n_steel * 0.88)
    m.add_element("Cr", n_steel * 0.12)

    m.set_density("sum")
    m.temperature = T_INTERIOR_K
    return m


# Mass density of a fresh U-10Zr driver interior; the fixed density used when
# rebuilding materials from evolved {nuclide: grams} compositions.
RHO_INTERIOR = (SMEAR_FUEL_METAL * _fuel_alloy_density(ZR_METAL_FRAC)
                + SMEAR_SODIUM * _RHO_NA + SMEAR_STEEL * _RHO_STEEL)


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
            )
        elif kind == "blanket":
            mats[mat_key] = make_smeared_blanket(
                name=mat_key,
                u235_frac=fp["u235_frac"],
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