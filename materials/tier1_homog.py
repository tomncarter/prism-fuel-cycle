# materials/tier1_homog.py
import openmc

# -----------------------------
# Tunable knobs (Tier-1)
# -----------------------------
# Blanket feed options (U-235 atom fractions)
NU_U235_FRAC = 0.00711   # Natural uranium
DU_U235_FRAC = 0.0020    # Depleted uranium (0.20 wt% tails assay, atom-fraction approx)

BLANKET_FEEDS = {"NU": NU_U235_FRAC, "DU": DU_U235_FRAC}

# Reactor-grade Pu isotopic vector (atom fractions within Pu)
PU_VEC = {
    "Pu239": 0.55,
    "Pu240": 0.25,
    "Pu241": 0.12,
    "Pu242": 0.08,
}

# Fuel metal composition (atom fractions within fuel metal portion)
# We'll set: Pu_total (within heavy metal) and Zr fraction within metal
ZR_METAL_FRAC = 0.10  # ~10% Zr in U-Pu-Zr (simple surrogate)

# Smearing (volume-ish) of the assembly interior region:
# interior = fuel metal + sodium + structure
SMEAR_FUEL_METAL = 0.35
SMEAR_SODIUM     = 0.55
SMEAR_STEEL      = 0.10

# Effective density for the smeared interior region (g/cm3)
# This is deliberately lower than metal fuel density due to Na voiding in smear.
RHO_INTERIOR = 5.0

# Pu/HM targets for inner/outer zones (start here, then tune)
PU_HM_INNER = 0.20
PU_HM_OUTER = 0.17

# Single scalar multiplier to tune criticality quickly
# Increase -> higher k. Start at 1.0 then try 1.1, 1.2, etc.
FISSILE_SCALE = 0.923

# -----------------------------
# Vessel / pool region smearing
# -----------------------------
# Volume fraction of steel (Fe/Cr) in each region; remainder is sodium.
# Increase -> more neutron absorption/scattering in that region.
IHX_STEEL_VOL_FRAC          = 0.30   # IHX tube bundles + structure in radial pool
UPPER_PLENUM_STEEL_VOL_FRAC = 0.15   # lighter above core (fission gas plenums, light structure)
LOWER_PLENUM_STEEL_VOL_FRAC = 0.35   # heavier below core (pumps, lower internals)

_RHO_NA    = 0.92   # sodium density (g/cm3)
_RHO_STEEL = 7.80   # HT9/steel density (g/cm3)


def _normalize(d):
    s = sum(d.values())
    return {k: v / s for k, v in d.items()}


def make_smeared_driver(name: str, pu_hm: float) -> openmc.Material:
    """
    Smeared interior material approximating: (fuel metal + Na + steel).
    Fuel metal is U-Pu-Zr where Pu is reactor-grade vector.
    """
    pu_hm = max(0.0, min(0.35, pu_hm * FISSILE_SCALE))  # clamp for sanity
    u_hm = 1.0 - pu_hm

    # Build fuel metal atom fractions
    # total metal atoms: (HM portion) + Zr
    hm_frac_in_metal = 1.0 - ZR_METAL_FRAC
    pu_total_in_metal = hm_frac_in_metal * pu_hm
    u_total_in_metal  = hm_frac_in_metal * u_hm
    zr_in_metal       = ZR_METAL_FRAC

    pu_vec = _normalize(PU_VEC)

    m = openmc.Material(name=name)

    # Smear components: scale nuclide amounts by smear fractions.
    # This is not a strict volume-fraction model, but it is a stable Tier-1 surrogate.
    # Fuel metal contribution:
    m.add_nuclide("U238", SMEAR_FUEL_METAL * u_total_in_metal)

    for iso, frac in pu_vec.items():
        m.add_nuclide(iso, SMEAR_FUEL_METAL * pu_total_in_metal * frac)

    m.add_element("Zr", SMEAR_FUEL_METAL * zr_in_metal)

    # Sodium contribution:
    m.add_element("Na", SMEAR_SODIUM * 1.0)

    # Steel contribution (simple Fe/Cr surrogate):
    m.add_element("Fe", SMEAR_STEEL * 0.88)
    m.add_element("Cr", SMEAR_STEEL * 0.12)

    m.set_density("g/cm3", RHO_INTERIOR)
    return m


def make_smeared_blanket(name: str, u235_frac: float) -> openmc.Material:
    """
    Smeared blanket interior: NU or DU metal + sodium + steel, no Pu.
    u235_frac: atom fraction of U-235 (remainder is U-238).
    """
    u238_frac = 1.0 - u235_frac

    m = openmc.Material(name=name)

    # Uranium metal contribution
    m.add_nuclide("U235", SMEAR_FUEL_METAL * u235_frac)
    m.add_nuclide("U238", SMEAR_FUEL_METAL * u238_frac)

    # Sodium coolant contribution
    m.add_element("Na", SMEAR_SODIUM * 1.0)

    # Steel structure contribution
    m.add_element("Fe", SMEAR_STEEL * 0.88)
    m.add_element("Cr", SMEAR_STEEL * 0.12)

    m.set_density("g/cm3", RHO_INTERIOR)
    return m


def make_blended_na_steel(name: str, steel_vol_frac: float) -> openmc.Material:
    """
    Homogenized sodium + steel (Fe/Cr) region.
    steel_vol_frac: volume fraction of steel; remainder is sodium.
    Composition amounts are volume-fraction weighted (OpenMC normalises internally).
    Density is volume-averaged.
    """
    na_frac = 1.0 - steel_vol_frac
    m = openmc.Material(name=name)
    m.add_element("Na", na_frac)
    m.add_element("Fe", steel_vol_frac * 0.88)
    m.add_element("Cr", steel_vol_frac * 0.12)
    m.set_density("g/cm3", na_frac * _RHO_NA + steel_vol_frac * _RHO_STEEL)
    return m


def make_materials(internal_feed: str = "NU", radial_feed: str = "NU"):
    mats = {}

    # Sodium coolant (bulk)
    sodium = openmc.Material(name="sodium")
    sodium.add_element("Na", 1.0)
    sodium.set_density("g/cm3", 0.92)
    sodium.depletable = False
    mats["sodium"] = sodium

    # HT9 / steel surrogate
    ht9 = openmc.Material(name="ht9")
    ht9.add_element("Fe", 0.88)
    ht9.add_element("Cr", 0.12)
    ht9.set_density("g/cm3", 7.8)
    ht9.depletable = False
    mats["ht9"] = ht9

    # Reflector/shield surrogates (keep simple for now)
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

    # Absorber surrogates
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

    # Driver zones (smeared U-Pu-Zr + Na + steel)
    mats["fuel_driver_inner_homog"] = make_smeared_driver("fuel_driver_inner_homog", PU_HM_INNER)
    mats["fuel_driver_outer_homog"] = make_smeared_driver("fuel_driver_outer_homog", PU_HM_OUTER)

    # Blankets: NU or DU feed, selected per call
    int_u235 = BLANKET_FEEDS[internal_feed]
    rad_u235 = BLANKET_FEEDS[radial_feed]
    mats["fuel_blanket_internal_homog"] = make_smeared_blanket(
        f"fuel_blanket_internal_{internal_feed.lower()}", int_u235
    )
    mats["fuel_blanket_radial_homog"] = make_smeared_blanket(
        f"fuel_blanket_radial_{radial_feed.lower()}", rad_u235
    )

    # Vessel / pool regions (Na + steel blends) — not depletable
    for key, name, frac in [
        ("pool_ihx",     "pool_ihx",     IHX_STEEL_VOL_FRAC),
        ("upper_plenum", "upper_plenum", UPPER_PLENUM_STEEL_VOL_FRAC),
        ("lower_plenum", "lower_plenum", LOWER_PLENUM_STEEL_VOL_FRAC),
    ]:
        m = make_blended_na_steel(name, frac)
        m.depletable = False
        mats[key] = m

    return openmc.Materials(list(mats.values())), mats