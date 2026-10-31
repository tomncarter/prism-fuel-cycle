# fuel_cycle/nfrc.py
"""
NFRC electrometallurgical processing model and fresh fuel fabrication.

Pure Python — no OpenMC dependency.  All compositions are {nuclide: grams} dicts
where nuclide names follow OpenMC conventions (e.g. "Pu239", "Am241", "Zr90").

Electrorefiner physics (simplified):
  - Spent fuel anodes dissolve in molten LiCl-KCl salt.
  - TRU elements (Np, Pu, Am, Cm, ...) are co-deposited at the cathode
    at recovery fraction tru_recovery.  Clean Pu separation is electrochemically
    impossible in the presence of uranium chloride — all TRU is co-recovered.
    This is simultaneously the proliferation resistance feature and the MA
    transmutation mechanism of the IFR/PRISM cycle.
  - Uranium is recovered at u_recovery efficiency and set aside (either
    recycled into fresh fuel or stored as clean DU equivalent).
  - Fission products and structural material (Zr from cladding/fuel, Na, Fe, Cr)
    partition to the salt or metallic waste streams.
"""
from __future__ import annotations

import re
from dataclasses import dataclass

# ---------------------------------------------------------------------------
# Element classification
# ---------------------------------------------------------------------------

URANIUM_ELEMENTS: frozenset[str] = frozenset({"U"})

# TRU = all actinides heavier than uranium relevant to PRISM fuel
TRU_ELEMENTS: frozenset[str] = frozenset({"Np", "Pu", "Am", "Cm", "Bk", "Cf"})

HEAVY_METAL_ELEMENTS: frozenset[str] = URANIUM_ELEMENTS | TRU_ELEMENTS

# Comprehensive list of actinide nuclides tracked through the fuel cycle
ACTINIDE_NUCLIDES: list[str] = [
    "U232", "U233", "U234", "U235", "U236", "U238",
    "Np236", "Np237", "Np238", "Np239",
    "Pu236", "Pu237", "Pu238", "Pu239", "Pu240", "Pu241", "Pu242", "Pu244",
    "Am241", "Am242", "Am242_m1", "Am243",
    "Cm241", "Cm242", "Cm243", "Cm244", "Cm245", "Cm246", "Cm247", "Cm248",
]

# Out-of-core decay steps applied during cooling: parent -> (daughter, half-life in years)
DECAY_CHAINS: dict[str, tuple[str, float]] = {
    "Np239": ("Pu239", 2.356 / 365.25),
    "Cm242": ("Pu238", 162.8 / 365.25),
    "Pu241": ("Am241", 14.29),
    "Cm244": ("Pu240", 18.11),
    "Pu238": ("U234", 87.7),
}

# Fission products of interest for decay-heat and safeguards tracking
FP_NUCLIDES_OF_INTEREST: list[str] = [
    "Cs137", "Sr90", "Nd148", "Xe135", "Sm149", "Tc99", "I129",
]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def nuclide_element(nuc: str) -> str:
    """
    Extract the element symbol from an OpenMC nuclide name.
    Examples: "Pu239" -> "Pu",  "Am241_m1" -> "Am",  "Zr90" -> "Zr"
    """
    return re.match(r"[A-Za-z]+", nuc).group(0)


def split_hm(composition: dict) -> tuple[dict, dict]:
    """
    Partition {nuclide: grams} into (heavy_metal, non_hm).
    heavy_metal contains U and TRU nuclides; non_hm contains everything else
    (fission products, Zr, structural Na/Fe/Cr).
    """
    hm: dict[str, float] = {}
    non_hm: dict[str, float] = {}
    for nuc, grams in composition.items():
        if nuclide_element(nuc) in HEAVY_METAL_ELEMENTS:
            hm[nuc] = grams
        else:
            non_hm[nuc] = grams
    return hm, non_hm


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

@dataclass
class NFRCConfig:
    """Operating parameters for one NFRC facility."""
    tru_recovery: float = 0.999   # fraction of TRU recovered to cathode product
    u_recovery: float = 0.95      # fraction of U recovered (remainder to waste salt)
    cooling_cycles: int = 1       # in-vessel cooling cycles before NFRC processing
    du_u235_frac: float = 0.0020  # DU U-235 atom fraction (0.20 wt% tails assay)
    external_feed: str = "reactor_grade"  # Pu vector for startup / shortfall TRU
    recycle_u: bool = True        # recovered U forms the fresh-fuel matrix, topped up with DU


# ---------------------------------------------------------------------------
# Core processing functions
# ---------------------------------------------------------------------------

def electrorefine(
    hm_composition: dict,
    tru_recovery: float = 0.999,
    u_recovery: float = 0.95,
) -> tuple[dict, dict, dict]:
    """
    Simplified electrorefiner model operating on HM nuclides only.

    Args:
        hm_composition: {nuclide: grams} — U and TRU nuclides from spent fuel
        tru_recovery:   fraction of each TRU nuclide recovered to cathode product
        u_recovery:     fraction of each U nuclide recovered to cathode

    Returns:
        tru_product:  {nuclide: grams} — cathode product (U + TRU co-deposit)
        u_product:    {nuclide: grams} — separately recovered uranium (EG23 mode)
        hm_waste:     {nuclide: grams} — unrecovered HM remaining in salt waste
    """
    tru_product: dict[str, float] = {}
    u_product: dict[str, float] = {}
    hm_waste: dict[str, float] = {}

    for nuc, grams in hm_composition.items():
        elem = nuclide_element(nuc)
        if elem in TRU_ELEMENTS:
            tru_product[nuc] = grams * tru_recovery
            hm_waste[nuc] = grams * (1.0 - tru_recovery)
        else:  # Uranium
            u_product[nuc] = grams * u_recovery
            hm_waste[nuc] = grams * (1.0 - u_recovery)

    return tru_product, u_product, hm_waste


def decay(hm_composition: dict, years: float) -> dict:
    """
    Out-of-core decay of the short actinide chains that matter for the recycled
    TRU vector (Pu-241 -> Am-241 above all).  Single parent -> daughter steps,
    applied in table order; daughter mass scales by the mass-number ratio, so
    alpha decays lose the He mass.
    """
    out = dict(hm_composition)
    for parent, (daughter, half_life_y) in DECAY_CHAINS.items():
        g = out.get(parent, 0.0)
        if g <= 0.0:
            continue
        decayed = g * (1.0 - 0.5 ** (years / half_life_y))
        out[parent] = g - decayed
        out[daughter] = out.get(daughter, 0.0) + (
            decayed * atomic_mass(daughter) / atomic_mass(parent)
        )
    return out


@dataclass
class FabricationResult:
    """One batch of fresh fuel and where its heavy metal came from."""
    hm_vec: dict           # {nuclide: atom fraction within HM}, same for every assembly
    from_tru_stock: dict   # {nuclide: grams} drawn from recovered TRU
    from_u_stock: dict     # {nuclide: grams} drawn from recovered U
    from_makeup_u: dict    # {nuclide: grams} of fresh DU/NU
    from_external: dict    # {nuclide: grams} of external fissile feed (startup / shortfall)

    @property
    def total_g(self) -> float:
        return sum(
            sum(d.values()) for d in (
                self.from_tru_stock, self.from_u_stock,
                self.from_makeup_u, self.from_external,
            )
        )


def fabricate(
    n_assemblies: int,
    hm_mol_per_assembly: float,
    tru_hm_frac: float,
    tru_stock: dict,
    u_stock: dict,
    makeup_u_vec: dict,
    external_vec: dict,
    recycle_u: bool = False,
) -> FabricationResult:
    """
    Fabricate a batch of identical fresh assemblies at a target TRU/HM atom
    fraction, drawing TRU from the recovered stock first and making up any
    shortfall from the external feed (RG-Pu by default).  The U matrix comes
    from recovered U when recycle_u is set, topped up with makeup U.

    tru_stock and u_stock are {nuclide: grams} and are reduced in place by
    what is drawn; whatever TRU is left over is the surplus.

    makeup_u_vec / external_vec: {nuclide: atom fraction} isotopic vectors.
    """
    total_mol = n_assemblies * hm_mol_per_assembly
    tru_mol = tru_hm_frac * total_mol
    u_mol = total_mol - tru_mol

    from_tru, tru_short = _draw(tru_stock, tru_mol)
    from_ext = _grams_from_vec(external_vec, tru_short)

    from_u: dict = {}
    u_short = u_mol
    if recycle_u:
        from_u, u_short = _draw(u_stock, u_mol)
    from_makeup = _grams_from_vec(makeup_u_vec, u_short)

    moles: dict[str, float] = {}
    for part in (from_tru, from_ext, from_u, from_makeup):
        for nuc, g in part.items():
            moles[nuc] = moles.get(nuc, 0.0) + g / atomic_mass(nuc)
    total = sum(moles.values())
    hm_vec = {nuc: mol / total for nuc, mol in moles.items()}

    return FabricationResult(
        hm_vec=hm_vec,
        from_tru_stock=from_tru,
        from_u_stock=from_u,
        from_makeup_u=from_makeup,
        from_external=from_ext,
    )


def mass_number(nuc: str) -> int:
    """Mass number from an OpenMC nuclide name: "Am242_m1" -> 242."""
    return int(re.match(r"[A-Za-z]+(\d+)", nuc).group(1))


def atomic_mass(nuc: str) -> float:
    """
    Atomic mass (g/mol).  Uses OpenMC's mass table when available so gram
    ledgers here agree with Material.get_mass(); falls back to the mass number.
    """
    if _atomic_mass_table is not None:
        try:
            return _atomic_mass_table(nuc)
        except (KeyError, ValueError):
            pass
    return float(mass_number(nuc))


try:   # optional: keeps this module importable without OpenMC
    from openmc.data import atomic_mass as _atomic_mass_table
except ImportError:
    _atomic_mass_table = None


def _draw(stock: dict, mol: float) -> tuple[dict, float]:
    """
    Take up to `mol` moles from a {nuclide: grams} stock, pro rata across its
    isotopes.  Mutates stock; returns (grams drawn, moles still needed).
    """
    stock_mol = sum(g / atomic_mass(n) for n, g in stock.items())
    if stock_mol <= 0.0 or mol <= 0.0:
        return {}, max(mol, 0.0)
    take = min(1.0, mol / stock_mol)
    drawn = {n: g * take for n, g in stock.items() if g > 0.0}
    for n, g in drawn.items():
        stock[n] -= g
    # A partial draw meets the demand exactly; don't let round-off book a shortfall
    return drawn, (mol - stock_mol) if take == 1.0 else 0.0


def _grams_from_vec(vec: dict, mol: float) -> dict:
    """{nuclide: grams} for `mol` moles of material with atom-fraction vector vec."""
    if mol <= 0.0:
        return {}
    s = sum(vec.values())
    return {n: mol * f / s * atomic_mass(n) for n, f in vec.items() if f > 0.0}


# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------

def compute_tru_metrics(hm_composition: dict) -> tuple[float, float]:
    """
    Compute total in-core TRU mass and MA/TRU weight fraction.

    MA (minor actinides) = Np + Am + Cm.
    Returns (tru_mass_g, ma_tru_frac).
    """
    tru_g = sum(
        g for nuc, g in hm_composition.items()
        if nuclide_element(nuc) in TRU_ELEMENTS
    )
    pu_g = sum(
        g for nuc, g in hm_composition.items()
        if nuclide_element(nuc) == "Pu"
    )
    ma_g = tru_g - pu_g
    ma_frac = ma_g / tru_g if tru_g > 0.0 else 0.0
    return tru_g, ma_frac


def aggregate_tru_feed(compositions: list[dict]) -> dict:
    """
    Sum TRU nuclides across a list of assembly HM compositions.
    Used to pool the TRU feed from multiple discharged assemblies before
    blending into fresh fuel slugs.
    """
    pooled: dict[str, float] = {}
    for comp in compositions:
        for nuc, grams in comp.items():
            if nuclide_element(nuc) in TRU_ELEMENTS:
                pooled[nuc] = pooled.get(nuc, 0.0) + grams
    return pooled
