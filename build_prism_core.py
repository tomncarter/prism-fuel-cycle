# build_prism_core.py
from __future__ import annotations

import math
import csv
import os
from collections import Counter
from contextlib import contextmanager
from typing import Dict, Tuple, List

import openmc
import openmc.deplete

# -----------------------------
# Depletion / run parameters
# -----------------------------
CHAIN_FILE      = "/chain/chain_endfb80_fast.xml"
POWER_MWth      = 840.0          # thermal power
CAPACITY_FACTOR = 1.0            # fraction of full power; multiply to get effective power
POWER_W         = POWER_MWth * 1e6 * CAPACITY_FACTOR
TIMESTEP_DAYS   = 30             # length of each depletion step
N_CYCLES        = 36             # total steps (36 × 30 d = ~3 years)

# -----------------------------
# Run profiles (--profile NAME)
# -----------------------------
# Transport statistics + depletion schedule per profile.  "standard" reproduces
# the historical defaults and writes to runs/<loadout>/; every other profile
# writes to runs/<loadout>__<profile>/ so results never mix (or resume) across
# profiles.  Depletion "timesteps" are in days.
RUN_PROFILES = {
    # Pipeline check only: seconds-to-minutes, k-eff noise of several hundred pcm.
    "smoke": dict(particles=2_000, batches=25, inactive=10,
                  timesteps=[30, 30]),
    # Screening / early data: ~1-year cycle in 6 steps.  Short leading steps
    # let Xe-135 / Sm-149 saturate before the long steps.
    "quick": dict(particles=10_000, batches=60, inactive=20,
                  timesteps=[5, 25, 60, 90, 90, 90]),
    # Historical defaults (36 × 30 d).
    "standard": dict(particles=10_000, batches=80, inactive=20,
                     timesteps=[TIMESTEP_DAYS] * N_CYCLES),
    # Per-assembly depletion: half the quick histories, since depletion tallies
    # run at ~2,000 p/s.  Reload searches keep their own, tighter k-eff.
    "cycle": dict(particles=5_000, batches=60, inactive=20,
                  timesteps=[5, 25, 60, 90, 90, 90]),
    # T5.1 target statistics (≤15 pcm); expensive.
    "production": dict(particles=5_000_000, batches=600, inactive=100,
                       timesteps=[TIMESTEP_DAYS] * N_CYCLES),
}


@contextmanager
def _working_dir(path: str):
    """Temporarily change working directory."""
    orig = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(orig)

from geometry.params import AssemblyGeom, Prism3D, VesselGeom
from geometry.universe_factory import UniverseFactory
from materials.tier1_homog import make_materials, make_materials_from_loadout

from geometry.hexmap import build_hexlattice_rings_from_assignment, full_hex_coords, hex_distance, axial_to_xy
from core_maps.abr_like_assignment import abr_like_assignment
from core_maps.loader import load_axial, load_core

Axial = Tuple[int, int]


# -----------------------------
# Debug artifacts
# -----------------------------
def write_core_assignment_csv(
    assignment: Dict[Axial, str],
    max_ring: int,
    start_angle_deg: float,
    out_csv: str = "core_assignment.csv",
) -> None:
    """
    Write per-position mapping: axial coords -> ring -> index -> universe name.
    Index is within the OpenMC ring ordering (clockwise, starting at top).
    """
    # Build rings exactly as OpenMC expects (outer->inner + center)
    ring_names, _ = build_hexlattice_rings_from_assignment(
        assignment=assignment,
        max_ring=max_ring,
        default_universe="u_reflector",
        start_angle_deg=start_angle_deg,
    )

    # Reconstruct which coordinates correspond to each ring element in that ordering
    # We do this by regenerating the coordinate rings and sorting them the same way.
    from geometry.hexmap import ring_angle_sort  # imported here to avoid circulars

    coords = full_hex_coords(max_ring)
    by_ring: Dict[int, List[Axial]] = {}
    for c in coords:
        d = hex_distance(c)
        by_ring.setdefault(d, []).append(c)

    ordered_coords_per_ring: List[List[Axial]] = []
    for ring in range(max_ring, 0, -1):
        ordered_coords_per_ring.append(ring_angle_sort(by_ring[ring], start_angle_deg=start_angle_deg))
    ordered_coords_per_ring.append([(0, 0)])  # center

    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["q", "r", "ring", "index_in_ring", "universe"])
        for ring_idx, (coords_ring, names_ring) in enumerate(zip(ordered_coords_per_ring, ring_names), start=1):
            # ring_idx here counts outer->inner + center as final; we also write the true hex distance ring
            for i, (c, name) in enumerate(zip(coords_ring, names_ring)):
                q, r = c
                w.writerow([q, r, hex_distance(c), i, name])


def write_debug_plot(
    assignment: Dict[Axial, str],
    pitch_cm: float,
    max_ring: int = 11,
    out_png: str = "core_debug_map.png",
) -> None:
    """
    Writes a debug plot showing axial hex map. Uses matplotlib.
    """
    import matplotlib.pyplot as plt
    from matplotlib.patches import RegularPolygon

    # Conservative grayscale/hatch styling
    styles = {
        "u_inner_fuel":     dict(facecolor="white", edgecolor="black", linewidth=0.5),
        "u_outer_fuel":     dict(facecolor="0.85", edgecolor="black", linewidth=0.5),
        "u_reflector":      dict(facecolor="0.65", edgecolor="black", linewidth=0.5),
        "u_shield":         dict(facecolor="white", edgecolor="black", linewidth=0.5, hatch="///"),
        "u_control":        dict(facecolor="black", edgecolor="black", linewidth=0.5),
        "u_ult_shutdown":   dict(facecolor="0.3", edgecolor="black", linewidth=0.5),
        "u_gem":            dict(facecolor="0.9", edgecolor="black", linewidth=0.5, hatch=".."),
        "u_source":         dict(facecolor="0.9", edgecolor="black", linewidth=0.5, hatch="xx"),
        "u_internal_blanket": dict(facecolor="0.75", edgecolor="black", linewidth=0.5),
        "u_radial_blanket":   dict(facecolor="0.75", edgecolor="black", linewidth=0.5),
    }

    coords = full_hex_coords(max_ring)

    fig = plt.figure(figsize=(18, 18))
    ax = plt.gca()
    ax.set_aspect("equal")
    ax.axis("off")

    # pointy-top hex; fill tightly (0.999 leaves a hairline gap to keep edges visible)
    hex_radius = pitch_cm / math.sqrt(3.0) * 0.999

    xs, ys = [], []
    for (q, r) in coords:
        un = assignment.get((q, r), "u_reflector")
        x, y = axial_to_xy(q, r, pitch=pitch_cm)
        xs.append(x); ys.append(y)

        kw = styles.get(un, dict(facecolor="white", edgecolor="black", linewidth=0.5))
        poly = RegularPolygon(
            (x, y),
            numVertices=6,
            radius=hex_radius,
            orientation=math.radians(30),
            **kw
        )
        ax.add_patch(poly)

    # Double-width outline on center hex
    cx, cy = axial_to_xy(0, 0, pitch=pitch_cm)
    ax.add_patch(RegularPolygon(
        (cx, cy), numVertices=6, radius=hex_radius,
        orientation=math.radians(30),
        facecolor="none", edgecolor="black", linewidth=2.0,
    ))

    ax.set_xlim(min(xs) - pitch_cm, max(xs) + pitch_cm)
    ax.set_ylim(min(ys) - pitch_cm, max(ys) + pitch_cm)
    ax.set_title("Core Debug Map (axial coords -> universe assignment)")
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close(fig)

def write_reference_plot(
    assignment,
    max_ring: int,
    pitch_cm: float,
    mode: str = "qr",  # "qr", "ring", "index", "both"
    out_png: str = "core_reference.png",
):
    import matplotlib.pyplot as plt
    from matplotlib.patches import RegularPolygon

    from geometry.hexmap import ring_angle_sort

    coords = full_hex_coords(max_ring)

    # Build OpenMC ordering to get ring + index
    ring_names, _ = build_hexlattice_rings_from_assignment(
        assignment=assignment,
        max_ring=max_ring,
        default_universe="u_reflector",
        start_angle_deg=90.0,
    )

    # Build coordinate ordering per ring
    by_ring = {}
    for c in coords:
        by_ring.setdefault(hex_distance(c), []).append(c)

    ordered_coords = []
    for ring in range(max_ring, 0, -1):
        ordered_coords.append(ring_angle_sort(by_ring[ring], start_angle_deg=90.0))
    ordered_coords.append([(0, 0)])

    fig = plt.figure(figsize=(30, 30))
    ax = plt.gca()
    ax.set_aspect("equal")
    ax.axis("off")

    hex_radius = pitch_cm / math.sqrt(3.0) * 0.999

    xs, ys = [], []

    for ring_idx, coords_ring in enumerate(ordered_coords):
        for i, (q, r) in enumerate(coords_ring):
            x, y = axial_to_xy(q, r, pitch=pitch_cm)
            xs.append(x); ys.append(y)

            poly = RegularPolygon(
                (x, y),
                numVertices=6,
                radius=hex_radius,
                orientation=math.radians(30),
                facecolor="white",
                edgecolor="black",
                linewidth=0.7,
            )
            ax.add_patch(poly)

            if mode == "qr":
                label = f"{q},{r}"
            elif mode == "ring":
                label = f"{hex_distance((q,r))}"
            elif mode == "index":
                label = f"{i}"
            elif mode == "both":
                label = f"{hex_distance((q,r))}:{i}"
            else:
                label = ""

            ax.text(
                x, y,
                label,
                ha="center", va="center",
                fontsize=14,
            )

    # Double-width outline on center hex
    cx, cy = axial_to_xy(0, 0, pitch=pitch_cm)
    ax.add_patch(RegularPolygon(
        (cx, cy), numVertices=6, radius=hex_radius,
        orientation=math.radians(30),
        facecolor="none", edgecolor="black", linewidth=2.5,
    ))

    margin = hex_radius * 1.1
    ax.set_xlim(min(xs) - margin, max(xs) + margin)
    ax.set_ylim(min(ys) - margin, max(ys) + margin)
    ax.set_title(f"Core reference map ({mode})", fontsize=20)
    plt.tight_layout(pad=0.3)
    plt.savefig(out_png, dpi=300)
    plt.close(fig)

def write_loading_map_with_legend(
    assignment: Dict[Axial, str],
    max_ring: int,
    pitch_cm: float,
    out_png: str = "core_loading_map.png",
) -> None:
    """
    Loading map by universe type with legend. Grayscale + hatches.
    """
    import matplotlib.pyplot as plt
    from matplotlib.patches import RegularPolygon, Patch

    coords = full_hex_coords(max_ring)

    styles = {
        "u_inner_fuel":       dict(facecolor="white", edgecolor="black", linewidth=0.5),
        "u_outer_fuel":       dict(facecolor="0.85", edgecolor="black", linewidth=0.5),
        "u_reflector":        dict(facecolor="0.65", edgecolor="black", linewidth=0.5),
        "u_shield":           dict(facecolor="white", edgecolor="black", linewidth=0.5, hatch="///"),
        "u_control":          dict(facecolor="black", edgecolor="black", linewidth=0.5),
        "u_ult_shutdown":     dict(facecolor="0.3", edgecolor="black", linewidth=0.5),
        "u_gem":              dict(facecolor="0.9", edgecolor="black", linewidth=0.5, hatch=".."),
        "u_source":           dict(facecolor="0.9", edgecolor="black", linewidth=0.5, hatch="xx"),
        "u_internal_blanket": dict(facecolor="0.75", edgecolor="black", linewidth=0.5, hatch=".."),
        "u_radial_blanket":   dict(facecolor="0.75", edgecolor="black", linewidth=0.5, hatch="xx"),
    }

    # Nice labels for legend
    labels = {
        "u_inner_fuel": "Inner fuel",
        "u_outer_fuel": "Outer fuel",
        "u_reflector": "Reflector",
        "u_shield": "Shield",
        "u_control": "Control",
        "u_ult_shutdown": "Ultimate shutdown",
        "u_gem": "GEM",
        "u_source": "Source",
        "u_internal_blanket": "Internal blanket",
        "u_radial_blanket": "Radial blanket",
    }

    fig = plt.figure(figsize=(20, 20))
    ax = plt.gca()
    ax.set_aspect("equal")
    ax.axis("off")

    hex_radius = pitch_cm / math.sqrt(3.0) * 0.999

    xs, ys = [], []
    used = set()

    for (q, r) in coords:
        un = assignment.get((q, r), "u_reflector")
        used.add(un)

        x, y = axial_to_xy(q, r, pitch=pitch_cm)
        xs.append(x); ys.append(y)

        kw = styles.get(un, dict(facecolor="white", edgecolor="black", linewidth=0.5))
        poly = RegularPolygon(
            (x, y),
            numVertices=6,
            radius=hex_radius,
            orientation=math.radians(30),
            **kw
        )
        ax.add_patch(poly)

    # Double-width outline on center hex
    cx, cy = axial_to_xy(0, 0, pitch=pitch_cm)
    ax.add_patch(RegularPolygon(
        (cx, cy), numVertices=6, radius=hex_radius,
        orientation=math.radians(30),
        facecolor="none", edgecolor="black", linewidth=2.0,
    ))

    # Legend entries in a stable order
    order = [
        "u_inner_fuel", "u_outer_fuel", "u_internal_blanket", "u_radial_blanket",
        "u_reflector", "u_shield", "u_control", "u_ult_shutdown", "u_gem", "u_source",
    ]
    legend_handles = []
    for key in order:
        if key in used and key in styles:
            s = styles[key]
            legend_handles.append(
                Patch(
                    facecolor=s.get("facecolor", "white"),
                    edgecolor=s.get("edgecolor", "black"),
                    hatch=s.get("hatch", None),
                    label=labels.get(key, key),
                )
            )

    ax.set_xlim(min(xs) - pitch_cm, max(xs) + pitch_cm)
    ax.set_ylim(min(ys) - pitch_cm, max(ys) + pitch_cm)
    ax.legend(handles=legend_handles, loc="upper right", frameon=True, fontsize=14)
    ax.set_title("Core loading map (by universe type)", fontsize=16)
    plt.tight_layout()
    plt.savefig(out_png, dpi=250)
    plt.close(fig)


# -----------------------------
# Build OpenMC model (3D)
# -----------------------------
def _entropy_mesh(vg: VesselGeom, axial: Prism3D) -> openmc.RegularMesh:
    """Shannon entropy mesh over the active core, for checking source convergence."""
    mesh = openmc.RegularMesh()
    mesh.dimension = (9, 9, 6)
    mesh.lower_left = (-vg.lattice_bound_r_cm, -vg.lattice_bound_r_cm, -axial.active_height_cm / 2.0)
    mesh.upper_right = (vg.lattice_bound_r_cm, vg.lattice_bound_r_cm, axial.active_height_cm / 2.0)
    return mesh


def build_model(loadout_path: str, profile: str = "standard"):
    """
    Build an OpenMC model from a core loadout JSON file.

    The loadout specifies which assembly type (and therefore which material)
    occupies each lattice position, and may set "active_height_cm".  All other
    geometry constants (pitch, duct dimensions, vessel radii, etc.) remain in Python.
    Transport statistics come from RUN_PROFILES[profile].
    """
    prof = RUN_PROFILES[profile]
    # Load core layout
    positions, asm_types = load_core(loadout_path)

    # Materials derived from loadout
    materials, mats_dict = make_materials_from_loadout(asm_types)

    # Geometry params
    ag = AssemblyGeom()
    axial = load_axial(loadout_path)
    vg = VesselGeom()

    # Universes
    uf = UniverseFactory(ag, mats_dict)
    univs = uf.build_all()

    # Flat assignment dict: (q, r) -> universe name  (used by lattice builder + debug plots)
    assignment = {(p.q, p.r): p.universe for p in positions}
    max_ring = max(hex_distance((p.q, p.r)) for p in positions)

    # Analytic material volumes (needed for depletion normalization).
    # Interior volume per assembly = hex area of duct inner bore × active height.
    # Tier-1: all assemblies of the same type share one material, so the volume
    # is the sum across all assemblies of that type.
    _f2f  = ag.duct_inner_flat_to_flat_cm
    _area = (math.sqrt(3) / 2.0) * _f2f ** 2   # regular hexagon, flat-to-flat F: A = (√3/2)F²
    _vol_per_assy = _area * axial.active_height_cm

    # Count positions by universe to compute aggregate volumes
    _univ_counts: Counter = Counter(p.universe for p in positions)
    _fuel_vol_map = {
        "fuel_driver_inner_homog":    "u_inner_fuel",
        "fuel_driver_outer_homog":    "u_outer_fuel",
        "fuel_blanket_internal_homog":"u_internal_blanket",
        "fuel_blanket_radial_homog":  "u_radial_blanket",
    }
    for mat_key, univ_name in _fuel_vol_map.items():
        if mat_key in mats_dict:
            mats_dict[mat_key].volume = _univ_counts[univ_name] * _vol_per_assy

    # Debug artifacts
    write_reference_plot(
        assignment, max_ring=max_ring, pitch_cm=ag.assembly_pitch_cm,
        mode="qr", out_png="core_qr_labels.png",
    )
    write_reference_plot(
        assignment, max_ring=max_ring, pitch_cm=ag.assembly_pitch_cm,
        mode="both", out_png="core_ring_index_map.png",
    )
    write_loading_map_with_legend(
        assignment, max_ring=max_ring, pitch_cm=ag.assembly_pitch_cm,
        out_png="core_loading_map.png",
    )
    write_core_assignment_csv(
        assignment, max_ring=max_ring, start_angle_deg=90.0,
        out_csv="core_assignment.csv",
    )

    # Build lattice ring lists from position assignment
    ring_names, _center = build_hexlattice_rings_from_assignment(
        assignment=assignment,
        max_ring=max_ring,
        default_universe="u_reflector",
        start_angle_deg=90.0,
    )

    ring_lengths = [len(r) for r in ring_names]
    print("Hex ring lengths:", ring_lengths)

    # Hex lattice
    lattice = openmc.HexLattice(name="core_hex_lattice")
    lattice.center = (0.0, 0.0)
    lattice.pitch = (ag.assembly_pitch_cm,)
    lattice.universes = [[univs[name] for name in ring] for ring in ring_names]

    # Outer universe (outside lattice): sodium, UNBOUNDED
    outer_u = openmc.Universe(name="outer_sodium")
    outer_u.add_cell(openmc.Cell(fill=mats_dict["sodium"]))
    lattice.outer = outer_u

    # -----------------------------
    # Vessel geometry (radial + axial layering)
    # -----------------------------
    H = axial.active_height_cm

    # Radial surfaces
    r_lat  = openmc.ZCylinder(r=vg.lattice_bound_r_cm)
    r_pool = openmc.ZCylinder(r=vg.ihx_pool_outer_r_cm)
    r_ves  = openmc.ZCylinder(r=vg.vessel_outer_r_cm, boundary_type="vacuum")

    # Axial surfaces
    z_bot_vac  = openmc.ZPlane(z0=-(H/2 + vg.lower_plenum_h_cm + vg.end_plate_cm), boundary_type="vacuum")
    z_bot_plen = openmc.ZPlane(z0=-(H/2 + vg.lower_plenum_h_cm))
    z_act_bot  = openmc.ZPlane(z0=-H/2)
    z_act_top  = openmc.ZPlane(z0=+H/2)
    z_top_plen = openmc.ZPlane(z0=+(H/2 + vg.upper_plenum_h_cm))
    z_top_vac  = openmc.ZPlane(z0=+(H/2 + vg.upper_plenum_h_cm + vg.end_plate_cm), boundary_type="vacuum")

    ht9          = mats_dict["ht9"]
    ihx_mat      = mats_dict["pool_ihx"]
    lo_plen_mat  = mats_dict["lower_plenum"]
    hi_plen_mat  = mats_dict["upper_plenum"]

    # Root cells — 9 cells covering all space within the vessel boundary
    root_cells = [
        # Active zone: lattice core | IHX annulus | vessel wall
        openmc.Cell(name="active_core",        fill=lattice,     region= -r_lat  & +z_act_bot & -z_act_top),
        openmc.Cell(name="ihx_pool",           fill=ihx_mat,     region= +r_lat  & -r_pool & +z_act_bot & -z_act_top),
        openmc.Cell(name="vessel_wall_active", fill=ht9,         region= +r_pool & -r_ves  & +z_act_bot & -z_act_top),
        # Lower plenum slab
        openmc.Cell(name="lower_plenum",       fill=lo_plen_mat, region= -r_pool & +z_bot_plen & -z_act_bot),
        openmc.Cell(name="vessel_wall_lower",  fill=ht9,         region= +r_pool & -r_ves  & +z_bot_plen & -z_act_bot),
        # Upper plenum slab
        openmc.Cell(name="upper_plenum",       fill=hi_plen_mat, region= -r_pool & +z_act_top & -z_top_plen),
        openmc.Cell(name="vessel_wall_upper",  fill=ht9,         region= +r_pool & -r_ves  & +z_act_top & -z_top_plen),
        # End plates (full r_ves width, HT9)
        openmc.Cell(name="lower_end_plate",    fill=ht9,         region= -r_ves  & +z_bot_vac & -z_bot_plen),
        openmc.Cell(name="upper_end_plate",    fill=ht9,         region= -r_ves  & +z_top_plen & -z_top_vac),
    ]

    root_univ = openmc.Universe(cells=root_cells)
    geometry = openmc.Geometry(root_univ)

    # Settings
    settings = openmc.Settings()
    settings.run_mode = "eigenvalue"
    settings.batches = prof["batches"]
    settings.inactive = prof["inactive"]
    settings.particles = prof["particles"]
    settings.temperature = {"method": "nearest"}
    settings.statepoint = {"batches": [settings.batches]}
    settings.statepoint["filename"] = "statepoint_boc.h5"
    settings.entropy_mesh = _entropy_mesh(vg, axial)

    # Source (0.15.3): keep only_fissionable=True
    bounds = [-250, -250, -axial.active_height_cm / 2.0, 250, 250, axial.active_height_cm / 2.0]
    settings.source = openmc.IndependentSource(
        space=openmc.stats.Box(bounds[:3], bounds[3:]),
        constraints={"fissionable": True},
    )

    # -----------------------------
    # Tallies (minimal QA set)
    # -----------------------------
    tallies = openmc.Tallies()

    # Collect interior cells by type.
    # UniverseFactory names them "{assy_name}_interior_region", e.g. "inner_fuel_interior_region".
    inner_cells = []
    outer_cells = []
    internal_blanket_cells = []
    radial_blanket_cells = []
    for cell in geometry.get_all_cells().values():
        if not (cell.name and cell.name.endswith("_interior_region")):
            continue
        if cell.name.startswith("inner_fuel"):
            inner_cells.append(cell)
        elif cell.name.startswith("outer_fuel"):
            outer_cells.append(cell)
        elif cell.name.startswith("internal_blanket"):
            internal_blanket_cells.append(cell)
        elif cell.name.startswith("radial_blanket"):
            radial_blanket_cells.append(cell)
    blanket_cells = internal_blanket_cells + radial_blanket_cells

    fuel_cells = inner_cells + outer_cells

    # (A) Fission rate split: inner vs outer driver
    if inner_cells:
        t_inner = openmc.Tally(name="fission_inner")
        t_inner.filters = [openmc.CellFilter(inner_cells)]
        t_inner.scores = ["fission"]
        tallies.append(t_inner)

    if outer_cells:
        t_outer = openmc.Tally(name="fission_outer")
        t_outer.filters = [openmc.CellFilter(outer_cells)]
        t_outer.scores = ["fission"]
        tallies.append(t_outer)

    # (B) Flux spectrum in driver fuel — coarse few-group fast-focused bins
    e_bins = [1e-3, 1e-1, 1.0, 1e2, 1e4, 1e5, 1e6, 2e7]
    if fuel_cells:
        t_flux = openmc.Tally(name="flux_spectrum_fuel")
        t_flux.filters = [
            openmc.CellFilter(fuel_cells),
            openmc.EnergyFilter(e_bins)
        ]
        t_flux.scores = ["flux"]
        tallies.append(t_flux)

    # (C) Blanket breeding tallies — one per blanket zone to avoid ordering ambiguity.
    #     capture = absorption - fission; capture on U-238 drives Pu-239 production.
    if internal_blanket_cells:
        t_int = openmc.Tally(name="blanket_breeding_internal")
        t_int.filters = [openmc.CellFilter(internal_blanket_cells)]
        t_int.scores = ["absorption", "fission"]
        tallies.append(t_int)

    if radial_blanket_cells:
        t_rad = openmc.Tally(name="blanket_breeding_radial")
        t_rad.filters = [openmc.CellFilter(radial_blanket_cells)]
        t_rad.scores = ["absorption", "fission"]
        tallies.append(t_rad)

        t_blanket_flux = openmc.Tally(name="flux_spectrum_blanket")
        t_blanket_flux.filters = [
            openmc.CellFilter(blanket_cells),
            openmc.EnergyFilter(e_bins)
        ]
        t_blanket_flux.scores = ["flux"]
        tallies.append(t_blanket_flux)

    # (C2) Conversion ratio = fissile production / fissile destruction:
    #      (U-238 + Pu-240 capture) / (U-235 + Pu-239 + Pu-241 absorption)
    #      Same definition as build_model_per_assembly; one bin per zone cell
    #      so driver and blanket contributions can be separated.
    cr_cells = fuel_cells + blanket_cells
    if cr_cells:
        t_cr_prod = openmc.Tally(name="cr_fissile_production")
        t_cr_prod.filters = [openmc.CellFilter(cr_cells)]
        t_cr_prod.nuclides = ["U238", "Pu240"]
        t_cr_prod.scores = ["(n,gamma)"]
        tallies.append(t_cr_prod)

        t_cr_dest = openmc.Tally(name="cr_fissile_destruction")
        t_cr_dest.filters = [openmc.CellFilter(cr_cells)]
        t_cr_dest.nuclides = ["U235", "Pu239", "Pu241"]
        t_cr_dest.scores = ["absorption"]
        tallies.append(t_cr_dest)

    # (D) Vessel leakage — radial + axial, total and energy-resolved
    t_rad_total = openmc.Tally(name="leakage_radial_total")
    t_rad_total.filters = [openmc.SurfaceFilter(r_ves)]
    t_rad_total.scores = ["current"]
    tallies.append(t_rad_total)

    t_rad_spectrum = openmc.Tally(name="leakage_radial_spectrum")
    t_rad_spectrum.filters = [openmc.SurfaceFilter(r_ves), openmc.EnergyFilter(e_bins)]
    t_rad_spectrum.scores = ["current"]
    tallies.append(t_rad_spectrum)

    t_ax_top = openmc.Tally(name="leakage_axial_top")
    t_ax_top.filters = [openmc.SurfaceFilter(z_top_vac)]
    t_ax_top.scores = ["current"]
    tallies.append(t_ax_top)

    t_ax_top_spectrum = openmc.Tally(name="leakage_axial_top_spectrum")
    t_ax_top_spectrum.filters = [openmc.SurfaceFilter(z_top_vac), openmc.EnergyFilter(e_bins)]
    t_ax_top_spectrum.scores = ["current"]
    tallies.append(t_ax_top_spectrum)

    t_ax_bot = openmc.Tally(name="leakage_axial_bot")
    t_ax_bot.filters = [openmc.SurfaceFilter(z_bot_vac)]
    t_ax_bot.scores = ["current"]
    tallies.append(t_ax_bot)

    t_ax_bot_spectrum = openmc.Tally(name="leakage_axial_bot_spectrum")
    t_ax_bot_spectrum.filters = [openmc.SurfaceFilter(z_bot_vac), openmc.EnergyFilter(e_bins)]
    t_ax_bot_spectrum.scores = ["current"]
    tallies.append(t_ax_bot_spectrum)

    model = openmc.Model(geometry, materials, settings, tallies=tallies)
    return model


def build_model_per_assembly(
    positions,               # list[AssemblyPosition] from core_maps.loader
    per_asm_fuel_mats: dict, # {material_id: openmc.Material} for fuel positions
    base_mats_dict: dict,    # structural materials: sodium, ht9, reflector, etc.
    settings_override: dict | None = None,
    axial: Prism3D | None = None,  # default Prism3D(); see core_maps.loader.load_axial
) -> openmc.Model:
    """
    Build an OpenMC model where each fuel assembly position has its own
    unique OpenMC Material — required for per-assembly burnup tracking,
    shuffling, and the NFRC coupling loop.

    Non-fuel positions (reflector, shield, control, etc.) continue to use
    shared materials from base_mats_dict via the standard UniverseFactory.

    Args:
        positions:          list of AssemblyPosition objects (from load_core)
        per_asm_fuel_mats:  {material_id: openmc.Material} for each fuel slot
        base_mats_dict:     structural material dict (same structure as mats_dict
                            returned by make_materials / make_materials_from_loadout)
        settings_override:  optional dict of Settings attribute overrides
                            (e.g. {"particles": 5_000_000, "batches": 600})
        axial:              axial geometry (active height)
    """
    ag = AssemblyGeom()
    axial = axial or Prism3D()
    vg = VesselGeom()

    uf = UniverseFactory(ag, base_mats_dict)
    # Build the shared non-fuel universes
    shared_univs = uf.build_all()

    # Build per-assembly universes for fuel positions; fall back to shared for others
    position_univs: dict = {}
    for p in positions:
        if p.material_id in per_asm_fuel_mats:
            position_univs[(p.q, p.r)] = uf.build_for_material(
                p.material_id, per_asm_fuel_mats[p.material_id]
            )
        else:
            position_univs[(p.q, p.r)] = shared_univs[p.universe]

    # Flat assignment for the ring builder: (q,r) -> unique universe name
    # We use f"u_pos_{material_id}" for fuel and the standard name for others
    assignment = {}
    all_univs = dict(shared_univs)  # start with shared universes
    for p in positions:
        if p.material_id in per_asm_fuel_mats:
            key = f"u_pos_{p.material_id}"
            assignment[(p.q, p.r)] = key
            all_univs[key] = position_univs[(p.q, p.r)]
        else:
            assignment[(p.q, p.r)] = p.universe

    max_ring = max(hex_distance((p.q, p.r)) for p in positions)

    ring_names, _center = build_hexlattice_rings_from_assignment(
        assignment=assignment,
        max_ring=max_ring,
        default_universe="u_reflector",
        start_angle_deg=90.0,
    )

    lattice = openmc.HexLattice(name="core_hex_lattice")
    lattice.center = (0.0, 0.0)
    lattice.pitch = (ag.assembly_pitch_cm,)
    lattice.universes = [[all_univs[name] for name in ring] for ring in ring_names]

    outer_u = openmc.Universe(name="outer_sodium")
    outer_u.add_cell(openmc.Cell(fill=base_mats_dict["sodium"]))
    lattice.outer = outer_u

    # Vessel geometry (identical to build_model)
    H = axial.active_height_cm
    r_lat  = openmc.ZCylinder(r=vg.lattice_bound_r_cm)
    r_pool = openmc.ZCylinder(r=vg.ihx_pool_outer_r_cm)
    r_ves  = openmc.ZCylinder(r=vg.vessel_outer_r_cm, boundary_type="vacuum")
    z_bot_vac  = openmc.ZPlane(z0=-(H/2 + vg.lower_plenum_h_cm + vg.end_plate_cm), boundary_type="vacuum")
    z_bot_plen = openmc.ZPlane(z0=-(H/2 + vg.lower_plenum_h_cm))
    z_act_bot  = openmc.ZPlane(z0=-H/2)
    z_act_top  = openmc.ZPlane(z0=+H/2)
    z_top_plen = openmc.ZPlane(z0=+(H/2 + vg.upper_plenum_h_cm))
    z_top_vac  = openmc.ZPlane(z0=+(H/2 + vg.upper_plenum_h_cm + vg.end_plate_cm), boundary_type="vacuum")

    root_cells = [
        openmc.Cell(name="active_core",        fill=lattice,                      region=-r_lat  & +z_act_bot & -z_act_top),
        openmc.Cell(name="ihx_pool",           fill=base_mats_dict["pool_ihx"],   region=+r_lat  & -r_pool & +z_act_bot & -z_act_top),
        openmc.Cell(name="vessel_wall_active", fill=base_mats_dict["ht9"],        region=+r_pool & -r_ves  & +z_act_bot & -z_act_top),
        openmc.Cell(name="lower_plenum",       fill=base_mats_dict["lower_plenum"], region=-r_pool & +z_bot_plen & -z_act_bot),
        openmc.Cell(name="vessel_wall_lower",  fill=base_mats_dict["ht9"],        region=+r_pool & -r_ves  & +z_bot_plen & -z_act_bot),
        openmc.Cell(name="upper_plenum",       fill=base_mats_dict["upper_plenum"], region=-r_pool & +z_act_top & -z_top_plen),
        openmc.Cell(name="vessel_wall_upper",  fill=base_mats_dict["ht9"],        region=+r_pool & -r_ves  & +z_act_top & -z_top_plen),
        openmc.Cell(name="lower_end_plate",    fill=base_mats_dict["ht9"],        region=-r_ves  & +z_bot_vac & -z_bot_plen),
        openmc.Cell(name="upper_end_plate",    fill=base_mats_dict["ht9"],        region=-r_ves  & +z_top_plen & -z_top_vac),
    ]
    geometry = openmc.Geometry(openmc.Universe(cells=root_cells))

    # Structural materials plus the per-assembly fuel.  The zone-default fuel
    # materials in base_mats_dict are depletable but unplaced (no volume), so
    # they must stay out of the list or the depletion operator rejects them.
    all_mats = [m for m in base_mats_dict.values() if not m.depletable]
    all_mats += list(per_asm_fuel_mats.values())
    materials = openmc.Materials(all_mats)

    settings = openmc.Settings()
    settings.run_mode = "eigenvalue"
    settings.batches = 600
    settings.inactive = 100
    settings.particles = 5_000_000
    settings.temperature = {"method": "nearest"}
    settings.entropy_mesh = _entropy_mesh(vg, axial)
    bounds = [-250, -250, -axial.active_height_cm / 2.0, 250, 250, axial.active_height_cm / 2.0]
    settings.source = openmc.IndependentSource(
        space=openmc.stats.Box(bounds[:3], bounds[3:]),
        constraints={"fissionable": True},
    )
    if settings_override:
        for attr, val in settings_override.items():
            setattr(settings, attr, val)
    # After the overrides, so the statepoint lands on the final batch
    settings.statepoint = {"batches": [settings.batches]}

    # Tallies — same as build_model plus CR and per-assembly fission tallies
    tallies = openmc.Tallies()
    e_bins = [1e-3, 1e-1, 1.0, 1e2, 1e4, 1e5, 1e6, 2e7]

    fuel_mats = list(per_asm_fuel_mats.values())
    if fuel_mats:
        t_flux = openmc.Tally(name="flux_spectrum_fuel")
        t_flux.filters = [openmc.MaterialFilter(fuel_mats), openmc.EnergyFilter(e_bins)]
        t_flux.scores = ["flux"]
        tallies.append(t_flux)

        # Conversion ratio = fissile production / fissile destruction:
        #   (U-238 + Pu-240 capture) / (U-235 + Pu-239 + Pu-241 absorption)
        t_cr_prod = openmc.Tally(name="cr_fissile_production")
        t_cr_prod.filters = [openmc.MaterialFilter(fuel_mats)]
        t_cr_prod.nuclides = ["U238", "Pu240"]
        t_cr_prod.scores = ["(n,gamma)"]
        tallies.append(t_cr_prod)

        t_cr_dest = openmc.Tally(name="cr_fissile_destruction")
        t_cr_dest.filters = [openmc.MaterialFilter(fuel_mats)]
        t_cr_dest.nuclides = ["U235", "Pu239", "Pu241"]
        t_cr_dest.scores = ["absorption"]
        tallies.append(t_cr_dest)

        # Per-assembly fission rate: power shape and peak burnup (T2.2)
        t_fis_asm = openmc.Tally(name="fission_by_assembly")
        t_fis_asm.filters = [openmc.MaterialFilter(fuel_mats)]
        t_fis_asm.scores = ["fission"]
        tallies.append(t_fis_asm)

    t_rad_total = openmc.Tally(name="leakage_radial_total")
    t_rad_total.filters = [openmc.SurfaceFilter(r_ves)]
    t_rad_total.scores = ["current"]
    tallies.append(t_rad_total)

    t_ax_top = openmc.Tally(name="leakage_axial_top")
    t_ax_top.filters = [openmc.SurfaceFilter(z_top_vac)]
    t_ax_top.scores = ["current"]
    tallies.append(t_ax_top)

    t_ax_bot = openmc.Tally(name="leakage_axial_bot")
    t_ax_bot.filters = [openmc.SurfaceFilter(z_bot_vac)]
    t_ax_bot.scores = ["current"]
    tallies.append(t_ax_bot)

    return openmc.Model(geometry, materials, settings, tallies=tallies)


def write_run_info(run_dir: str, loadout_path: str, profile: str, deplete: bool) -> None:
    """Record loadout + profile + tuning constants so postprocess can attach them to metadata."""
    import json
    import materials.tier1_homog as t1

    with open(loadout_path) as f:
        loadout = json.load(f)
    prof = RUN_PROFILES[profile]
    info = {
        "loadout":        loadout_path,
        "loadout_name":   loadout.get("name"),
        "profile":        profile,
        "particles":      prof["particles"],
        "batches":        prof["batches"],
        "inactive":       prof["inactive"],
        "timesteps_days": prof["timesteps"] if deplete else [],
        "power_MWth":     POWER_MWth,
        "active_height_cm": load_axial(loadout_path).active_height_cm,
        "fissile_scale":  t1.FISSILE_SCALE,
        "fuel_params":    {k: at["fuel_params"] for k, at in loadout["assembly_types"].items()
                           if "fuel_params" in at},
    }
    with open(os.path.join(run_dir, "run_info.json"), "w") as f:
        json.dump(info, f, indent=2)


if __name__ == "__main__":
    import argparse
    import sys

    ap = argparse.ArgumentParser(description="Build and run a PRISM Tier-1 core from a JSON loadout.")
    ap.add_argument("loadout", help="e.g. core_loadouts/prism_unf_recycle.json")
    ap.add_argument("--deplete", action="store_true", help="run depletion instead of a single transport solve")
    ap.add_argument("--profile", choices=RUN_PROFILES, default="standard",
                    help="transport statistics + depletion schedule (default: standard)")
    opts = ap.parse_args()

    loadout_path = opts.loadout
    deplete = opts.deplete
    profile = opts.profile
    timesteps = RUN_PROFILES[profile]["timesteps"]

    loadout_name = os.path.splitext(os.path.basename(loadout_path))[0]
    run_tag = loadout_name if profile == "standard" else f"{loadout_name}__{profile}"
    run_dir = os.path.join("runs", run_tag)
    os.makedirs(run_dir, exist_ok=True)

    print(f"\n=== Loadout: {loadout_path} [{profile}] → {run_dir} ===")
    model = build_model(loadout_path, profile=profile)
    model.export_to_xml(directory=run_dir)
    write_run_info(run_dir, loadout_path, profile, deplete)

    if deplete:
        n_steps = len(timesteps)
        dep_h5 = os.path.join(run_dir, "depletion_results.h5")
        prev_results = None
        completed_steps = 0
        if os.path.exists(dep_h5):
            prev_results = openmc.deplete.Results(dep_h5)
            completed_steps = len(prev_results) - 1  # len includes initial state
            if completed_steps >= n_steps:
                print(f"  Already complete ({completed_steps}/{n_steps} steps). Nothing to do.")
                sys.exit(0)
            print(f"  Resuming from step {completed_steps}/{n_steps} "
                  f"({sum(timesteps[:completed_steps])} → {sum(timesteps)} d)")
        else:
            print(f"  Running depletion: {n_steps} steps / {sum(timesteps)} d @ {POWER_MWth} MWth "
                  f"(CF={CAPACITY_FACTOR})")
        remaining = timesteps[completed_steps:]
        with _working_dir(run_dir):
            operator = openmc.deplete.CoupledOperator(model, CHAIN_FILE, prev_results=prev_results)
            integrator = openmc.deplete.CECMIntegrator(
                operator,
                timesteps=remaining,
                power=POWER_W,
                timestep_units="d",
            )
            integrator.integrate()
        print(f"  Depletion complete → {run_dir}/depletion_results.h5")
    else:
        sp = model.run(cwd=run_dir)
        print(f"  statepoint: {sp}")