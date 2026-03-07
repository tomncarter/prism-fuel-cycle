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
from materials.tier1_homog import make_materials

from geometry.hexmap import build_hexlattice_rings_from_assignment, full_hex_coords, hex_distance, axial_to_xy
from core_maps.abr_like_assignment import abr_like_assignment

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

    fig = plt.figure(figsize=(10, 10))
    ax = plt.gca()
    ax.set_aspect("equal")
    ax.axis("off")

    # pointy-top hex; RegularPolygon radius ~ pitch/sqrt(3) looks right
    hex_radius = pitch_cm / math.sqrt(3.0) * 0.98

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

    ax.set_xlim(min(xs) - pitch_cm * 2, max(xs) + pitch_cm * 2)
    ax.set_ylim(min(ys) - pitch_cm * 2, max(ys) + pitch_cm * 2)
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

    fig = plt.figure(figsize=(14, 14))
    ax = plt.gca()
    ax.set_aspect("equal")
    ax.axis("off")

    hex_radius = pitch_cm / math.sqrt(3.0) * 0.98

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
                linewidth=0.4,
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
                fontsize=6,
            )

    ax.set_xlim(min(xs) - pitch_cm * 2, max(xs) + pitch_cm * 2)
    ax.set_ylim(min(ys) - pitch_cm * 2, max(ys) + pitch_cm * 2)
    ax.set_title(f"Core reference map ({mode})")
    plt.tight_layout()
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

    fig = plt.figure(figsize=(12, 12))
    ax = plt.gca()
    ax.set_aspect("equal")
    ax.axis("off")

    hex_radius = pitch_cm / math.sqrt(3.0) * 0.98

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

    ax.set_xlim(min(xs) - pitch_cm * 2, max(xs) + pitch_cm * 2)
    ax.set_ylim(min(ys) - pitch_cm * 2, max(ys) + pitch_cm * 2)
    ax.legend(handles=legend_handles, loc="upper right", frameon=True)
    ax.set_title("Core loading map (by universe type)")
    plt.tight_layout()
    plt.savefig(out_png, dpi=250)
    plt.close(fig)


# -----------------------------
# Build OpenMC model (3D)
# -----------------------------
def build_model(max_ring: int = 11, internal_feed: str = "NU", radial_feed: str = "NU"):
    # Materials
    materials, mats_dict = make_materials(internal_feed=internal_feed, radial_feed=radial_feed)

    # Geometry params
    ag = AssemblyGeom()
    axial = Prism3D()
    vg = VesselGeom()

    # Universes
    uf = UniverseFactory(ag, mats_dict)
    univs = uf.build_all()

    # Coordinate-based assignment (Tier-1 ABR-like zoning)
    assignment = abr_like_assignment()

    # Analytic material volumes (needed for depletion normalization).
    # Interior volume per assembly = hex area of duct inner bore × active height.
    # Tier-1 note: all assemblies of the same type share one material, so volume
    # is the sum across all assemblies of that type.  For Tier-2 (pin-level),
    # switch to model.calculate_volumes() stochastic estimate.
    _assy_counts = Counter(assignment.values())
    _f2f  = ag.duct_inner_flat_to_flat_cm
    _area = (math.sqrt(3) / 2.0) * _f2f ** 2   # regular hexagon, flat-to-flat F: A = (√3/2)F²
    _vol_per_assy = _area * axial.active_height_cm
    mats_dict["fuel_driver_inner_homog"].volume  = _assy_counts["u_inner_fuel"]       * _vol_per_assy
    mats_dict["fuel_driver_outer_homog"].volume  = _assy_counts["u_outer_fuel"]       * _vol_per_assy
    mats_dict["fuel_blanket_internal_homog"].volume = _assy_counts["u_internal_blanket"] * _vol_per_assy
    mats_dict["fuel_blanket_radial_homog"].volume   = _assy_counts["u_radial_blanket"]   * _vol_per_assy

    # Debug artifacts
    write_reference_plot(
        assignment,
        max_ring=max_ring,
        pitch_cm=ag.assembly_pitch_cm,
        mode="qr",
        out_png="core_qr_labels.png",
    )

    write_reference_plot(
        assignment,
        max_ring=max_ring,
        pitch_cm=ag.assembly_pitch_cm,
        mode="both",
        out_png="core_ring_index_map.png",
    )
    write_loading_map_with_legend(assignment, max_ring=max_ring, pitch_cm=ag.assembly_pitch_cm, out_png="core_loading_map.png")
    write_core_assignment_csv(assignment, max_ring=max_ring, start_angle_deg=90.0, out_csv="core_assignment.csv")

    # Build complete rings (11 rings + center) from coordinates
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
    settings.batches = 80
    settings.inactive = 20
    settings.particles = 10000
    settings.temperature = {"method": "nearest"}
    settings.statepoint = {"batches": [settings.batches]}
    settings.statepoint["filename"] = "statepoint_boc.h5"

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


if __name__ == "__main__":
    import itertools
    import sys

    deplete = "--deplete" in sys.argv   # pass --deplete flag to run full depletion
                                        # default: single transport solve only

    feeds = ["NU", "DU"]
    for int_feed, rad_feed in itertools.product(feeds, feeds):
        tag = f"{int_feed}int_{rad_feed}rad"
        run_dir = os.path.join("runs", tag)
        os.makedirs(run_dir, exist_ok=True)
        print(f"\n=== Blanket permutation: internal={int_feed}, radial={rad_feed} ===")
        model = build_model(max_ring=11, internal_feed=int_feed, radial_feed=rad_feed)
        model.export_to_xml(directory=run_dir)

        if deplete:
            print(f"  Running depletion: {N_CYCLES} × {TIMESTEP_DAYS} d @ {POWER_MWth} MWth "
                  f"(CF={CAPACITY_FACTOR})")
            with _working_dir(run_dir):
                operator = openmc.deplete.CoupledOperator(model, CHAIN_FILE)
                integrator = openmc.deplete.CECMIntegrator(
                    operator,
                    timesteps=[TIMESTEP_DAYS] * N_CYCLES,
                    power=POWER_W,
                    timestep_units="d",
                )
                integrator.integrate()
            print(f"  Depletion complete → {run_dir}/depletion_results.h5")
        else:
            sp = model.run(cwd=run_dir)
            print(f"  statepoint: {sp}")