#!/usr/bin/env python3
"""
Geometry figures for a core loadout: 2-D slices and 3-D ray-traced views.

    ./run_openmc.sh python visualize_geometry.py                      # breeder core (default)
    ./run_openmc.sh python visualize_geometry.py core_loadouts/prism_unf_recycle.json --pixels 2400

Writes PNGs to figures/geometry/<loadout>/ :
    core_map               radial midplane slice with a legend and assembly counts
    xy_midplane            the same slice without annotation
    xy_assembly_zoom       centre assembly and neighbours (duct, inter-duct gap, interior)
    yz_axial               axial slice through the core centre (x = 0)
    raytrace_solid_cutaway 3-D view with sodium, vessel and plena hidden
    raytrace_wireframe     3-D view: faint fuel/blankets, outlined absorber and GEM columns

Colours follow the cell (assembly type), not the material: withdrawn control
positions, GEMs and the source are all sodium-filled ducts and would otherwise
look identical.  Plotting runs OpenMC in plot mode; no transport is done.
"""
import argparse
import os
import shutil
import tempfile

import openmc

from build_prism_core import build_model
from core_maps.loader import load_core, type_counts
from geometry.params import AssemblyGeom, VesselGeom

DEFAULT_LOADOUT = "core_loadouts/prism_breeder.json"

# (cell-name prefix, legend label, RGB).  First match wins, so order matters.
ASSEMBLY_STYLE = [
    ("inner_fuel_interior",        "Inner driver",              (214, 69, 65)),
    ("outer_fuel_interior",        "Outer driver",              (240, 140, 50)),
    ("internal_blanket_interior",  "Internal blanket",          (90, 140, 200)),
    ("radial_blanket_interior",    "Radial blanket",            (45, 75, 150)),
    ("control_rod_interior",       "Control (B4C inserted)",    (70, 70, 70)),
    ("ult_shutdown_interior",      "Ultimate shutdown (B4C)",   (20, 20, 20)),
    ("control_withdrawn_interior", "Control / USS (withdrawn)", (150, 110, 190)),
    ("gem_interior",               "Gas expansion module",      (60, 175, 165)),
    ("source_interior",            "Neutron source",            (205, 55, 150)),
    ("reflector_interior",         "Reflector",                 (150, 120, 90)),
    ("shield_interior",            "Shield",                    (110, 150, 100)),
]
STRUCTURE_STYLE = [
    ("_duct",       (170, 170, 175)),
    ("vessel_wall", (120, 120, 128)),
    ("end_plate",   (120, 120, 128)),
    ("ihx_pool",    (190, 210, 230)),
    ("plenum",      (215, 225, 235)),
]
SODIUM = (250, 230, 160)
POOL = (190, 210, 230)

# universe name -> cell prefix of its interior, for legend counts
UNIVERSE_PREFIX = {
    "u_inner_fuel": "inner_fuel_interior", "u_outer_fuel": "outer_fuel_interior",
    "u_internal_blanket": "internal_blanket_interior", "u_radial_blanket": "radial_blanket_interior",
    "u_control": "control_rod_interior", "u_control_rod": "control_rod_interior",
    "u_ult_shutdown": "ult_shutdown_interior", "u_ultimate_shutdown": "ult_shutdown_interior",
    "u_control_withdrawn": "control_withdrawn_interior", "u_gem": "gem_interior",
    "u_source": "source_interior", "u_reflector": "reflector_interior", "u_shield": "shield_interior",
}


def cell_colors(model):
    colors = {}
    for cell in model.geometry.get_all_cells().values():
        name = cell.name or ""
        rgb = next((c for p, _, c in ASSEMBLY_STYLE if name.startswith(p)), None)
        if rgb is None:
            rgb = next((c for p, c in STRUCTURE_STYLE if p in name), None)
        if rgb is None and cell.fill_type == "material":
            rgb = SODIUM          # assembly background coolant, empty slots, outer pool
        if rgb is not None:
            colors[cell] = rgb
    return colors


def cells_named(model, *prefixes):
    return [c for c in model.geometry.get_all_cells().values()
            if any((c.name or "").startswith(p) for p in prefixes)]


def make_plots(model, pixels):
    ag, vg = AssemblyGeom(), VesselGeom()
    colors = cell_colors(model)
    r_out = vg.ihx_pool_outer_r_cm + vg.vessel_wall_cm
    plots = []

    def slice_plot(name, basis, width):
        p = openmc.Plot(name=name)
        p.filename = name
        p.basis = basis
        p.origin = (0, 0, 0)
        p.width = width
        p.pixels = (pixels, max(1, int(round(pixels * width[1] / width[0]))))
        p.color_by = "cell"
        p.colors = colors
        p.background = (255, 255, 255)
        return p

    full = 2 * r_out + 10
    plots.append(slice_plot("xy_midplane", "xy", (full, full)))
    zw = 3.4 * ag.assembly_pitch_cm
    plots.append(slice_plot("xy_assembly_zoom", "xy", (zw, zw)))
    # Cut along x = 0: through assembly centres and across flats.  A y = 0 cut runs
    # along shared duct edges and shows the inter-duct sodium as wide stripes.
    zs = [s.z0 for s in model.geometry.get_all_surfaces().values() if isinstance(s, openmc.ZPlane)]
    plots.append(slice_plot("yz_axial", "yz", (full, max(zs) - min(zs) + 10)))

    cam = (2.2 * r_out, -1.6 * r_out, 1.4 * r_out)
    raster = (pixels, int(pixels * 0.75))

    def camera(p):
        p.camera_position = cam
        p.look_at = (0, 0, 0)
        p.up = (0, 0, 1)
        p.horizontal_field_of_view = 40
        p.pixels = raster
        p.color_by = "cell"
        p.colors = colors
        p.background = (255, 255, 255)
        return p

    solid = camera(openmc.SolidRayTracePlot(name="raytrace_solid_cutaway"))
    solid.filename = "raytrace_solid_cutaway"
    solid.opaque_domains = cells_named(model, *(p for p, _, _ in ASSEMBLY_STYLE))
    solid.light_position = (r_out, -2 * r_out, 3 * r_out)
    solid.diffuse_fraction = 0.25
    plots.append(solid)

    wire = camera(openmc.WireframeRayTracePlot(name="raytrace_wireframe"))
    wire.filename = "raytrace_wireframe"
    wire.set_transparent(model.geometry)
    for c in cells_named(model, "inner_fuel_interior", "outer_fuel_interior",
                         "internal_blanket_interior", "radial_blanket_interior"):
        wire.xs[c] = 0.004
    highlighted = cells_named(model, "control_", "ult_shutdown_interior",
                              "gem_interior", "source_interior")
    for c in highlighted:
        wire.xs[c] = 0.05
    wire.wireframe_domains = highlighted
    wire.wireframe_color = (0, 0, 0)
    wire.wireframe_thickness = 1
    plots.append(wire)
    return plots


def annotate_core_map(png_in, png_out, loadout, positions):
    """Midplane slice plus a legend of the assembly types present, with counts."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    n_by_prefix = {}
    for p in positions:
        pre = UNIVERSE_PREFIX.get(p.universe)
        if pre:
            n_by_prefix[pre] = n_by_prefix.get(pre, 0) + 1

    def patch(rgb, label):
        return Patch(facecolor=[v / 255 for v in rgb], edgecolor="0.3", label=label)

    handles = [patch(rgb, f"{label}  ({n_by_prefix[pre]})")
               for pre, label, rgb in ASSEMBLY_STYLE if pre in n_by_prefix]
    handles += [patch(SODIUM, "Sodium"), patch(POOL, "IHX / pool region")]

    n_assy = sum(v for k, v in type_counts(positions).items() if k != "empty")
    fig, ax = plt.subplots(figsize=(13, 9), dpi=150)
    ax.imshow(plt.imread(png_in))
    ax.axis("off")
    ax.set_title(f"{os.path.splitext(os.path.basename(loadout))[0]}: "
                 f"{n_assy} assemblies, core midplane", fontsize=13)
    ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.0, 0.5),
              frameon=False, fontsize=11)
    fig.tight_layout()
    fig.savefig(png_out, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("loadout", nargs="?", default=DEFAULT_LOADOUT,
                    help=f"core loadout JSON (default {DEFAULT_LOADOUT})")
    ap.add_argument("--pixels", type=int, default=1600, help="horizontal resolution")
    ap.add_argument("--out", default=None, help="output dir (default figures/geometry/<loadout>)")
    opts = ap.parse_args()

    name = os.path.splitext(os.path.basename(opts.loadout))[0]
    out_dir = opts.out or os.path.join("figures", "geometry", name)
    os.makedirs(out_dir, exist_ok=True)

    model = build_model(opts.loadout, profile="smoke")
    model.plots = openmc.Plots(make_plots(model, opts.pixels))

    with tempfile.TemporaryDirectory() as work:
        model.plot_geometry(cwd=work)
        for f in sorted(os.listdir(work)):
            if f.endswith(".png"):
                shutil.copy(os.path.join(work, f), os.path.join(out_dir, f))
                print(f"  wrote {os.path.join(out_dir, f)}")

    positions, _ = load_core(opts.loadout)
    out = os.path.join(out_dir, "core_map.png")
    annotate_core_map(os.path.join(out_dir, "xy_midplane.png"), out, opts.loadout, positions)
    print(f"  wrote {out}")


if __name__ == "__main__":
    main()
