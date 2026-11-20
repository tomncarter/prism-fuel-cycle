"""Generate a Q,R coordinate reference map for the PRISM core.

Hexes are coloured by zone type so the core layout is immediately readable.
Labels show (q,r) in each hex.

No OpenMC required — pure matplotlib + geometry modules.

Usage:
    python _gen_qr_map.py
    python _gen_qr_map.py --out figures/core_qr_labels.png --max-ring 11
"""
import argparse
import math
import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from geometry.hexmap import full_hex_coords, hex_distance, ring_angle_sort
from geometry.params import AssemblyGeom
from core_maps.abr_like_assignment import abr_like_assignment
from core_maps.prism_cr072_assignment import prism_cr072_assignment

ASSIGNMENTS = {
    "abr":       abr_like_assignment,
    "prism_cr072": prism_cr072_assignment,
}


def flat_axial_to_xy(q: int, r: int, size: float = 1.0):
    """Flat-top hex layout: flat edges parallel to X axis, columns staggered."""
    x = size * 1.5 * q
    y = size * (math.sqrt(3) * r + math.sqrt(3) / 2.0 * q)
    return x, y

ZONE_STYLE = {
    "u_inner_fuel":       {"color": "#4a90d9", "label": "Inner fuel (rings 0–3)"},
    "u_outer_fuel":       {"color": "#7bc67e", "label": "Outer fuel (rings 5–7)"},
    "u_internal_blanket": {"color": "#f5c842", "label": "Internal blanket (ring 4)"},
    "u_radial_blanket":   {"color": "#f4a261", "label": "Radial blanket (ring 8)"},
    "u_reflector":        {"color": "#c0c0c0", "label": "Reflector (rings 9–10)"},
    "u_shield":           {"color": "#707070", "label": "Shield (ring 11)"},
    "u_control_rod":      {"color": "#d62728", "label": "Control rod"},
    "u_ultimate_shutdown": {"color": "#9467bd", "label": "Ultimate shutdown"},
    "u_gem":              {"color": "#17becf", "label": "GEM"},
    "u_empty":            {"color": None,      "label": "Empty"},
}
DEFAULT_STYLE = {"color": "#eeeeee", "label": "Unknown"}


def generate(out_path: str = "figures/core_qr_labels.png", max_ring: int = 11,
             design: str = "abr") -> None:
    ag = AssemblyGeom()
    pitch = ag.assembly_pitch_cm
    assignment = ASSIGNMENTS.get(design, abr_like_assignment)()

    coords = full_hex_coords(max_ring)
    by_ring: dict = {}
    for c in coords:
        by_ring.setdefault(hex_distance(c), []).append(c)

    all_centers = []
    for ring in range(max_ring, 0, -1):
        for q, r in ring_angle_sort(by_ring[ring], start_angle_deg=90.0):
            x, y = flat_axial_to_xy(q, r, size=pitch)
            all_centers.append((x, y, q, r))
    x0, y0 = flat_axial_to_xy(0, 0, size=pitch)
    all_centers.append((x0, y0, 0, 0))

    # axial_to_xy places adjacent hex centers sqrt(3)*pitch apart, so circumradius = pitch
    R = pitch

    n_across = 2 * max_ring + 1          # number of hexes across the widest row
    px_per_hex = 160                      # target pixels per hex
    total_px = n_across * px_per_hex
    dpi = 150
    figsize = total_px / dpi             # inches

    fontsize = max(4, px_per_hex * 0.18)

    fig = plt.figure(figsize=(figsize, figsize * 0.92), facecolor="#e8e8e8")
    ax = fig.add_subplot(111)
    ax.set_facecolor("#e8e8e8")
    ax.set_aspect("equal")
    ax.axis("off")

    for x, y, q, r in all_centers:
        is_origin = (q == 0 and r == 0)
        if is_origin:
            facecolor = "#e040fb"
        else:
            zone = assignment.get((q, r), "u_reflector")
            style = ZONE_STYLE.get(zone, DEFAULT_STYLE)
            if style["color"] is None:
                continue  # empty position — don't draw
            facecolor = style["color"]
        ax.add_patch(RegularPolygon(
            (x, y), numVertices=6, radius=R,
            orientation=math.radians(30),
            facecolor=facecolor,
            edgecolor="#e8e8e8", linewidth=1.0,
        ))
        label = f"{q},{r}"
        fs = fontsize if len(label) <= 4 else fontsize * 0.72
        ax.text(x, y, label, ha="center", va="center",
                fontsize=fs, color="white" if is_origin else "black",
                fontweight="bold", fontfamily="monospace")

    xs = [c[0] for c in all_centers]
    ys = [c[1] for c in all_centers]
    margin = R * 1.5
    ax.set_xlim(min(xs) - margin, max(xs) + margin)
    ax.set_ylim(min(ys) - margin, max(ys) + margin)

    ax.set_title("PRISM core — (q, r) reference map", fontsize=fontsize * 2.2,
                 pad=10, fontweight="bold")

    pathlib.Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout(pad=0.5)
    plt.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    w = int(figsize * dpi)
    print(f"Written: {out_path}  (~{w}px wide)")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default="figures/core_qr_labels.png")
    parser.add_argument("--max-ring", type=int, default=11)
    parser.add_argument("--design", default="abr", choices=list(ASSIGNMENTS))
    args = parser.parse_args()
    generate(out_path=args.out, max_ring=args.max_ring, design=args.design)
