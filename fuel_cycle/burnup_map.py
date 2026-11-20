# fuel_cycle/burnup_map.py
"""
Per-assembly burnup for a cycle run, as a table and a core map.

    python -m fuel_cycle.burnup_map runs/prism_breakeven_k103__cycles__cycle            # last cycle
    python -m fuel_cycle.burnup_map runs/prism_breakeven_k103__cycles__cycle --cycle 12
    python -m fuel_cycle.burnup_map runs/<run> --all                                   # every cycle

Writes <run>/burnup/:
    burnup.parquet            one row per (cycle, assembly) for every cycle in the run
    cycle_NNN_burnup.png      core map of end-of-cycle burnup for the chosen cycle(s)

Method.  Every fuel and blanket assembly in a cycle run has its own depletable
material, named <type>_q<q>_r<r>_c<cycle loaded>; assemblies never move.  For
each cycle the core's thermal energy (sum of source_rate x step length in
depletion_results.h5) is shared among materials in proportion to their fissions,
taken as the drop in heavy-metal (Z >= 90) atoms over the cycle: a fission
removes one HM atom, capture and alpha decay do not.  Burnup is that energy over
the assembly's heavy-metal mass when it was loaded, summed over its cycles.
"""
from __future__ import annotations

import argparse
import json
import math
import re
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pandas as pd

AVOGADRO = 6.02214076e23
NAME_RE = re.compile(r"^(?P<type>.+)_q(?P<q>[+-]?\d+)_r(?P<r>[+-]?\d+)_c(?P<c>\d+)$")

# Sequential single-hue ramp (light -> dark blue) for magnitude
SEQ_RAMP = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7", "#3987e5",
            "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"]
EMPTY_FILL = "#ececec"
INK = "#1f1f1f"
MUTED = "#6b6b6b"


def _cycle_dirs(run: Path) -> list[tuple[int, Path]]:
    out = []
    for d in sorted(run.glob("cycle_*")):
        m = re.fullmatch(r"cycle_(\d+)", d.name)
        if m and (d / "depletion_results.h5").exists() and (d / "materials.xml").exists():
            out.append((int(m.group(1)), d))
    return out


def _material_names(materials_xml: Path) -> dict[str, str]:
    """material id (str) -> name, without parsing every nuclide into openmc objects."""
    names = {}
    for _, el in ET.iterparse(materials_xml, events=("start",)):
        if el.tag == "material":
            names[el.get("id")] = el.get("name", "")
    return names


def _cycle_table(cycle: int, d: Path) -> pd.DataFrame:
    import openmc.data
    import openmc.deplete

    res = openmc.deplete.Results(str(d / "depletion_results.h5"))
    first, last = res[0], res[-1]
    times = np.asarray(res.get_times(time_units="s"))
    power = np.array([s.source_rate for s in res], dtype=float)
    power = power.reshape(len(res), -1)[:, 0]
    energy_mwd = float(np.sum(power[:-1] * np.diff(times))) / 1e6 / 86400.0

    nucs = sorted(first.index_nuc, key=first.index_nuc.get)
    hm_idx, hm_amu = [], []
    for n in nucs:
        z, a, _ = openmc.data.zam(n)
        if z >= 90:
            hm_idx.append(first.index_nuc[n])
            hm_amu.append(openmc.data.atomic_mass(n))
    hm_idx = np.array(hm_idx)
    hm_amu = np.array(hm_amu)

    def atoms(step):
        a = np.asarray(step.data)
        return a[0] if a.ndim == 3 else a   # (n_mat, n_nuc)

    a0, a1 = atoms(first)[:, hm_idx], atoms(last)[:, hm_idx]
    names = _material_names(d / "materials.xml")

    rows = []
    for mat_id, i in first.index_mat.items():
        name = names.get(str(mat_id), "")
        m = NAME_RE.match(name)
        if not m:
            continue
        rows.append(dict(
            cycle=cycle, material=name, assembly_type=m["type"],
            q=int(m["q"]), r=int(m["r"]), cycle_loaded=int(m["c"]),
            hm_boc_kg=float(a0[i] @ hm_amu) / AVOGADRO / 1000.0,
            hm_eoc_kg=float(a1[i] @ hm_amu) / AVOGADRO / 1000.0,
            fissions=float(a0[i].sum() - a1[i].sum()),
        ))
    df = pd.DataFrame(rows)
    total = df["fissions"].clip(lower=0).sum()
    df["energy_mwd"] = energy_mwd * df["fissions"].clip(lower=0) / total if total > 0 else 0.0
    df["core_energy_mwd"] = energy_mwd
    df["cycle_days"] = (times[-1] - times[0]) / 86400.0
    return df


def burnup_table(run: Path) -> pd.DataFrame:
    """One row per (cycle, assembly) with cycle and cumulative end-of-cycle burnup."""
    cycles = _cycle_dirs(run)
    if not cycles:
        raise SystemExit(f"{run}: no cycle_NNN/ directories with depletion results")
    df = pd.concat([_cycle_table(c, d) for c, d in cycles], ignore_index=True)

    # HM mass at loading: BOC mass in the cycle it was loaded, else its first cycle here
    first_seen = df.sort_values("cycle").groupby("material").first()
    loaded = df[df["cycle"] == df["cycle_loaded"]].set_index("material")["hm_boc_kg"]
    hm_load = first_seen["hm_boc_kg"].copy()
    hm_load.update(loaded)
    df["hm_loaded_kg"] = df["material"].map(hm_load)
    df["partial_history"] = df["material"].map(
        first_seen["cycle"] > first_seen["cycle_loaded"])   # loaded before the first cycle on disk

    df = df.sort_values(["material", "cycle"])
    df["burnup_cycle_mwd_kg"] = df["energy_mwd"] / df["hm_loaded_kg"]
    df["burnup_mwd_kg"] = df.groupby("material")["burnup_cycle_mwd_kg"].cumsum()
    df["fima_pct"] = 100.0 * (1.0 - df["hm_eoc_kg"] / df["hm_loaded_kg"])
    df["cycles_in_core"] = df.groupby("material").cumcount() + 1

    # Discharged at the end of this cycle: absent from the next cycle, or for the last
    # cycle on disk absent from the run's current in-core list
    present = set(zip(df["cycle"], df["material"]))
    last_cycle = df["cycle"].max()
    in_core_next = None
    state = run / "cycle_state.json"
    if state.exists():
        st = json.loads(state.read_text())
        if st.get("current_cycle") == last_cycle + 1:
            in_core_next = {e["material_id"] for e in st["in_core"].values()}

    def discharged(row):
        if row.cycle < last_cycle:
            return (row.cycle + 1, row.material) not in present
        return None if in_core_next is None else row.material not in in_core_next

    df["discharged"] = [discharged(r) for r in df.itertuples()]
    df.insert(0, "run", run.name)
    return df.sort_values(["cycle", "q", "r"]).reset_index(drop=True)


def plot_cycle(df: pd.DataFrame, cycle: int, out_png: Path, labels: bool = True) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap, Normalize
    from matplotlib.patches import RegularPolygon

    from geometry.hexmap import axial_to_xy, full_hex_coords, hex_distance

    c = df[df["cycle"] == cycle]
    max_ring = max(11, int(max(hex_distance((q, r)) for q, r in zip(c.q, c.r))))
    cmap = LinearSegmentedColormap.from_list("seq_blue", SEQ_RAMP)
    vmax = float(np.ceil(df["burnup_mwd_kg"].max() / 20.0) * 20.0)   # common scale for all cycles
    norm = Normalize(0.0, vmax)
    radius = 1.0 / math.sqrt(3.0)

    fig, ax = plt.subplots(figsize=(11, 10), dpi=150)
    ax.set_aspect("equal")
    ax.axis("off")
    burn = {(q, r): row for q, r, row in zip(c.q, c.r, c.itertuples())}
    for q, r in full_hex_coords(max_ring):
        x, y = axial_to_xy(q, r, 1.0)
        row = burn.get((q, r))
        if row is None:
            ax.add_patch(RegularPolygon((x, y), 6, radius=radius * 0.97, orientation=math.pi / 6,
                                        facecolor=EMPTY_FILL, edgecolor="white", linewidth=0.6))
            continue
        val = row.burnup_mwd_kg
        ax.add_patch(RegularPolygon((x, y), 6, radius=radius * 0.97, orientation=math.pi / 6,
                                    facecolor=cmap(norm(val)), edgecolor="white", linewidth=0.6))
        if row.discharged:
            ax.add_patch(RegularPolygon((x, y), 6, radius=radius * 0.80, orientation=math.pi / 6,
                                        facecolor="none", edgecolor=INK, linewidth=1.4))
        if labels:
            ax.text(x, y, f"{val:.0f}", ha="center", va="center", fontsize=5.5,
                    color="white" if norm(val) > 0.5 else INK)
    lim = (max_ring + 1) * 0.9
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim - 0.6, lim)

    drv = c[~c["assembly_type"].str.contains("blanket")]
    dis = drv[drv["discharged"] == True]  # noqa: E712  (None = unknown)
    stats = [f"cycle {cycle} ({c['cycle_days'].iloc[0]:.0f} d), end-of-cycle burnup, MWd/kg HM"]
    if len(dis):
        stats.append(f"driver discharge: mean {dis.burnup_mwd_kg.mean():.1f}, peak {dis.burnup_mwd_kg.max():.1f} "
                     f"({len(dis)} assemblies)")
    stats.append(f"driver in-core mean {drv.burnup_mwd_kg.mean():.1f}, peak {drv.burnup_mwd_kg.max():.1f}")
    bl = c[c["assembly_type"].str.contains("blanket")]
    if len(bl):
        stats.append(f"blanket mean {bl.burnup_mwd_kg.mean():.1f}, peak {bl.burnup_mwd_kg.max():.1f}")
    ax.set_title(f"{df['run'].iloc[0]}\n" + "  |  ".join(stats[:2]) + "\n" + "  |  ".join(stats[2:]),
                 fontsize=10, color=INK)

    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    cb = fig.colorbar(sm, ax=ax, shrink=0.6, pad=0.01)
    cb.set_label("Burnup (MWd/kg HM)", color=INK)
    cb.outline.set_visible(False)
    ax.text(-lim, -lim - 0.4, "Outlined: discharged at end of cycle.  Grey: not depleted "
            "(reflector, shield, rods, GEMs, source).", fontsize=8, color=MUTED, va="center")
    if df["partial_history"].any():
        ax.text(-lim, -lim - 1.0, "Some assemblies were loaded before the first cycle on disk; "
                "their burnup is a lower bound.", fontsize=8, color=MUTED, va="center")
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description="Per-assembly burnup table and core map for a cycle run.")
    ap.add_argument("run_dir")
    ap.add_argument("--cycle", type=int, help="cycle to plot (default: last)")
    ap.add_argument("--all", action="store_true", help="plot every cycle")
    ap.add_argument("--no-labels", action="store_true", help="omit the value printed in each assembly")
    opts = ap.parse_args()

    run = Path(opts.run_dir)
    df = burnup_table(run)
    out = run / "burnup"
    out.mkdir(exist_ok=True)
    df.to_parquet(out / "burnup.parquet", index=False)
    print(f"  wrote {out / 'burnup.parquet'} ({len(df)} rows, cycles {df.cycle.min()}-{df.cycle.max()})")

    cycles = sorted(df.cycle.unique()) if opts.all else [opts.cycle if opts.cycle is not None else df.cycle.max()]
    for cyc in cycles:
        png = out / f"cycle_{cyc:03d}_burnup.png"
        plot_cycle(df, int(cyc), png, labels=not opts.no_labels)
        print(f"  wrote {png}")


if __name__ == "__main__":
    main()
