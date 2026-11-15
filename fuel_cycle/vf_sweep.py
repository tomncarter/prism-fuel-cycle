# fuel_cycle/vf_sweep.py
"""
Fuel-volume-fraction sweep: BOC loading search for each core at each
Till & Chang Table 14-3 pin design, then TRU and conversion ratio at the
searched point.  Calibration aid: which pin design reproduces each Triplett
core's TRU content, and does its CR follow?

    # one core, all designs (run several cores as parallel processes)
    python -m fuel_cycle.vf_sweep run core_loadouts/prism_breeder.json
    # collect everything under runs/vf_sweep/ into a table and chart
    python -m fuel_cycle.vf_sweep summarize

Only the volume fractions change; the assembly pitch, duct and layout stay at
the loadout's values (the book's pitch grows with pin diameter, ours does not).
"""
from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path

ROOT = Path("runs/vf_sweep")
DESIGNS = ["tc_d090", "tc_d100", "tc_d110", "tc_d120"]

# Triplett 2012 Table II (driver TRU wt%, CR); UNF recycle TRU is inner/outer
PAPER = {
    "prism_unf_recycle": dict(tru=(18.9, 22.7), cr=0.72),
    "prism_breakeven":   dict(tru=(21.1,), cr=1.06),
    "prism_breeder":     dict(tru=(20.9,), cr=1.23),
    "prism_weapons_pu":  dict(tru=(18.3,), cr=1.00),
}


def run(loadout: str, designs: list[str], k_target: float, profile: str) -> None:
    from fuel_cycle.loading_search import search_boc_loading

    src = Path(loadout)
    doc = json.loads(src.read_text())
    for design in designs:
        d = ROOT / f"{src.stem}__{design}"
        d.mkdir(parents=True, exist_ok=True)
        if (d / "search" / "search_result.json").exists():
            print(f"  {d}: done, skipping")
            continue
        trial = copy.deepcopy(doc)
        trial["pin_design"] = design
        for at in trial["assembly_types"].values():       # top-level design applies to all
            fp = at.get("fuel_params")
            if fp:
                fp.pop("pin_design", None)
                fp.pop("volume_fractions", None)
        (d / "loadout.json").write_text(json.dumps(trial, indent=2))
        search_boc_loading(str(d / "loadout.json"), k_target=k_target, profile=profile,
                           run_dir=str(d / "search"))


def _cr_at(statepoint: Path) -> tuple[float, float, float]:
    """(CR, driver part, blanket part) from the cr_fissile_* tallies."""
    import openmc
    with openmc.StatePoint(str(statepoint)) as sp:
        tp = sp.get_tally(name="cr_fissile_production")
        td = sp.get_tally(name="cr_fissile_destruction")
        cells = list(tp.filters[0].bins)
        names = {c: sp.summary.geometry.get_all_cells()[c].name for c in cells}
    p = tp.mean.reshape(len(cells), -1).sum(1)
    dest = td.mean.reshape(len(cells), -1).sum(1).sum()
    blanket = sum(p[i] for i, c in enumerate(cells) if "blanket" in names[c])
    return p.sum() / dest, (p.sum() - blanket) / dest, blanket / dest


def summarize() -> "pd.DataFrame":
    import pandas as pd
    from materials.tier1_homog import PIN_DESIGNS

    rows = []
    for d in sorted(ROOT.glob("*__tc_d*")):
        res_file = d / "search" / "search_result.json"
        if not res_file.exists():
            continue
        core, design = d.name.rsplit("__", 1)
        res = json.loads(res_file.read_text())
        r = res["result"]
        i = min(range(len(r["history"])), key=lambda j: abs(r["history"][j]["x"] - r["x"]))
        sp = sorted((d / "search" / f"iter_{i:02d}").glob("statepoint.*.h5"))[-1]
        cr, cr_drv, cr_bl = _cr_at(sp)
        at = res["loadout"]["assembly_types"]
        tru = [100 * at[t]["fuel_params"]["pu_hm_frac"] for t in ("inner_fuel", "outer_fuel") if t in at]
        pd_ = PIN_DESIGNS[design]
        rows.append(dict(
            core=core, pin_design=design, pin_od_cm=pd_["pin_od_cm"], fuel_vf_cell=pd_["fuel"],
            tru_inner_pct=tru[0], tru_outer_pct=tru[-1], keff=r["k"], keff_unc=r["k_unc"],
            converged=r["converged"], cr=cr, cr_driver=cr_drv, cr_blanket=cr_bl,
            paper_tru_inner_pct=PAPER.get(core, {}).get("tru", (None,))[0],
            paper_tru_outer_pct=PAPER.get(core, {}).get("tru", (None,))[-1],
            paper_cr=PAPER.get(core, {}).get("cr"),
        ))
    df = pd.DataFrame(rows).sort_values(["core", "fuel_vf_cell"]).reset_index(drop=True)
    df.to_csv(ROOT / "summary.csv", index=False)
    print(df.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    _plot(df, ROOT / "vf_sweep.png")
    print(f"  wrote {ROOT / 'summary.csv'} and {ROOT / 'vf_sweep.png'}")
    return df


def _plot(df, out: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # categorical order (fixed): one hue per core
    colors = {"prism_unf_recycle": "#2a78d6", "prism_breakeven": "#eb6834",
              "prism_breeder": "#1baf7a", "prism_weapons_pu": "#eda100"}
    labels = {"prism_unf_recycle": "UNF recycle", "prism_breakeven": "Breakeven",
              "prism_breeder": "Breeder", "prism_weapons_pu": "Weapons Pu"}
    ink, muted, grid = "#1f1f1f", "#6b6b6b", "#e4e4e4"

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), dpi=150)
    for ax, col, paper_col, ylab in (
        (axes[0], "tru_outer_pct", "paper_tru_outer_pct", "Driver TRU for k = 1.03 (wt% HM)"),
        (axes[1], "cr", "paper_cr", "BOC conversion ratio"),
    ):
        for core, g in df.groupby("core"):
            c = colors.get(core, "#888888")
            ax.plot(g.fuel_vf_cell, g[col], color=c, lw=2, marker="o", ms=6,
                    label=labels.get(core, core))
            if g[paper_col].notna().any():
                ax.axhline(g[paper_col].iloc[0], color=c, lw=1, ls=(0, (4, 3)))
        ax.set_xlabel("Fuel volume fraction of the lattice cell (Table 14-3)", color=ink)
        ax.set_ylabel(ylab, color=ink)
        ax.grid(True, color=grid, lw=0.8)
        ax.set_axisbelow(True)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        ax.tick_params(colors=muted)
        top = ax.secondary_xaxis("top", functions=(lambda v: v, lambda v: v))
        top.set_xticks(sorted(df.fuel_vf_cell.unique()))
        top.set_xticklabels([f"{od:.1f} cm" for od in sorted(df.pin_od_cm.unique())], color=muted)
        top.set_xlabel("Pin OD (169 pins)", color=muted)
    axes[1].legend(frameon=False, loc="best")
    fig.suptitle("Fuel volume fraction sweep: solid = model at k = 1.03 (BOC, quick), "
                 "dashed = Triplett 2012 Table II.  Development check, not a result.",
                 fontsize=10, color=ink)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("loadout")
    r.add_argument("--designs", nargs="+", default=DESIGNS)
    r.add_argument("--k-target", type=float, default=1.03)
    r.add_argument("--profile", default="quick")
    sub.add_parser("summarize")
    opts = ap.parse_args()
    if opts.cmd == "run":
        run(opts.loadout, opts.designs, opts.k_target, opts.profile)
    else:
        summarize()


if __name__ == "__main__":
    main()
