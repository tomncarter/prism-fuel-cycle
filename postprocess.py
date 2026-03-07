"""
postprocess.py — Extract OpenMC depletion results to Parquet files.

Usage:
    python postprocess.py [runs_dir] [results_dir]

Defaults:
    runs_dir    = ./runs
    results_dir = ./results

Outputs (all in results_dir):
    keff.parquet        — k-eff vs time for all permutations
    tallies.parquet     — final-step tally snapshots (fission split, breeding, leakage)
    spectra.parquet     — energy-resolved tally snapshots at final step
    inventories.parquet — key nuclide masses vs time for all depletable materials
"""
from __future__ import annotations

import sys
import datetime
from pathlib import Path
from xml.etree import ElementTree

import numpy as np
import pandas as pd
import openmc
import openmc.deplete

# Nuclides to extract from depletion results.
# The fast chain may not carry every nuclide; missing ones are silently skipped.
KEY_NUCLIDES = [
    "U235", "U238",
    "Pu238", "Pu239", "Pu240", "Pu241", "Pu242",
    "Am241", "Am243",
    "Np237",
    "Cs137", "Sr90",   # major fission products (may not be in fast chain)
]

# Tally names that are scalar (no energy bins) — go into tallies.parquet
SCALAR_TALLIES = [
    "fission_inner", "fission_outer",
    "blanket_breeding_internal", "blanket_breeding_radial",
    "leakage_radial_total", "leakage_axial_top", "leakage_axial_bot",
]

# Tally names that are energy-binned — go into spectra.parquet
SPECTRUM_TALLIES = [
    "flux_spectrum_fuel", "flux_spectrum_blanket",
    "leakage_radial_spectrum",
    "leakage_axial_top_spectrum", "leakage_axial_bot_spectrum",
]


def _find_final_statepoint(run_dir: Path) -> Path | None:
    """Return the highest-numbered statepoint.*.h5 in run_dir."""
    sps = sorted(run_dir.glob("statepoint.*.h5"),
                 key=lambda p: int(p.stem.split(".")[1]))
    return sps[-1] if sps else None


def extract_keff(results: openmc.deplete.Results, run_tag: str) -> pd.DataFrame:
    times, eigenvalues = results.get_keff(time_units="d")
    return pd.DataFrame({
        "run_tag":   run_tag,
        "time_days": times,
        "keff":      eigenvalues[:, 0],
        "keff_unc":  eigenvalues[:, 1],
    })


def extract_inventories(
    results: openmc.deplete.Results,
    run_tag: str,
    mat_keys: dict[str, str],   # {mats_dict_key: display_name}
) -> pd.DataFrame:
    """
    Extract key nuclide masses (g) and atoms for each depletable material over time.
    mat_keys maps the material name string to a human-readable label.
    """
    rows = []
    for mat_name, mat_label in mat_keys.items():
        for nuc in KEY_NUCLIDES:
            try:
                times, masses = results.get_mass(mat_name, nuc,
                                                 mass_units="g", time_units="d")
                _, atoms = results.get_atoms(mat_name, nuc,
                                             nuc_units="atoms", time_units="d")
                for t, m, a in zip(times, masses, atoms):
                    rows.append({
                        "run_tag":   run_tag,
                        "time_days": t,
                        "material":  mat_label,
                        "nuclide":   nuc,
                        "mass_g":    m,
                        "atoms":     a,
                    })
            except (KeyError, ValueError):
                pass  # nuclide not in chain or material not found
    return pd.DataFrame(rows)


def extract_tallies(sp: openmc.StatePoint, run_tag: str, time_days: float) -> pd.DataFrame:
    rows = []
    for tname in SCALAR_TALLIES:
        try:
            t = sp.get_tally(name=tname)
            df = t.get_pandas_dataframe()
            df["run_tag"]   = run_tag
            df["time_days"] = time_days
            df["tally"]     = tname
            rows.append(df)
        except Exception:
            pass
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def extract_spectra(sp: openmc.StatePoint, run_tag: str, time_days: float) -> pd.DataFrame:
    rows = []
    for tname in SPECTRUM_TALLIES:
        try:
            t = sp.get_tally(name=tname)
            df = t.get_pandas_dataframe()
            df["run_tag"]   = run_tag
            df["time_days"] = time_days
            df["tally"]     = tname
            rows.append(df)
        except Exception:
            pass
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def extract_metadata(
    run_dir: Path,
    run_tag: str,
    sp: openmc.StatePoint | None,
    results: openmc.deplete.Results | None,
) -> dict:
    """Collect run provenance into a flat dict (one row in metadata.parquet)."""
    meta: dict = {
        "run_tag":           run_tag,
        "postprocessed_at":  datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
    }

    # Parse internal/radial feed from run_tag (e.g. "NUint_DUrad")
    parts = run_tag.split("_")
    meta["internal_feed"] = parts[0].replace("int", "") if len(parts) >= 1 else None
    meta["radial_feed"]   = parts[1].replace("rad", "") if len(parts) >= 2 else None

    # Statepoint modification time — proxy for when the transport run finished
    sp_file = _find_final_statepoint(run_dir)
    if sp_file and sp_file.exists():
        mtime = sp_file.stat().st_mtime
        meta["statepoint_file"]  = sp_file.name
        meta["statepoint_mtime"] = datetime.datetime.fromtimestamp(mtime, datetime.timezone.utc).isoformat(timespec="seconds")
    else:
        meta["statepoint_file"]  = None
        meta["statepoint_mtime"] = None

    # Settings from settings.xml
    settings_xml = run_dir / "settings.xml"
    if settings_xml.exists():
        root = ElementTree.parse(settings_xml).getroot()
        meta["n_particles"] = int(root.findtext("particles") or 0)
        meta["n_batches"]   = int(root.findtext("batches")   or 0)
        meta["n_inactive"]  = int(root.findtext("inactive")  or 0)
    else:
        meta["n_particles"] = meta["n_batches"] = meta["n_inactive"] = None

    # OpenMC version from statepoint
    if sp is not None:
        meta["openmc_version"] = ".".join(str(v) for v in sp.version)
    else:
        meta["openmc_version"] = None

    # Depletion info
    if results is not None:
        times, _ = results.get_keff(time_units="d")
        meta["depleted"]          = True
        meta["n_depletion_steps"] = len(times) - 1   # steps = intervals
        meta["final_time_days"]   = float(times[-1])
        dep_file = run_dir / "depletion_results.h5"
        dep_mtime = dep_file.stat().st_mtime
        meta["depletion_mtime"] = datetime.datetime.fromtimestamp(dep_mtime, datetime.timezone.utc).isoformat(timespec="seconds")
    else:
        meta["depleted"]          = False
        meta["n_depletion_steps"] = 0
        meta["final_time_days"]   = 0.0
        meta["depletion_mtime"]   = None

    return meta


def process_run(run_dir: Path, run_tag: str):
    print(f"  Processing {run_tag} ...")

    dep_file = run_dir / "depletion_results.h5"
    sp_file  = _find_final_statepoint(run_dir)

    keff_df = inv_df = tally_df = spec_df = None
    results = None

    if dep_file.exists():
        results = openmc.deplete.Results(str(dep_file))

        keff_df = extract_keff(results, run_tag)

        # Discover depletable materials from materials.xml.
        # Results.get_mass/get_atoms key by str(material_id), not name.
        dep_mat_names = {}
        mats_xml = run_dir / "materials.xml"
        if mats_xml.exists():
            for mat in openmc.Materials.from_xml(str(mats_xml)):
                if mat.depletable:
                    dep_mat_names[str(mat.id)] = mat.name

        if dep_mat_names:
            inv_df = extract_inventories(results, run_tag, dep_mat_names)

    else:
        print(f"    WARNING: no depletion_results.h5 found in {run_dir}")

    sp_obj = None
    if sp_file and sp_file.exists():
        sp_obj = openmc.StatePoint(str(sp_file))
        final_time = float(keff_df["time_days"].iloc[-1]) if keff_df is not None and len(keff_df) else 0.0
        tally_df = extract_tallies(sp_obj, run_tag, final_time)
        spec_df  = extract_spectra(sp_obj, run_tag, final_time)
    else:
        print(f"    WARNING: no statepoint found in {run_dir}")

    meta = extract_metadata(run_dir, run_tag, sp_obj, results)
    if sp_obj is not None:
        sp_obj.close()

    return keff_df, inv_df, tally_df, spec_df, meta


def main(runs_dir: Path, results_dir: Path):
    results_dir.mkdir(parents=True, exist_ok=True)

    all_keff     = []
    all_inv      = []
    all_tallies  = []
    all_spectra  = []
    all_metadata = []

    run_dirs = sorted(runs_dir.iterdir()) if runs_dir.exists() else []
    if not run_dirs:
        print(f"No run directories found in {runs_dir}")
        return

    for run_dir in run_dirs:
        if not run_dir.is_dir():
            continue
        run_tag = run_dir.name
        k, inv, t, s, meta = process_run(run_dir, run_tag)
        if k    is not None: all_keff.append(k)
        if inv  is not None: all_inv.append(inv)
        if t    is not None and not t.empty: all_tallies.append(t)
        if s    is not None and not s.empty: all_spectra.append(s)
        all_metadata.append(meta)

    def _save(frames, name):
        if frames:
            df = pd.concat(frames, ignore_index=True)
            path = results_dir / name
            df.to_parquet(path, index=False)
            print(f"  Wrote {path}  ({len(df):,} rows)")
        else:
            print(f"  No data for {name}")

    print("\nWriting Parquet files ...")
    _save(all_keff,    "keff.parquet")
    _save(all_inv,     "inventories.parquet")
    _save(all_tallies, "tallies.parquet")
    _save(all_spectra, "spectra.parquet")

    meta_df = pd.DataFrame(all_metadata)
    meta_path = results_dir / "metadata.parquet"
    meta_df.to_parquet(meta_path, index=False)
    print(f"  Wrote {meta_path}  ({len(meta_df):,} rows)")
    print("Done.")


if __name__ == "__main__":
    _runs    = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("runs")
    _results = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("results")
    main(_runs, _results)
