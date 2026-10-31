# fuel_cycle/loading_search.py
"""
BOC fissile-loading search for a fresh core (Phase 1 / Phase 2).

Scales a loadout's driver fissile loading (see fuel_cycle.search) until the
BOC k-eff of the all-fresh core matches the target, using transport-only
solves of the standard Tier-1 model.

Writes runs/<loadout>__search[__<profile>]/:
    iter_NN/             one transport solve per trial scale
    search_result.json   history, accepted point, and the searched loadout

Usage:
    python -m fuel_cycle.loading_search core_loadouts/prism_breakeven.json \
        [--k-target 1.005] [--profile quick] [--write-loadout PATH]
"""
from __future__ import annotations

import copy
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import openmc

from build_prism_core import build_model, _working_dir
from fuel_cycle.search import (
    apply_loading_scale, loading_knob, scale_bounds, secant_search,
)

K_TARGET = 1.005    # T1.1 nominal BOEC k-eff
K_MAX = 1.010       # T7.4 excursion limit


def search_boc_loading(
    loadout_path: str,
    k_target: float = K_TARGET,
    profile: str = "quick",
    run_dir: str | None = None,
    tol_pcm: float = 100.0,
    max_iter: int = 8,
) -> dict:
    """Run the search; return the search_result.json contents."""
    loadout_path = os.path.abspath(loadout_path)
    with open(loadout_path) as f:
        doc = json.load(f)
    name = Path(loadout_path).stem
    run = Path(run_dir or (f"runs/{name}__search" if profile == "standard"
                           else f"runs/{name}__search__{profile}"))
    run.mkdir(parents=True, exist_ok=True)

    knob = loading_knob(doc["assembly_types"])
    bounds = scale_bounds(doc["assembly_types"])
    n_iter = 0

    def evaluate(scale: float) -> tuple[float, float]:
        nonlocal n_iter
        d = (run / f"iter_{n_iter:02d}").resolve()
        n_iter += 1
        d.mkdir(exist_ok=True)
        trial = copy.deepcopy(doc)
        trial["assembly_types"] = apply_loading_scale(doc["assembly_types"], scale)
        trial_path = d / "loadout.json"
        trial_path.write_text(json.dumps(trial, indent=2))
        with _working_dir(str(d)):   # build_model drops its plots/CSV in cwd
            model = build_model(str(trial_path), profile=profile)
        with openmc.StatePoint(model.run(cwd=str(d))) as sp:
            k = sp.keff
        print(f"  [{knob}] scale {scale:.5f} -> k = {k.n:.5f} ± {k.s:.5f}")
        return k.n, k.s

    print(f"\n=== Loading search: {name} [{knob}, {profile}] target k = {k_target} ===")
    res = secant_search(evaluate, k_target, tol_pcm=tol_pcm, max_iter=max_iter, bounds=bounds)

    searched = copy.deepcopy(doc)
    searched["assembly_types"] = apply_loading_scale(doc["assembly_types"], res.x)
    searched["loading_search"] = {
        "source_loadout": os.path.relpath(loadout_path),
        "knob": knob,
        "scale": res.x,
        "k_target": k_target,
        "keff": res.k,
        "keff_unc": res.k_unc,
        "converged": res.converged,
        "profile": profile,
        "searched_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }
    result = {"result": res.to_dict(), "loadout": searched}
    (run / "search_result.json").write_text(json.dumps(result, indent=2))

    status = "converged" if res.converged else "NOT converged"
    print(f"  {status}: scale {res.x:.5f}, k = {res.k:.5f} ± {res.k_unc:.5f} "
          f"after {len(res.history)} solves -> {run / 'search_result.json'}")
    if res.k > K_MAX:
        print(f"  WARNING: k above the T7.4 limit {K_MAX}")
    return result


if __name__ == "__main__":
    import argparse
    import sys

    from build_prism_core import RUN_PROFILES

    ap = argparse.ArgumentParser(description="BOC fissile-loading search for a fresh core.")
    ap.add_argument("loadout")
    ap.add_argument("--k-target", type=float, default=K_TARGET)
    ap.add_argument("--profile", choices=RUN_PROFILES, default="quick")
    ap.add_argument("--tol-pcm", type=float, default=100.0)
    ap.add_argument("--max-iter", type=int, default=8)
    ap.add_argument("--run-dir")
    ap.add_argument("--write-loadout", help="also write the searched loadout JSON here")
    opts = ap.parse_args()

    result = search_boc_loading(
        opts.loadout, k_target=opts.k_target, profile=opts.profile,
        run_dir=opts.run_dir, tol_pcm=opts.tol_pcm, max_iter=opts.max_iter,
    )
    if opts.write_loadout:
        Path(opts.write_loadout).write_text(json.dumps(result["loadout"], indent=2) + "\n")
        print(f"  searched loadout -> {opts.write_loadout}")
    sys.exit(0 if result["result"]["converged"] else 1)
