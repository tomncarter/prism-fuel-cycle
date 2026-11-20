# fuel_cycle/search.py
"""
Fissile-loading search: the noise-tolerant root finder and the loadout knobs.

Pure Python — no OpenMC dependency.  The transport side lives in
fuel_cycle.loading_search (fresh cores) and fuel_cycle.cycle_runner (reloads).

The search variable is always a dimensionless scale on the startup loading:
  - Pu-fuelled drivers:  pu_hm_frac_zone = scale × startup pu_hm_frac_zone
                         (zone-to-zone ratio preserved)
  - U-fuelled drivers:   enriched components of u_vector scaled by `scale`,
                         depleted U fills the remainder of the U blend
"""
from __future__ import annotations

import copy
from dataclasses import dataclass, field, asdict
from typing import Callable

# Driver assembly universes -> zone-default material key (see materials.tier1_homog)
DRIVER_UNIVERSES = {
    "u_inner_fuel": "fuel_driver_inner_homog",
    "u_outer_fuel": "fuel_driver_outer_homog",
}

# u_vector components treated as fertile filler rather than fissile feed
FILLER_U_VECTORS = frozenset({"depleted", "natural", "u238"})

PU_HM_MAX = 0.35   # same clamp as materials.tier1_homog.make_smeared_driver


# ---------------------------------------------------------------------------
# Root finder
# ---------------------------------------------------------------------------

@dataclass
class SearchPoint:
    x: float
    k: float
    k_unc: float


@dataclass
class SearchResult:
    x: float
    k: float
    k_unc: float
    converged: bool
    history: list = field(default_factory=list)   # list[SearchPoint]

    def to_dict(self) -> dict:
        return asdict(self)


def secant_search(
    evaluate: Callable[[float], tuple[float, float]],
    target: float,
    x0: float = 1.0,
    tol_pcm: float = 100.0,
    max_iter: int = 8,
    bounds: tuple[float, float] = (0.2, 2.0),
    worth: float = 0.6,
) -> SearchResult:
    """
    Find x with k(x) = target, where evaluate(x) -> (k, 1σ) is a Monte Carlo
    eigenvalue.  A point is accepted when |k - target| ≤ max(tol_pcm, 2σ).

    The first step assumes a fissile worth of dk/k = worth × dx/x (~0.5–0.7 in
    fast cores); later steps are secants through the two latest points.  A
    non-positive secant slope (statistical noise swamping a small step) falls
    back to the worth estimate.

    Returns the accepted point, or the closest point with converged=False if
    max_iter runs out or the step is pinned at a bound.
    """
    lo, hi = bounds
    clamp = lambda x: min(max(x, lo), hi)
    history: list[SearchPoint] = []

    def run(x: float) -> SearchPoint:
        k, s = evaluate(x)
        p = SearchPoint(x, k, s)
        history.append(p)
        return p

    def accepted(p: SearchPoint) -> bool:
        return abs(p.k - target) <= max(tol_pcm * 1e-5, 2.0 * p.k_unc)

    def worth_step(p: SearchPoint) -> float:
        return clamp(p.x * (1.0 + (target - p.k) / (p.k * worth)))

    def done(p: SearchPoint, ok: bool) -> SearchResult:
        return SearchResult(p.x, p.k, p.k_unc, ok, history)

    prev = run(clamp(x0))
    if accepted(prev):
        return done(prev, True)
    x = worth_step(prev)

    for _ in range(max_iter - 1):
        if x == prev.x:   # pinned at a bound
            break
        cur = run(x)
        if accepted(cur):
            return done(cur, True)
        slope = (cur.k - prev.k) / (cur.x - prev.x)
        x = clamp(cur.x + (target - cur.k) / slope) if slope > 0.0 else worth_step(cur)
        prev = cur

    best = min(history, key=lambda p: abs(p.k - target))
    return done(best, False)


# ---------------------------------------------------------------------------
# Loadout knobs
# ---------------------------------------------------------------------------

def driver_types(asm_types: dict) -> list[str]:
    """Assembly type names whose universe is a driver fuel zone."""
    return sorted(t for t, at in asm_types.items()
                  if at.get("universe") in DRIVER_UNIVERSES)


def loading_knob(asm_types: dict) -> str:
    """'tru' if any driver carries Pu, else 'u_fissile'."""
    pu = [asm_types[t].get("fuel_params", {}).get("pu_hm_frac", 0.0)
          for t in driver_types(asm_types)]
    return "tru" if any(f > 0.0 for f in pu) else "u_fissile"


def _u_components(fp: dict) -> dict:
    uv = fp.get("u_vector", "u238")
    return {uv: 1.0} if isinstance(uv, str) else dict(uv)


def scale_bounds(asm_types: dict) -> tuple[float, float]:
    """Physical range of the scale: Pu/HM ≤ PU_HM_MAX, or enriched U weight ≤ 1."""
    limits = []
    for t in driver_types(asm_types):
        fp = asm_types[t].get("fuel_params", {})
        if loading_knob(asm_types) == "tru":
            base = fp.get("pu_hm_frac", 0.0)
            if base > 0.0:
                limits.append(PU_HM_MAX / base)
        else:
            fissile = sum(w for n, w in _u_components(fp).items() if n not in FILLER_U_VECTORS)
            if fissile > 0.0:
                limits.append(1.0 / fissile)
    hi = min([2.0] + limits)
    return (0.2, hi)


def apply_loading_scale(asm_types: dict, scale: float) -> dict:
    """
    Copy of asm_types with every driver's fissile loading scaled.  Non-driver
    types (blankets, reflector, ...) are returned unchanged.
    """
    out = copy.deepcopy(asm_types)
    knob = loading_knob(asm_types)
    for t in driver_types(asm_types):
        fp = out[t].setdefault("fuel_params", {})
        if knob == "tru":
            if fp.get("pu_hm_frac", 0.0) > 0.0:
                fp["pu_hm_frac"] = fp["pu_hm_frac"] * scale
            continue
        comps = _u_components(fp)
        fissile = {n: w * scale for n, w in comps.items() if n not in FILLER_U_VECTORS}
        total = sum(fissile.values())
        if total > 1.0 + 1e-12:
            raise ValueError(f"{t}: scale {scale:.4f} puts enriched U above 100% of the blend")
        filler = next((n for n in comps if n in FILLER_U_VECTORS), "depleted")
        fissile[filler] = max(1.0 - total, 0.0)
        fp["u_vector"] = fissile
    return out
