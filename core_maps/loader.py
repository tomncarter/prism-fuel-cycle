# core_maps/loader.py
"""
Core loadout loader.

Reads a JSON loadout file and returns a list of AssemblyPosition objects — one
per lattice cell — plus the raw assembly_types dict for downstream use by the
material factory.

Two layout modes are supported:

  "rings"                — ring-range notation only
  "rings_with_overrides" — rings, then per-position overrides (rods, GEMs, notches);
                           an override can carry "sym": 3 or 6 to place all its
                           rotational images at once
  "positions"            — explicit per-position list (Tier-2 / post-shuffle)

Ring mode expands to the same internal representation as position mode, so all
downstream code works with AssemblyPosition lists regardless of input format.
"""
from __future__ import annotations

import copy
import json
from dataclasses import dataclass, field
from typing import Dict, List, Tuple

from geometry.hexmap import full_hex_coords, hex_distance
from geometry.params import Prism3D

Axial = Tuple[int, int]


@dataclass
class AssemblyPosition:
    q: int
    r: int
    type: str        # key into assembly_types
    universe: str    # resolved from assembly_types[type]["universe"]
    material_id: str # stable per-assembly identifier (Tier-2 compatible)
    fuel_params: dict = field(default_factory=dict)  # {} for non-fuel types


def _auto_material_id(type_name: str, q: int, r: int) -> str:
    """Stable, human-readable material ID derived from type and coordinates."""
    return f"{type_name}_q{q:+d}_r{r:+d}"


def load_axial(path: str) -> Prism3D:
    """Axial geometry of a loadout: its "active_height_cm", else the Prism3D default."""
    with open(path) as f:
        doc = json.load(f)
    if "active_height_cm" in doc:
        return Prism3D(active_height_cm=float(doc["active_height_cm"]))
    return Prism3D()


def rotate60(q: int, r: int, k: int = 1) -> Axial:
    """Rotate an axial coordinate by k x 60 degrees about the core centre."""
    for _ in range(k % 6):
        q, r = -r, q + r
    return q, r


def _expand_override(e: dict) -> list[Axial]:
    """Positions named by one override entry.

    {"q": 0, "r": 3}             -> that one position
    {"q": 0, "r": 3, "sym": 3}   -> it and its 120/240-degree images (3 positions)
    {"q": 0, "r": 7, "sym": 6}   -> all six 60-degree images
    {"positions": [[q, r], ...]} -> an explicit list (each may also take "sym")
    """
    sym = int(e.get("sym", 1))
    if sym not in (1, 2, 3, 6):
        raise ValueError(f"override 'sym' must be 1, 2, 3 or 6, got {sym}: {e}")
    seeds = [tuple(p) for p in e["positions"]] if "positions" in e else [(e["q"], e["r"])]
    out: list[Axial] = []
    for q, r in seeds:
        for i in range(sym):
            p = rotate60(q, r, i * 6 // sym)
            if p not in out:
                out.append(p)
    return out


def _ring_types(rings: list[dict]) -> dict[int, str]:
    ring_type: dict[int, str] = {}
    for entry in rings:
        lo, hi = entry["range"]
        for ring in range(lo, hi + 1):
            if ring in ring_type:
                raise ValueError(f"ring {ring} assigned twice ({ring_type[ring]!r}, {entry['type']!r})")
            ring_type[ring] = entry["type"]
    missing = sorted(set(range(max(ring_type) + 1)) - set(ring_type))
    if missing:
        raise ValueError(f"rings {missing} have no type")
    return ring_type


def resolve_layout(layout: dict) -> dict[Axial, str]:
    """(q, r) -> assembly-type name for a "rings" / "rings_with_overrides" layout.

    Rings give every position a base type; overrides then replace individual
    positions in file order.  An override may not land outside the ring map or
    hit a position an earlier override already set.
    """
    ring_type = _ring_types(layout["rings"])
    max_ring = max(ring_type)
    types = {(q, r): ring_type[hex_distance((q, r))] for q, r in full_hex_coords(max_ring)}

    set_by: dict[Axial, dict] = {}
    for e in layout.get("overrides", []):
        if "type" not in e:          # {"comment": ...} entries
            continue
        for pos in _expand_override(e):
            if pos not in types:
                raise ValueError(f"override {pos} is outside ring {max_ring}: {e}")
            if pos in set_by:
                raise ValueError(f"position {pos} overridden twice: {set_by[pos]} and {e}")
            set_by[pos] = e
            types[pos] = e["type"]
    return types


def type_counts(positions: list[AssemblyPosition]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for p in positions:
        counts[p.type] = counts.get(p.type, 0) + 1
    return counts


PIN_KEYS = ("pin_design", "volume_fractions")


def _with_pin_design(doc: dict) -> dict:
    """assembly_types with a top-level "pin_design" / "volume_fractions" copied into
    every fuel_params that doesn't set its own (fuel and blanket types only)."""
    asm_types = copy.deepcopy(doc["assembly_types"])
    top = {k: doc[k] for k in PIN_KEYS if k in doc}
    if top:
        for at in asm_types.values():
            fp = at.get("fuel_params")
            if fp is not None and not any(k in fp for k in PIN_KEYS):
                fp.update(copy.deepcopy(top))
    return asm_types


def load_core(path: str) -> tuple[list[AssemblyPosition], dict]:
    """
    Load a core loadout JSON file.

    Returns:
        positions  — list of AssemblyPosition, one per lattice cell
        asm_types  — raw assembly_types dict (passed to material factory)

    If the loadout has "expected_counts" ({type: n}), the loaded core must match
    it exactly; this pins each loadout to its source (e.g. Triplett Table III).
    """
    with open(path) as f:
        doc = json.load(f)

    asm_types: dict = _with_pin_design(doc)
    layout: dict = doc["layout"]
    mode: str = layout["mode"]

    if mode in ("rings", "rings_with_overrides"):
        if mode == "rings" and layout.get("overrides"):
            raise ValueError("'rings' mode has overrides; use mode 'rings_with_overrides'")
        entries = [
            {"q": q, "r": r, "type": t} for (q, r), t in resolve_layout(layout).items()
        ]
    elif mode == "positions":
        entries = layout["positions"]
    else:
        raise ValueError(
            f"Unknown layout mode: {mode!r}. Expected 'rings', 'rings_with_overrides' or 'positions'."
        )

    positions: list[AssemblyPosition] = []
    for entry in entries:
        q, r, type_name = entry["q"], entry["r"], entry["type"]
        if type_name not in asm_types:
            raise ValueError(f"{path}: type {type_name!r} at ({q}, {r}) is not in assembly_types")
        at = asm_types[type_name]
        positions.append(AssemblyPosition(
            q=q, r=r,
            type=type_name,
            universe=at["universe"],
            material_id=entry.get("material_id", _auto_material_id(type_name, q, r)),
            fuel_params=at.get("fuel_params", {}),
        ))

    expected = doc.get("expected_counts")
    if expected:
        got = type_counts(positions)
        diff = {t: (expected.get(t, 0), got.get(t, 0))
                for t in set(expected) | set(got) if expected.get(t, 0) != got.get(t, 0)}
        if diff:
            raise ValueError(f"{path}: assembly counts differ from expected_counts "
                             f"(type: expected, got): {diff}")

    return positions, asm_types
