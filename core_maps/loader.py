# core_maps/loader.py
"""
Core loadout loader.

Reads a JSON loadout file and returns a list of AssemblyPosition objects — one
per lattice cell — plus the raw assembly_types dict for downstream use by the
material factory.

Two layout modes are supported:

  "rings"     — compact ring-range notation, auto-generates material_ids
  "positions" — explicit per-position list (Tier-2 / post-shuffle)

Ring mode expands to the same internal representation as position mode, so all
downstream code works with AssemblyPosition lists regardless of input format.
"""
from __future__ import annotations

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


def load_core(path: str) -> tuple[list[AssemblyPosition], dict]:
    """
    Load a core loadout JSON file.

    Returns:
        positions  — list of AssemblyPosition, one per lattice cell
        asm_types  — raw assembly_types dict (passed to material factory)
    """
    with open(path) as f:
        doc = json.load(f)

    asm_types: dict = doc["assembly_types"]
    layout: dict = doc["layout"]
    mode: str = layout["mode"]

    positions: list[AssemblyPosition] = []

    if mode == "rings":
        # Build ring -> type lookup from the range entries
        ring_type: dict[int, str] = {}
        for entry in layout["rings"]:
            lo, hi = entry["range"]
            t = entry["type"]
            for ring in range(lo, hi + 1):
                ring_type[ring] = t

        max_ring = max(ring_type)
        for q, r in full_hex_coords(max_ring):
            d = hex_distance((q, r))
            type_name = ring_type.get(d)
            if type_name is None:
                continue
            at = asm_types[type_name]
            positions.append(AssemblyPosition(
                q=q, r=r,
                type=type_name,
                universe=at["universe"],
                material_id=_auto_material_id(type_name, q, r),
                fuel_params=at.get("fuel_params", {}),
            ))

    elif mode == "rings_with_overrides":
        # Base ring assignment
        ring_type: dict[int, str] = {}
        for entry in layout["rings"]:
            lo, hi = entry["range"]
            t = entry["type"]
            for ring in range(lo, hi + 1):
                ring_type[ring] = t

        # Per-position overrides: (q, r) -> type_name
        # Entries without "q"/"r" are treated as comments and skipped.
        overrides: dict[Axial, str] = {
            (e["q"], e["r"]): e["type"]
            for e in layout.get("overrides", [])
            if "q" in e and "r" in e
        }

        max_ring = max(ring_type)
        for q, r in full_hex_coords(max_ring):
            d = hex_distance((q, r))
            type_name = overrides.get((q, r), ring_type.get(d))
            if type_name is None:
                continue
            at = asm_types[type_name]
            positions.append(AssemblyPosition(
                q=q, r=r,
                type=type_name,
                universe=at["universe"],
                material_id=_auto_material_id(type_name, q, r),
                fuel_params=at.get("fuel_params", {}),
            ))

    elif mode == "positions":
        for entry in layout["positions"]:
            q, r = entry["q"], entry["r"]
            type_name = entry["type"]
            at = asm_types[type_name]
            positions.append(AssemblyPosition(
                q=q, r=r,
                type=type_name,
                universe=at["universe"],
                material_id=entry.get(
                    "material_id", _auto_material_id(type_name, q, r)
                ),
                fuel_params=at.get("fuel_params", {}),
            ))

    else:
        raise ValueError(f"Unknown layout mode: {mode!r}. Expected 'rings' or 'positions'.")

    return positions, asm_types
