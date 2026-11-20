# core_maps/abr1000_breakeven_assignment.py
"""
ABR-1000 Breakeven core (CR ≈ 1.0).

Three fuel zones, no blankets:
  Rings 0-2   — inner core
  Rings 3-5   — middle core
  Rings 6-8   — outer core
  Rings 9-10  — reflector
  Ring 11     — shield

Control rods, ultimate shutdown rods, and GEMs are at the same (q, r)
positions as the PRISM CR=0.72 design.

Corner notch (each 60° segment):
  Ring 11 corner vertices are empty (no assembly):
      (0, 11)  (11, 0)  (11, -11)  (0, -11)  (-11, 0)  (-11, 11)
  The directly-inward ring-10 vertices are promoted to shield:
      (0, 10)  (10, 0)  (10, -10)  (0, -10)  (-10, 0)  (-10, 10)
  Ring-8 corner vertices (outer fuel by default) are promoted to reflector:
      (0, 8)  (8, 0)  (8, -8)  (0, -8)  (-8, 0)  (-8, 8)
"""
from __future__ import annotations
from typing import Dict, Tuple

Axial = Tuple[int, int]

CONTROL_ROD_POSITIONS: list[Axial] = [
    # Ring 3 — safety (primary shutdown) rods
    (0, 3), (3, -3), (-3, 0),
    # Ring 6 — primary control rods
    (-3, 6), (3, 3), (6, -3), (3, -6), (-6, 3), (-3, -3),
]

ULTIMATE_SHUTDOWN_POSITIONS: list[Axial] = [
    # Ring 3 — ultimate shutdown rods (interleaved with safety rods)
    (-3, 3), (3, 0), (0, -3),
]

GEM_POSITIONS: list[Axial] = [
    # Ring 7 — Gas Expansion Modules
    (0, 7), (7, 0), (7, -7), (0, -7), (-7, 0), (-7, 7),
]

# Ring-11 corner vertices — no assembly at these positions
EMPTY_CORNER_POSITIONS: list[Axial] = [
    (0, 11), (11, 0), (11, -11), (0, -11), (-11, 0), (-11, 11),
]

# Ring-10 vertices directly inward from empty corners — promoted to shield
CORNER_SHIELD_POSITIONS: list[Axial] = [
    (0, 10), (10, 0), (10, -10), (0, -10), (-10, 0), (-10, 10),
]

# Ring-8 corner vertices — promoted to reflector (outer fuel by default)
CORNER_REFLECTOR_POSITIONS: list[Axial] = [
    (0, 8), (8, 0), (8, -8), (0, -8), (-8, 0), (-8, 8),
]


def abr1000_breakeven_assignment() -> Dict[Axial, str]:
    """ABR-1000 breakeven zoning with 3 fuel zones, no blankets."""
    def dist(c: Axial) -> int:
        q, r = c
        return max(abs(q), abs(r), abs(-q - r))

    assignment: Dict[Axial, str] = {}

    for q in range(-11, 12):
        for r in range(-11, 12):
            if not (-11 <= q + r <= 11):
                continue
            d = dist((q, r))
            if d == 11:
                assignment[(q, r)] = "u_shield"
            elif d in (9, 10):
                assignment[(q, r)] = "u_reflector"
            elif d in (6, 7, 8):
                assignment[(q, r)] = "u_outer_fuel"
            elif d in (3, 4, 5):
                assignment[(q, r)] = "u_middle_fuel"
            else:  # d in (0, 1, 2)
                assignment[(q, r)] = "u_inner_fuel"

    for pos in CONTROL_ROD_POSITIONS:
        assignment[pos] = "u_control_rod"
    for pos in ULTIMATE_SHUTDOWN_POSITIONS:
        assignment[pos] = "u_ultimate_shutdown"
    for pos in GEM_POSITIONS:
        assignment[pos] = "u_gem"
    for pos in CORNER_SHIELD_POSITIONS:
        assignment[pos] = "u_shield"
    for pos in CORNER_REFLECTOR_POSITIONS:
        assignment[pos] = "u_reflector"
    for pos in EMPTY_CORNER_POSITIONS:
        assignment[pos] = "u_empty"

    return assignment
