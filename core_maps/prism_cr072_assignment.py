# core_maps/prism_cr072_assignment.py
"""
PRISM UNF disposal core, CR=0.72.

Starts from the ABR-like ring zoning and overrides specific positions with
control rod / GEM types as the design is built up.

Control rod locations (Tier-2 design):
  Ring 3 (inner driver)  — 3 safety/shutdown rods:
      (0, 3)  (3, -3)  (-3, 0)
  Ring 3 (inner driver)  — 3 ultimate shutdown rods (interleaved):
      (-3, 3)  (3, 0)  (0, -3)
  Ring 6 (outer driver)  — 6 primary control rods:
      (-3, 6)  (3, 3)  (6, -3)  (3, -6)  (-6, 3)  (-3, -3)
  Ring 7 (outer driver)  — 6 GEMs:
      (0, 7)  (7, 0)  (7, -7)  (0, -7)  (-7, 0)  (-7, 7)

Corner notch (each 60° segment):
  Ring 11 corner vertices are empty (no assembly):
      (0, 11)  (11, 0)  (11, -11)  (0, -11)  (-11, 0)  (-11, 11)
  The directly-inward ring-10 vertices are promoted to shield:
      (0, 10)  (10, 0)  (10, -10)  (0, -10)  (-10, 0)  (-10, 10)
"""
from __future__ import annotations
from typing import Dict, Tuple

from core_maps.abr_like_assignment import abr_like_assignment

Axial = Tuple[int, int]

CONTROL_ROD_POSITIONS: list[Axial] = [
    # Ring 3 — safety rods
    (0, 3), (3, -3), (-3, 0),
    # Ring 6 — control rods
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

# Ring-8 vertices — promoted to reflector (reflector moves inward one ring at corners)
CORNER_REFLECTOR_POSITIONS: list[Axial] = [
    (0, 8), (8, 0), (8, -8), (0, -8), (-8, 0), (-8, 8),
]


def prism_cr072_assignment() -> Dict[Axial, str]:
    """ABR-like zoning with CR=0.72 control rod, ultimate shutdown, GEM, and corner notch."""
    assignment = abr_like_assignment()
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
