# geometry/hexmap.py
from __future__ import annotations
import math
from collections import defaultdict
from typing import Dict, Tuple, List

Axial = Tuple[int, int]

def hex_distance(qr: Axial) -> int:
    q, r = qr
    return max(abs(q), abs(r), abs(-q - r))

def axial_to_xy(q: int, r: int, pitch: float = 1.0) -> Tuple[float, float]:
    # pointy-top axial to x,y
    x = pitch * (math.sqrt(3) * (q + r/2.0))
    y = pitch * (1.5 * r)
    return x, y

def ring_angle_sort(coords: List[Axial], start_angle_deg: float = 90.0) -> List[Axial]:
    # clockwise, starting at "top" (90°)
    def ang(c: Axial) -> float:
        x, y = axial_to_xy(c[0], c[1], 1.0)
        a = math.degrees(math.atan2(y, x))
        a = (a + 360.0) % 360.0
        a = (a - start_angle_deg + 360.0) % 360.0
        return a
    return sorted(coords, key=ang)

def full_hex_coords(max_ring: int) -> List[Axial]:
    # full hex of radius max_ring (includes center), axial constraint form
    coords = []
    R = max_ring
    for q in range(-R, R+1):
        for r in range(-R, R+1):
            if -R <= q + r <= R:
                coords.append((q, r))
    return coords

def build_hexlattice_rings_from_assignment(
    assignment: Dict[Axial, str],
    max_ring: int,
    default_universe: str = "u_reflector",
    start_angle_deg: float = 90.0,
) -> Tuple[List[List[str]], Axial]:
    """
    Returns:
      - rings: list of rings outer->inner, plus center as last list (length 1)
      - center coord
    assignment maps (q,r) -> universe-name string (e.g., 'u_inner_fuel').

    Any missing coordinates in the full hex are filled with default_universe.
    """
    coords = full_hex_coords(max_ring)

    by_ring = defaultdict(list)
    for c in coords:
        by_ring[hex_distance(c)].append(c)

    rings: List[List[str]] = []
    for ring in range(max_ring, 0, -1):
        ordered = ring_angle_sort(by_ring[ring], start_angle_deg=start_angle_deg)
        rings.append([assignment.get(c, default_universe) for c in ordered])

    center = (0, 0)
    rings.append([assignment.get(center, default_universe)])
    return rings, center
