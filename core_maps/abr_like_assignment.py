# core_maps/abr_like_assignment.py
from __future__ import annotations
from typing import Dict, Tuple

Axial = Tuple[int, int]

def abr_like_assignment() -> Dict[Axial, str]:
    """
    S-PRISM-like radial zoning for Tier-1 waste-elimination study.

    Radial layout (assembly counts):
      Ring 0        :  1  inner fuel       (center)
      Rings 1-3     : 36  inner fuel
      Ring 4        : 24  internal blanket  (NU/DU, between driver zones)
      Rings 5-7     :108  outer fuel
      Ring 8        : 48  radial blanket    (NU/DU, captures radial leakage)
      Rings 9-10    :114  reflector
      Ring 11       : 66  shield
    """
    def dist(c: Axial) -> int:
        q, r = c
        return max(abs(q), abs(r), abs(-q-r))

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
            elif d == 8:
                assignment[(q, r)] = "u_radial_blanket"
            elif d in (5, 6, 7):
                assignment[(q, r)] = "u_outer_fuel"
            elif d == 4:
                assignment[(q, r)] = "u_internal_blanket"
            else:  # d in (0, 1, 2, 3)
                assignment[(q, r)] = "u_inner_fuel"

    return assignment
