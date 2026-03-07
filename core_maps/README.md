# core_maps/

Core loading maps: functions that return `Dict[(q,r), universe_name]` — a mapping
from every axial hex coordinate to the assembly type that should occupy it.

`build_prism_core.py` calls one of these functions and passes the result to
`geometry/hexmap.py::build_hexlattice_rings_from_assignment()`, which sorts it into
the ring-ordered list OpenMC expects.

## abr_like_assignment.py

S-PRISM / ABR-like radial zoning with 12 rings (center + rings 1–11):

```
Ring 0        :   1  inner fuel        (center assembly)
Rings 1–3     :  36  inner fuel        (37 total inner driver)
Ring 4        :  24  internal blanket  (NU or DU; surrounded by driver on both sides)
Rings 5–7     : 108  outer fuel        (outer driver zone)
Ring 8        :  48  radial blanket    (NU or DU; captures outward-leaking neutrons)
Rings 9–10    : 114  reflector         (HT9 steel)
Ring 11       :  66  shield            (HT9 steel)
              ──────
Total         : 397  assemblies
```

The internal blanket (ring 4) sits between the inner and outer driver zones, so it
sees high flux from both sides and accumulates Pu-239 quickly. The radial blanket
(ring 8) sees softer, lower flux and accumulates Pu-239 more slowly.

## Adding new maps

Add a new `*_assignment.py` file following the same interface:

```python
def my_assignment() -> Dict[Tuple[int,int], str]:
    ...
    return assignment  # keys are (q,r) axial coords, values are universe name strings
```

Universe name strings must match keys returned by `UniverseFactory.build_all()`:
`u_inner_fuel`, `u_outer_fuel`, `u_internal_blanket`, `u_radial_blanket`,
`u_reflector`, `u_shield`, `u_control`, `u_ult_shutdown`, `u_gem`, `u_source`.

Any coordinate in the full hex that is missing from the dict defaults to
`u_reflector` (configurable in `build_hexlattice_rings_from_assignment`).
