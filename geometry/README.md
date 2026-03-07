# geometry/

Three modules that together convert a ring-based core map into an OpenMC hex lattice.

## params.py

Frozen dataclasses holding all geometric constants. Nothing is hardcoded elsewhere —
`build_prism_core.py` instantiates these and passes them down.

- **`AssemblyGeom`** — assembly pitch, inter-duct gap, duct wall thickness. Derived
  properties compute flat-to-flat and edge lengths for both the outer duct surface and
  the inner duct surface (the boundary of the fuel/blanket interior region).
- **`Prism3D`** — active core height.
- **`VesselGeom`** — radial and axial extent of the vessel structure (lattice boundary
  cylinder, IHX pool annulus, vessel wall, upper/lower plenums, end plates). Also
  computes the vacuum boundary z-positions as a function of active height.

## universe_factory.py

`UniverseFactory` builds the reusable OpenMC `Universe` objects that fill each lattice
position. Takes an `AssemblyGeom` and a materials dict as constructor arguments.

Each assembly universe has three concentric cells:

```
  sodium background    (outside duct_outer, fills lattice voids)
  HT9 duct wall        (between duct_outer and duct_inner)
  interior region      (inside duct_inner — fuel, blanket, reflector, etc.)
```

`build_all()` returns a dict keyed by universe name string (e.g. `"u_inner_fuel"`),
matching the keys used in `core_maps/`. Currently defined types:

| Key | Interior material |
|-----|------------------|
| `u_inner_fuel` | Smeared U-Pu-Zr + Na + steel (inner Pu fraction) |
| `u_outer_fuel` | Smeared U-Pu-Zr + Na + steel (outer Pu fraction) |
| `u_internal_blanket` | Smeared NU or DU + Na + steel |
| `u_radial_blanket` | Smeared NU or DU + Na + steel |
| `u_reflector` | HT9/steel |
| `u_shield` | HT9/steel |
| `u_control` | B4C (normal shutdown) |
| `u_ult_shutdown` | B4C (ultimate shutdown, higher density) |
| `u_gem` | Sodium (GEM placeholder — not yet activated) |
| `u_source` | Sodium (neutron source placeholder) |

## hexmap.py

Coordinate utilities for pointy-top axial hex grids and lattice construction.

- **`hex_distance(q, r)`** — Chebyshev distance from origin; equals the ring number.
- **`axial_to_xy(q, r, pitch)`** — converts axial coordinates to Cartesian x, y.
- **`full_hex_coords(max_ring)`** — all axial coordinates within a hex of radius
  `max_ring`.
- **`ring_angle_sort(coords, start_angle_deg=90)`** — sorts a ring's coordinates
  clockwise starting from `start_angle_deg` (default top = 90°). Matches OpenMC's
  `HexLattice.universes` convention (outermost ring first, clockwise from top).
- **`build_hexlattice_rings_from_assignment(assignment, max_ring, ...)`** — top-level
  function called by `build_prism_core.py`. Takes the `(q,r) → universe_name` dict
  from `core_maps/`, applies ring sorting, and returns the nested list of universe
  name strings that OpenMC expects.

### Coordinate convention

Pointy-top axial: `x = √3·(q + r/2)`, `y = 3/2·r`. Ring 2 positions in clockwise
order from top: `(-1,2)` TOP → `(0,2)` → `(1,1)` → `(2,0)` RIGHT → `(2,-1)` →
`(2,-2)` → `(1,-2)` BOT → `(0,-2)` → `(-1,-1)` → `(-2,0)` LEFT → `(-2,1)` →
`(-2,2)` → back to top.
