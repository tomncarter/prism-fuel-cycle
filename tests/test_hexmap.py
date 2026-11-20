"""The (q, r) frame used by loadouts must match where OpenMC places each lattice cell."""
import openmc
import pytest

from geometry.hexmap import axial_to_xy, build_hexlattice_rings_from_assignment, full_hex_coords

PITCH = 16.142
MAX_RING = 5


def _lattice_with_unique_universes():
    """One distinct universe per (q, r), so a point lookup identifies the slot."""
    coords = full_hex_coords(MAX_RING)
    univs = {}
    assignment = {}
    for q, r in coords:
        name = f"u_{q}_{r}"
        cell = openmc.Cell(name=name, fill=None)
        univs[name] = openmc.Universe(name=name, cells=[cell])
        assignment[(q, r)] = name
    rings, _ = build_hexlattice_rings_from_assignment(assignment, MAX_RING)
    lat = openmc.HexLattice()
    lat.center = (0.0, 0.0)
    lat.pitch = (PITCH,)
    lat.universes = [[univs[n] for n in ring] for ring in rings]
    root = openmc.Cell(fill=lat, region=-openmc.ZCylinder(r=PITCH * (MAX_RING + 1)))
    return openmc.Geometry(openmc.Universe(cells=[root])), coords


def test_every_slot_lands_at_axial_to_xy():
    geom, coords = _lattice_with_unique_universes()
    wrong = []
    for q, r in coords:
        x, y = axial_to_xy(q, r, PITCH)
        found = geom.find((x, y, 0.0))[-1].name
        if found != f"u_{q}_{r}":
            wrong.append(((q, r), found))
    assert not wrong, f"{len(wrong)} slots misplaced, e.g. {wrong[:5]}"


@pytest.mark.parametrize("qr, xy", [
    ((0, 3), (0.0, 3.0)),                      # straight up
    ((0, -3), (0.0, -3.0)),                    # straight down
    ((3, 0), (3 * 0.8660254, 1.5)),            # upper-right ring corner
    ((-3, 0), (-3 * 0.8660254, -1.5)),         # lower-left ring corner
])
def test_frame_orientation(qr, xy):
    x, y = axial_to_xy(*qr, 1.0)
    assert x == pytest.approx(xy[0]) and y == pytest.approx(xy[1])
