# geometry/params.py
from dataclasses import dataclass
import math
from typing import Tuple

@dataclass(frozen=True)
class Prism3D:
    active_height_cm: float = 101.60  # Tier-1 basis; later swap to 66.0 / 106.7 PRISM missions

@dataclass(frozen=True)
class AssemblyGeom:
    # From ABR/S-PRISM-style basis:
    assembly_pitch_cm: float = 16.142
    inter_duct_gap_cm: float = 0.432
    duct_thickness_cm: float = 0.394

    @property
    def duct_flat_to_flat_cm(self) -> float:
        # Pitch = duct_flat_to_flat + gap
        return self.assembly_pitch_cm - self.inter_duct_gap_cm

    @property
    def duct_inner_flat_to_flat_cm(self) -> float:
        return self.duct_flat_to_flat_cm - 2.0 * self.duct_thickness_cm

    @property
    def duct_edge_cm(self) -> float:
        # edge length a where flat-to-flat F = sqrt(3)*a
        return self.duct_flat_to_flat_cm / math.sqrt(3.0)

    @property
    def duct_inner_edge_cm(self) -> float:
        return self.duct_inner_flat_to_flat_cm / math.sqrt(3.0)


@dataclass(frozen=True)
class VesselGeom:
    # Radial
    lattice_bound_r_cm: float = 195.0    # cylinder fully enclosing hex lattice (~187 cm outer edge + margin)
    ihx_pool_outer_r_cm: float = 240.0   # inner face of HT9 vessel wall (IHX/pool fills lattice_bound → here)
    vessel_wall_cm: float = 4.0          # HT9 cylindrical vessel wall thickness

    # Axial plenums (outside active zone)
    lower_plenum_h_cm: float = 50.0      # below active zone (more hardware: pumps, structure)
    upper_plenum_h_cm: float = 35.0      # above active zone
    end_plate_cm: float = 2.0            # HT9 top/bottom end plates

    @property
    def vessel_outer_r_cm(self) -> float:
        return self.ihx_pool_outer_r_cm + self.vessel_wall_cm

    def axial_bounds(self, active_height_cm: float) -> Tuple[float, float]:
        """Vacuum boundary z positions (bottom, top)."""
        H = active_height_cm
        z_bot = -(H / 2.0 + self.lower_plenum_h_cm + self.end_plate_cm)
        z_top = +(H / 2.0 + self.upper_plenum_h_cm + self.end_plate_cm)
        return z_bot, z_top
