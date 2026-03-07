# geometry/universe_factory.py

import openmc
from geometry.params import AssemblyGeom

class UniverseFactory:
    def __init__(self, geom: AssemblyGeom, mats: dict):
        self.g = geom
        self.m = mats

    def build_all(self):
        return {
            "u_inner_fuel": self._fuel_assy("inner_fuel", self.m["fuel_driver_inner_homog"]),
            "u_outer_fuel": self._fuel_assy("outer_fuel", self.m["fuel_driver_outer_homog"]),
            "u_internal_blanket": self._fuel_assy("internal_blanket", self.m["fuel_blanket_internal_homog"]),
            "u_radial_blanket": self._fuel_assy("radial_blanket", self.m["fuel_blanket_radial_homog"]),
            "u_reflector": self._solid_assy("reflector", self.m["reflector"]),
            "u_shield": self._solid_assy("shield", self.m["shield"]),
            "u_control": self._solid_assy("control", self.m["b4c"]),
            "u_ult_shutdown": self._solid_assy("ult_shutdown", self.m["b4c_ult"]),
            "u_gem": self._solid_assy("gem", self.m["sodium"]),
            "u_source": self._solid_assy("source", self.m["sodium"]),
        }

    def _fuel_assy(self, name: str, interior_mat: openmc.Material) -> openmc.Universe:
        u = openmc.Universe(name=f"assy_{name}")

        duct_outer = openmc.model.HexagonalPrism(edge_length=self.g.duct_edge_cm, orientation="x")
        duct_inner = openmc.model.HexagonalPrism(edge_length=self.g.duct_inner_edge_cm, orientation="x")

        background = openmc.Cell(
            name=f"{name}_background_coolant",
            fill=self.m["sodium"],
            region=+duct_outer
        )

        duct = openmc.Cell(
            name=f"{name}_duct",
            fill=self.m["ht9"],
            region=-duct_outer & +duct_inner
        )

        interior = openmc.Cell(
            name=f"{name}_interior_region",
            fill=interior_mat,
            region=-duct_inner
        )

        u.add_cells([background, duct, interior])
        return u

    def _solid_assy(self, name: str, fill_mat: openmc.Material) -> openmc.Universe:
        u = openmc.Universe(name=f"assy_{name}")

        duct_outer = openmc.model.HexagonalPrism(edge_length=self.g.duct_edge_cm, orientation="x")
        duct_inner = openmc.model.HexagonalPrism(edge_length=self.g.duct_inner_edge_cm, orientation="x")

        background = openmc.Cell(
            name=f"{name}_background_coolant",
            fill=self.m["sodium"],
            region=+duct_outer
        )

        duct = openmc.Cell(
            name=f"{name}_duct",
            fill=self.m["ht9"],
            region=-duct_outer & +duct_inner
        )

        interior = openmc.Cell(
            name=f"{name}_interior_region",
            fill=fill_mat,
            region=-duct_inner
        )

        u.add_cells([background, duct, interior])
        return u
