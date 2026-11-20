# geometry/universe_factory.py

import openmc
from geometry.params import AssemblyGeom

class UniverseFactory:
    def __init__(self, geom: AssemblyGeom, mats: dict):
        self.g = geom
        self.m = mats

    def build_all(self):
        # Build shared universe objects so aliases point to the same instance
        _cr  = self._solid_assy("control_rod",   self.m["b4c"])
        _us  = self._solid_assy("ult_shutdown",  self.m["b4c_ult"])
        _gem = self._solid_assy("gem",           self.m["sodium"])
        # Fuel universes only for materials the loadout defined (e.g. no-blanket cores)
        fuel = {
            univ: self._fuel_assy(name, self.m[mat_key])
            for univ, name, mat_key in [
                ("u_inner_fuel",       "inner_fuel",       "fuel_driver_inner_homog"),
                ("u_outer_fuel",       "outer_fuel",       "fuel_driver_outer_homog"),
                ("u_internal_blanket", "internal_blanket", "fuel_blanket_internal_homog"),
                ("u_radial_blanket",   "radial_blanket",   "fuel_blanket_radial_homog"),
            ]
            if mat_key in self.m
        }
        return {
            **fuel,
            "u_reflector":        self._solid_assy("reflector", self.m["reflector"]),
            "u_shield":            self._solid_assy("shield",    self.m["shield"]),
            # control rods — two names, one universe
            "u_control":           _cr,
            "u_control_rod":       _cr,
            # ultimate shutdown rods — two names, one universe
            "u_ult_shutdown":      _us,
            "u_ultimate_shutdown": _us,
            # GEMs (sodium-flooded tube, normal operating condition)
            "u_gem":               _gem,
            "u_source":            self._solid_assy("source", self.m["sodium"]),
            # Control / shutdown position with the absorber withdrawn above the
            # core: a sodium-filled duct over the active height (rods-out state)
            "u_control_withdrawn": self._solid_assy("control_withdrawn", self.m["sodium"]),
            # Empty corner notch positions — no duct, just pool sodium
            "u_empty":             self._sodium_assy("empty"),
        }

    def build_for_material(
        self, position_name: str, interior_mat: openmc.Material
    ) -> openmc.Universe:
        """
        Build a fuel assembly universe for a specific position and material.
        Used by cycle_runner for per-assembly Tier-2 model construction.
        position_name should be the assembly's material_id for unique naming.
        """
        return self._fuel_assy(position_name, interior_mat)

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

    def _sodium_assy(self, name: str) -> openmc.Universe:
        """No-duct sodium universe for empty (absent) lattice positions."""
        u = openmc.Universe(name=f"assy_{name}")
        u.add_cell(openmc.Cell(name=f"{name}_coolant", fill=self.m["sodium"]))
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
