# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.2.0 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

import dolfinx
import numpy as np
import ufl
from dolfinx import fem, mesh

from z3st.core.diagnostic import log


class MeshManager:
    """Handles Dolfinx mesh topology, tagging, and geometry utilities."""

    def __init__(
        self,
        mesh_obj: mesh.Mesh,
        cell_tags: mesh.MeshTags,
        facet_tags: mesh.MeshTags,
        geometry: dict = None,
    ):
        self.mesh = mesh_obj
        self.cell_tags = cell_tags
        self.facet_tags = facet_tags
        self.geometry = geometry or {}

        # --. Basic topology --..
        self.tdim = self.mesh.topology.dim
        self.fdim = self.tdim - 1
        log.info(f"Mesh topology dimension d={self.tdim}")

        # Ensure connectivities exist (needed for BCs, dof search, etc.)
        self.mesh.topology.create_connectivity(self.fdim, self.tdim)
        self.mesh.topology.create_connectivity(self.tdim, self.tdim)

        # --.. Boundary facets --..
        self.boundary_facets = dolfinx.mesh.exterior_facet_indices(self.mesh.topology)
        log.debug(f"Boundary facets: {len(self.boundary_facets)}")

        # --. Volume tags --..
        log.info("\nAvailable volume tags (dx):")
        tag_values = self.cell_tags.values
        unique_tags = sorted(set(tag_values))
        for tag in unique_tags:
            log.info(f"  Tag ID: {tag}")

        # --. Facet tags --..
        if self.facet_tags is not None:
            unique_facets = np.unique(self.facet_tags.values)
            log.info(f"\nUnique tags found in facet data: {unique_facets}")
        else:
            log.warning("No facet tags found in mesh.")

        # --. Label map --..
        self.label_map = self.geometry.get("labels", {})
        if self.label_map:
            log.info(f"Label map loaded from geometry:")
            for label, tag in sorted(self.label_map.items(), key=lambda kv: kv[1]):
                log.info(f"  {label:<12} → {tag}")
        else:
            log.warning("No label map found in geometry; defaulting to empty dict.")
            self.label_map = {}

        # --. Geometry attributes --..
        self.geometry_type = self.geometry.get("geometry_type", "").lower()
        self.normal = ufl.FacetNormal(self.mesh)

        self._init_geometry_parameters()

    @staticmethod
    def _first_present(g: dict, keys: list, default=None):
        """Return g[k] for the first key in keys that g contains, else default."""
        for k in keys:
            if k in g:
                return g[k]
        return default

    def _init_geometry_parameters(self):
        """Compute useful geometric quantities (like area, perimeter, etc.)."""
        g = self.geometry
        self.Lz = float(g.get("Lz", 0.0))
        log.info(f"  Lz = {self.Lz:.3f} m")

        if self.geometry_type == "rect":
            Lx_keys = ["Lx", "length_x"]
            Ly_keys = ["Ly", "length_y"]
            self.Lx = float(self._first_present(g, Lx_keys, None))
            self.Ly = float(self._first_present(g, Ly_keys, None))
            self.perimeter = (self.Lx + self.Ly) * 2.0
            self.area = self.Lx * self.Ly
            log.info(f"  Lx = {self.Lx:.3f} m, Ly = {self.Ly:.3f} m")

        elif self.geometry_type in ["cyl", "cylinder"]:
            inner_radius_keys = ["inner_radius", "inner_radius_1", "Ri"]
            outer_radius_keys = ["outer_radius", "outer_radius_1", "Ro"]
            self.inner_radius = self._first_present(g, inner_radius_keys, 0.0)
            self.outer_radius = self._first_present(g, outer_radius_keys, None)
            self.perimeter = 2.0 * np.pi * self.outer_radius
            self.area = np.pi * (self.outer_radius**2 - self.inner_radius**2)
            log.info(f"  Ri = {self.inner_radius:.3e} m, Ro = {self.outer_radius:.3e} m")

        elif self.geometry_type == "cyl-cyl":
            inner_1_keys = ["inner_radius_1", "Ri_1"]
            outer_1_keys = ["outer_radius_1", "Ro_1"]
            inner_2_keys = ["inner_radius_2", "Ri_2"]
            outer_2_keys = ["outer_radius_2", "Ro_2"]
            self.inner_radius_1 = self._first_present(g, inner_1_keys, None)
            self.outer_radius_1 = self._first_present(g, outer_1_keys, None)
            self.inner_radius_2 = self._first_present(g, inner_2_keys, None)
            self.outer_radius_2 = self._first_present(g, outer_2_keys, None)
            self.perimeter = 2.0 * np.pi * self.outer_radius_1
            self.area = np.pi * (self.outer_radius_1**2 - self.inner_radius_1**2)
            log.info(
                f"  inner_radius_1 = {self.inner_radius_1:.2e} m, outer_radius_1 = {self.outer_radius_1:.2e} m"
            )
            log.info(
                f"  inner_radius_2 = {self.inner_radius_2:.2e} m, outer_radius_2 = {self.outer_radius_2:.2e} m"
            )

        elif self.geometry_type == "sphere":
            self.inner_radius = g.get("Ri", None)
            self.outer_radius = g.get("Ro", None)
            self.perimeter = 2.0 * np.pi * self.outer_radius
            self.area = np.pi * (self.outer_radius**2 - self.inner_radius**2)
            log.info(f"  Ri = {self.inner_radius:.2e} m, Ro = {self.outer_radius:.2e} m")

        else:
            self.area = float(g.get("area", 0.0))
            self.perimeter = float(g.get("perimeter", 0.0))

        log.info(f"  area = {self.area:.3e} m², perimeter = {self.perimeter:.3e} m")

    def locate_facets_dofs(self, label: int, V: fem.FunctionSpace):
        """Locate DOFs on facets by label."""
        facets = self.facet_tags.find(label)
        return fem.locate_dofs_topological(V, self.fdim, facets)

    def locate_domain_dofs(self, label: int, V: fem.FunctionSpace):
        """Locate DOFs in domain by label."""
        cells = self.cell_tags.find(label)
        return fem.locate_dofs_topological(V, self.tdim, cells)

    def summary(self):
        log.info("=== Mesh summary ===")
        log.info(f"  Topology dim: {self.tdim}")
        log.info(f"  Facet dim: {self.fdim}")
        log.info(f"  Num cells: {self.mesh.topology.index_map(self.tdim).size_global}")
        log.info(f"  Cell tags: {set(self.cell_tags.values)}")
        log.info(f"  Facet tags: {set(self.facet_tags.values)}")
        log.info(f"  Geometry type: {self.geometry.get('geometry_type', 'Not specified')}")
