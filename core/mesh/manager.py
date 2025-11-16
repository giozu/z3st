# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.1.0 (2025)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

import numpy as np
import dolfinx
import ufl
from dolfinx import fem, mesh
from mpi4py import MPI
from core.diagnostic import log

class MeshManager:
    """Handles Dolfinx mesh topology, tagging, and geometry utilities."""

    def __init__(self, mesh_obj: mesh.Mesh, cell_tags: mesh.MeshTags, facet_tags: mesh.MeshTags, geometry: dict = None):
        self.mesh = mesh_obj
        self.cell_tags = cell_tags
        self.facet_tags = facet_tags
        self.geometry = geometry or {}

        # --- Basic topology ---
        self.tdim = self.mesh.topology.dim
        self.fdim = self.tdim - 1
        log.info(f"Mesh topology dimension d={self.tdim}")

        # Ensure connectivities exist (needed for BCs, dof search, etc.)
        self.mesh.topology.create_connectivity(self.fdim, self.tdim)
        self.mesh.topology.create_connectivity(self.tdim, self.tdim)

        # --- Boundary facets ---
        self.boundary_facets = dolfinx.mesh.exterior_facet_indices(self.mesh.topology)
        log.debug(f"Boundary facets: {len(self.boundary_facets)}")

        # --- Volume tags ---
        log.info("\nAvailable volume tags (dx):")
        tag_values = self.cell_tags.values
        unique_tags = sorted(set(tag_values))
        for tag in unique_tags:
            log.info(f"  Tag ID: {tag}")

        # --- Facet tags ---
        unique_facets = np.unique(self.facet_tags.values)
        log.info(f"\nUnique tags found in facet data: {unique_facets}")

        # --- Label map ---
        self.label_map = self.geometry.get("labels", {})
        if self.label_map:
            log.info(f"Label map loaded from geometry:")
            for label, tag in sorted(self.label_map.items(), key=lambda kv: kv[1]):
                log.info(f"  {label:<12} → {tag}")
        else:
            log.warning("No label map found in geometry; defaulting to empty dict.")
            self.label_map = {}

        # --- Geometry attributes ---
        self.geometry_type = self.geometry.get("geometry_type", "").lower()
        self.normal = ufl.FacetNormal(self.mesh)

        self._init_geometry_parameters()

    def _init_geometry_parameters(self):
        """Compute derived geometry quantities (area, perimeter, etc.)."""
        g = self.geometry
        self.Lz = float(g.get("Lz", 0.0))
        log.info(f"  Lz = {self.Lz:.3f} m")

        if self.geometry_type == "rect":
            self.Lx = float(g.get("Lx", None))
            self.Ly = float(g.get("Ly", None))
            self.perimeter = (self.Lx + self.Ly) * 2.0
            self.area = self.Lx * self.Ly
            log.info(f"  Lx = {self.Lx:.3f} m, Ly = {self.Ly:.3f} m")

        elif self.geometry_type in ["cyl", "cylinder"]:
            self.inner_radius = g.get("Ri", 0.0)
            self.outer_radius = g.get("Ro", None)
            self.perimeter = 2.0 * np.pi * self.outer_radius
            self.area = np.pi * (self.outer_radius**2 - self.inner_radius**2)
            log.info(f"  Ri = {self.inner_radius:.3e} m, Ro = {self.outer_radius:.3e} m")

        elif self.geometry_type == "cyl-cyl":
            self.inner_radius_1 = g.get("inner_radius_1", None)
            self.outer_radius_1 = g.get("outer_radius_1", None)
            self.inner_radius_2 = g.get("inner_radius_2", None)
            self.outer_radius_2 = g.get("outer_radius_2", None)
            self.perimeter = 2. * np.pi * self.outer_radius_1
            self.area = np.pi * (self.outer_radius_1**2 - self.inner_radius_1**2)
            log.info(f"  inner_radius_1 = {self.inner_radius_1:.2e} m, outer_radius_1 = {self.outer_radius_1:.2e} m")
            log.info(f"  inner_radius_2 = {self.inner_radius_2:.2e} m, outer_radius_2 = {self.outer_radius_2:.2e} m")

        elif self.geometry_type == "sphere":
            self.inner_radius = g.get("Ri", None)
            self.outer_radius = g.get("Ro", None)
            self.perimeter = 2. * np.pi * self.outer_radius
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
        log.info(f"  Geometry type: {self.geometry_type}")
