# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.1.0 (2025)
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

        # --. Geometry attributes ----
        self.geometry_type = self.geometry.get("geometry_type", "").lower()
        self.normal = ufl.FacetNormal(self.mesh)
        self.subdomain_areas = {}

        self._init_geometry_parameters()

    def _init_geometry_parameters(self):
        """Compute derived geometry quantities (area, perimeter, etc.)."""
        g = self.geometry
        self.Lz = float(g.get("Lz", 0.0))
        log.info(f"  Lz = {self.Lz:.3f} m")
        
        def safe_float(key, default=0.0):
            val = g.get(key, default)
            return float(val) if val is not None else default

        if self.geometry_type == "rect":
            self.Lx = safe_float("Lx")
            self.Ly = safe_float("Ly")
            self.perimeter = (self.Lx + self.Ly) * 2.0
            self.area = self.Lx * self.Ly
            log.info(f"  Lx = {self.Lx:.3f} m, Ly = {self.Ly:.3f} m")

        elif self.geometry_type in ["cyl", "cylinder"]:
            # Check if it's a contact case with two rectangular domains (Cartesian 2D)
            if "R_mid" in g and "Ly_rect" in g:
                r_int = safe_float("R_int")
                r_mid = safe_float("R_mid")
                r_o = safe_float("R_o")
                ly = safe_float("Ly_rect")
                gap = safe_float("gap")
                
                # steel_in: from R_int to R_mid - gap/2
                area_in = (r_mid - gap/2.0 - r_int) * ly
                # steel_o: from R_mid + gap/2 to R_o
                area_o = (r_o - (r_mid + gap/2.0)) * ly
                
                self.subdomain_areas = {
                    "steel_in": area_in,
                    "steel_o": area_o
                }
                self.area = area_in + area_o
                self.perimeter = ly # Use height as reference for heat flux
                log.info(f"  Contact areas (2D): steel_in={area_in:.3e} m², steel_o={area_o:.3e} m²")
            else:
                self.inner_radius = safe_float("Ri", 0.0)
                self.outer_radius = safe_float("Ro", 0.0)
                self.perimeter = 2.0 * np.pi * self.outer_radius
                self.area = np.pi * (self.outer_radius**2 - self.inner_radius**2)
                log.info(f"  Ri = {self.inner_radius:.3e} m, Ro = {self.outer_radius:.3e} m")

        elif self.geometry_type == "cyl-cyl":
            self.inner_radius_1 = safe_float("inner_radius_1")
            self.outer_radius_1 = safe_float("outer_radius_1")
            self.inner_radius_2 = safe_float("inner_radius_2")
            self.outer_radius_2 = safe_float("outer_radius_2")
            self.perimeter = 2.0 * np.pi * self.outer_radius_1
            self.area = np.pi * (self.outer_radius_1**2 - self.inner_radius_1**2)
            log.info(f"  inner_radius_1 = {self.inner_radius_1:.2e} m, outer_radius_1 = {self.outer_radius_1:.2e} m")
            log.info(f"  inner_radius_2 = {self.inner_radius_2:.2e} m, outer_radius_2 = {self.outer_radius_2:.2e} m")

        elif self.geometry_type == "sphere":
            self.inner_radius = safe_float("Ri")
            self.outer_radius = safe_float("Ro")
            self.perimeter = 2.0 * np.pi * self.outer_radius
            self.area = np.pi * (self.outer_radius**2 - self.inner_radius**2)
            log.info(f"  Ri = {self.inner_radius:.2e} m, Ro = {self.outer_radius:.2e} m")

        else:
            self.area = float(g.get("area", 0.0))
            self.perimeter = float(g.get("perimeter", 0.0))

        if self.area <= 0:
            log.warning(f"Computed area is {self.area}. Check your geometry_type and dimensions in geometry.yaml.")
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
    
    def get_contact_facets(self, label_map):
        """
        Récupère les facettes de contact entre deux sous-domaines.
        Args:
            label_map: Dictionnaire des labels (ex: {"steel_in_right": 4, "steel_o_left": 8}).
        Returns:
            Tuple de (facettes_steel_in, facettes_steel_o).
        """
        # Récupérer les tags pour les surfaces de contact
        steel_in_right_tag = label_map.get("steel_in_right")
        steel_o_left_tag = label_map.get("steel_o_left")

        if steel_in_right_tag is None or steel_o_left_tag is None:
            raise ValueError("Tags pour steel_in_right ou steel_o_left non trouvés dans label_map.")

        # Trouver les facettes correspondantes
        steel_in_facets = self.facet_tags.find(steel_in_right_tag)
        steel_o_facets = self.facet_tags.find(steel_o_left_tag)

        return steel_in_facets, steel_o_facets