# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.1.0 (2025)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---


import dolfinx
import numpy as np
from petsc4py import PETSc


class GapModel:
    def __init__(self):
        print("[GapModel] initializer")

    def set_gap_conductance(self, T_i):

        # Gap condutance
        if self.gap_model == "Fixed":
            h_gap_value = self.h_gap_value
            print(f"\n[INFO] Applying fixed gap conductance (h_gap = {h_gap_value:.2f} W/m²K)")

        elif self.gap_model == "Gas":

            self.set_gap_temperature(T_i)

            gas_thermal_conductivity = (
                self.h_gap_value * 1e-4 * self.gap_temperature**0.79
            )  # (W/m-K)

            gap_size = self.average_gap_distance(
                self.mesh,
                self.facet_tags,
                label_a=self.label_map["lateral_1"],
                label_b=self.label_map["inner_2"],
            )

            h_gap_value = gas_thermal_conductivity / gap_size  # (W/m2-K)

            print(f"  → k_gas       = {gas_thermal_conductivity:.2f} W/m·K")
            print(f"  → gap_size    = {gap_size*1e3:.3f} mm")
            print(f"  → h_gap       = {h_gap_value:.2f} W/m²K")

        else:
            h_gap_value = 0.0

        # Contact-coupled conductance (Todreas & Kazimi, Nuclear Systems I,
        # 3rd ed., Eq. 8.141/8.142): on gap closure the pellet-clad contact
        # pressure adds a contact term to the open-gap value.
        if self.on.get("contact", False) and getattr(self, "gap_contact_coupling", False):
            h_contact = self.contact_conductance()
            if h_contact > 0.0:
                print(
                    f"  → h_open      = {h_gap_value:.2f} W/m²K, "
                    f"h_contact = {h_contact:.2f} W/m²K (gap closed)"
                )
                h_gap_value += h_contact

        h_gap = dolfinx.fem.Constant(self.mesh, PETSc.ScalarType(h_gap_value))

        return h_gap

    def contact_conductance(self):
        """
        Solid-contact gap conductance on gap closure, Todreas & Kazimi,
        *Nuclear Systems I*, 3rd ed., Eq. 8.141 (Ross-Stoute form):

            h_contact = C * (2 k_f k_c)/(k_f + k_c) * P_i / (H * sqrt(delta_g)),

        with the pellet-clad contact pressure P_i taken from the penalty
        :class:`ContactModel`. The empirical constant C = 10 ft^-1/2 is here
        expressed in SI as C_SI = 18.11 m^-1/2 so that, with k in W/(m.K),
        delta_g in m and the dimensionless ratio P_i/H, the result is W/(m^2.K).
        Returns 0 when the surfaces are not in contact (P_i = 0).
        """
        P_i = float(getattr(self, "_last_pressure", 0.0))  # Pa, from ContactModel
        if P_i <= 0.0:
            return 0.0

        # harmonic mean of the two solid conductivities
        ks = [
            float(m["k"])
            for m in self.materials.values()
            if isinstance(m.get("k"), (int, float))
        ]
        if len(ks) < 2:
            return 0.0
        k_f, k_c = ks[0], ks[1]
        k_harm = 2.0 * k_f * k_c / (k_f + k_c)

        H = self.gap_meyer_hardness            # Pa (Meyer hardness, softer solid)
        delta_g = self.gap_contact_thickness   # m (roughness-based gas space)
        C_SI = 18.11                           # m^-1/2  (= 10 ft^-1/2, Eq. 8.141)

        return C_SI * k_harm * (P_i / H) / (delta_g ** 0.5)

    def set_gap_temperature(self, T_i):

        for label in self.materials:
            for bc_info in self.robin_thermal[label]:
                region_id = bc_info["id"]
                pair_region = bc_info["pair"]

                dofs_here = self.mgr.locate_facets_dofs(region_id, self.V_t)
                dofs_other = self.mgr.locate_facets_dofs(self.label_map[pair_region], self.V_t)

                T_here = T_i.x.array[dofs_here]
                T_other = T_i.x.array[dofs_other]

                self.gap_temperature = 0.5 * (T_here.mean() + T_other.mean())

                print(
                    f"  → Average gap temperature between {label} and {pair_region}: {self.gap_temperature:.2f} K"
                )

                break
            break

    def average_gap_distance(self, mesh, ft, label_a, label_b):

        facets_a = ft.find(label_a)  # e.g., facets of cyl_1_interface
        facets_b = ft.find(label_b)  # e.g., facets of cyl_2_interface

        x = mesh.geometry.x
        topology = mesh.topology
        fdim = topology.dim - 1

        topology.create_connectivity(fdim, 0)  # facet → vertex
        connectivity = topology.connectivity(fdim, 0)

        def facet_centroids(facets):
            centroids = []
            for facet in facets:
                vertex_indices = connectivity.links(facet)
                coords = x[vertex_indices]
                centroids.append(coords.mean(axis=0))
            return np.array(centroids)

        centroids_a = facet_centroids(facets_a)
        centroids_b = facet_centroids(facets_b)

        # min_len = min(len(centroids_a), len(centroids_b))
        # distances = np.linalg.norm(centroids_a[:min_len] - centroids_b[:min_len], axis=1)

        from scipy.spatial import cKDTree

        tree_b = cKDTree(centroids_b)
        distances, _ = tree_b.query(centroids_a, k=1)

        # print(f"Average distance cyl_1 - cyl_2: {distances.mean()*1e3:.3e} mm")

        return distances.mean()
