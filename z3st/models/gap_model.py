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

        h_gap = dolfinx.fem.Constant(self.mesh, PETSc.ScalarType(h_gap_value))

        return h_gap

    def set_gap_temperature(self, T_i):

        for label in self.materials:
            for bc_info in self.robin_thermal[label]:
                region_id = bc_info["id"]
                pair_region = bc_info["pair"]

                dofs_here = self.locateFacetDofs(region_id, self.V_t)
                dofs_other = self.locateFacetDofs(self.label_map[pair_region], self.V_t)

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
