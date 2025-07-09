import sys, numpy as np
import dolfinx, ufl
from petsc4py import PETSc

class GapModel:
    def __init__(self):
        print("__GapModel initializer__")

    def set_gap_conductance(self):

        # Gap condutance
        if self.gap_model == "Fixed":
            h_gap_value = self.h_gap_value
            print(f"\n[INFO] Applying fixed gap conductance (h_gap = {h_gap_value:.2f} W/m²K)")

        elif self.gap_model == "Gas":
            gas_thermal_conductivity = self.h_gap_value * 1e-4 * self.gap_temperature**0.79 # (W/m-K)

            gap_size = self.average_gap_distance(self.mesh, self.ft,
                                    label_a=self.label_map["fuel_interface"],
                                    label_b=self.label_map["clad_interface"])

            h_gap_value = gas_thermal_conductivity / gap_size # (W/m2-K)
            print(f"\n[INFO] Applying gas gap conductance (h_gap = {h_gap_value:.2f} W/m²K)")

        else:
            h_gap_value = 0.0

        h_gap = dolfinx.fem.Constant(self.mesh, PETSc.ScalarType(h_gap_value))

        return h_gap
    
    def set_gap_temperature(self):
        self.gap_temperature = 10
    
    def average_gap_distance(self, mesh, ft, label_a, label_b):

        facets_a = ft.find(label_a) # e.g., facets of fuel_interface
        facets_b = ft.find(label_b) # e.g., facets of clad_interface

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

        min_len = min(len(centroids_a), len(centroids_b))
        distances = np.linalg.norm(centroids_a[:min_len] - centroids_b[:min_len], axis=1)

        print(f"Average distance fuel-cladding: {distances.mean()*1e3:.3e} mm")

        return distances.mean()