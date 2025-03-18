import numpy as np

class GeometryManager:
    def __init__(self, mesh, fuel_outer_radius, cladding_inner_radius, cladding_outer_radius):
        self.mesh = mesh
        self.fuel_outer_radius = fuel_outer_radius
        self.cladding_inner_radius = cladding_inner_radius
        self.cladding_outer_radius = cladding_outer_radius

        # To be populated in initialize_geometry()
        self.idx_fuel_outer = None
        self.idx_clad_inner = None
        self.idx_clad_outer = None

    def initialize_geometry(self, atol=1e-6):
        """
        Store original node indices for fuel outer, cladding inner, and outer surfaces.
        """
        mesh_geometry = self.mesh.geometry.x
        r = np.sqrt(mesh_geometry[:, 0]**2 + mesh_geometry[:, 1]**2)

        self.idx_fuel_outer = np.where(np.isclose(r, self.fuel_outer_radius, atol=atol))[0]
        self.idx_clad_inner = np.where(np.isclose(r, self.cladding_inner_radius, atol=atol))[0]
        self.idx_clad_outer = np.where(np.isclose(r, self.cladding_outer_radius, atol=atol))[0]

        print(f"Original indices stored: fuel_outer={len(self.idx_fuel_outer)}, \
              clad_inner={len(self.idx_clad_inner)}, clad_outer={len(self.idx_clad_outer)}")

    def apply_displacement(self, displacement):
        """
        Apply displacement field to the mesh geometry.
        """
        u = displacement.x.array.reshape(self.mesh.geometry.x.shape)
        self.mesh.geometry.x[:] += u

    def update_deformed_dimensions(self):
        """
        Recalculate and update fuel/cladding radii from deformed geometry.
        """
        mesh_geometry = self.mesh.geometry.x
        r_new = np.sqrt(mesh_geometry[:, 0]**2 + mesh_geometry[:, 1]**2)

        self.fuel_outer_radius = np.max(r_new[self.idx_fuel_outer])
        self.cladding_inner_radius = np.min(r_new[self.idx_clad_inner])
        self.cladding_outer_radius = np.max(r_new[self.idx_clad_outer])

        print("\n✅ Updated Dimensions after Deformation:")
        print(f"Fuel Outer Radius (deformed): {self.fuel_outer_radius:.6e} m")
        print(f"Cladding Inner Radius (deformed): {self.cladding_inner_radius:.6e} m")
        print(f"Cladding Outer Radius (deformed): {self.cladding_outer_radius:.6e} m")

    def radial_displacement_summary(self, displacement):
        u = displacement.x.array.reshape(self.mesh.geometry.x.shape)
        disp_r = np.sqrt(u[:, 0]**2 + u[:, 1]**2)

        mesh_geometry = self.mesh.geometry.x
        r_new = np.sqrt(mesh_geometry[:, 0]**2 + mesh_geometry[:, 1]**2)

        print("\nRadial Displacement Summary:")
        print(f"Max radial displacement: {np.max(disp_r):.3e} m")
        print(f"Avg at fuel outer: {np.mean(disp_r[self.idx_fuel_outer]):.3e} m")
        print(f"Avg at cladding inner: {np.mean(disp_r[self.idx_clad_inner]):.3e} m")
        print(f"Avg at cladding outer: {np.mean(disp_r[self.idx_clad_outer]):.3e} m")

        print(f"Radial displacement at fuel outer radius (max): {np.max(disp_r[self.idx_fuel_outer]):.3e} m")
        print(f"Radial displacement at cladding inner radius (max): {np.max(disp_r[self.idx_clad_inner]):.3e} m")
        print(f"Radial displacement at cladding outer radius (max): {np.max(disp_r[self.idx_clad_outer]):.3e} m")

        print("\nSpread of deformed radial positions:")
        print(f"Fuel outer radius (min/avg/max): {np.min(r_new[self.idx_fuel_outer]):.6e} / {np.mean(r_new[self.idx_fuel_outer]):.6e} / {np.max(r_new[self.idx_fuel_outer]):.6e}")
        print(f"Cladding inner radius (min/avg/max): {np.min(r_new[self.idx_clad_inner]):.6e} / {np.mean(r_new[self.idx_clad_inner]):.6e} / {np.max(r_new[self.idx_clad_inner]):.6e}")
        print(f"Cladding outer radius (min/avg/max): {np.min(r_new[self.idx_clad_outer]):.6e} / {np.mean(r_new[self.idx_clad_outer]):.6e} / {np.max(r_new[self.idx_clad_outer]):.6e}")