import numpy as np

class GeometryManager:
    def __init__(self, mesh, fuel_outer_radius, cladding_inner_radius, cladding_outer_radius, height):
        self.mesh = mesh
        self.fuel_outer_radius = fuel_outer_radius
        self.cladding_inner_radius = cladding_inner_radius
        self.cladding_outer_radius = cladding_outer_radius
        self.height = height

        # To be populated in initialize_geometry()
        self.idx_fuel_outer = None
        self.idx_clad_inner = None
        self.idx_clad_outer = None
        self.idx_height = None

    def initialize_geometry(self, atol=1e-6):
        """
        Store original node indices for fuel outer, cladding inner, and outer surfaces.
        """
        mesh_geometry = self.mesh.geometry.x
        r = np.sqrt(mesh_geometry[:, 0]**2 + mesh_geometry[:, 1]**2)
        z = mesh_geometry[:, 2]

        self.idx_fuel_outer = np.where(np.isclose(r, self.fuel_outer_radius, atol=atol))[0]
        self.idx_clad_inner = np.where(np.isclose(r, self.cladding_inner_radius, atol=atol))[0]
        self.idx_clad_outer = np.where(np.isclose(r, self.cladding_outer_radius, atol=atol))[0]
        self.idx_height = np.where(np.isclose(z, self.height, atol=atol))[0]

        # print(f"Original indices stored: fuel_outer={len(self.idx_fuel_outer)}, clad_inner={len(self.idx_clad_inner)}, clad_outer={len(self.idx_clad_outer)}")

    def geometry_summary(self, displacement):
        """
        Print a detailed summary of radial and axial geometry:
        - Input (reference) radii
        - Current mesh radii (undeformed)
        - Radial displacements
        - Updated estimated radii
        - Axial height and displacement
        """

        # Extract undeformed mesh geometry
        mesh_geometry = self.mesh.geometry.x
        r_initial = np.sqrt(mesh_geometry[:, 0]**2 + mesh_geometry[:, 1]**2)
        z_initial = mesh_geometry[:, 2]

        # Displacement field
        u = displacement.x.array.reshape(mesh_geometry.shape)
        u_r = np.sqrt(u[:, 0]**2 + u[:, 1]**2)
        u_z = u[:, 2]

        # Estimated deformed radial positions
        r_updated = r_initial + u_r
        z_updated = z_initial + u_z

        # Axial info
        height = np.max(z_initial) - np.min(z_initial)
        avg_axial_disp = np.mean(u_z)
        max_axial_disp = np.max(u_z)
        min_axial_disp = np.min(u_z)

        print("\n📐 Radial & Axial Geometry Summary")
        print("────────────────────────────────────────────")

        print("📎 Input Radii (m) (from FuelRod):")
        print(f"   Fuel outer       : {self.fuel_outer_radius:.6e}")
        print(f"   Cladding inner   : {self.cladding_inner_radius:.6e}")
        print(f"   Cladding outer   : {self.cladding_outer_radius:.6e}")

        print("\n📍 Mesh Radii (m) (undeformed):")
        print(f"   Fuel outer (avg) : {np.mean(r_initial[self.idx_fuel_outer]):.6e}")
        print(f"   Cladding inner   : {np.mean(r_initial[self.idx_clad_inner]):.6e}")
        print(f"   Cladding outer   : {np.mean(r_initial[self.idx_clad_outer]):.6e}")

        print("\n🔁 Radial Displacements (m):")
        print(f"   Fuel outer (avg) : {np.mean(u_r[self.idx_fuel_outer]):.3e}")
        print(f"   Cladding inner   : {np.mean(u_r[self.idx_clad_inner]):.3e}")
        print(f"   Cladding outer   : {np.mean(u_r[self.idx_clad_outer]):.3e}")

        print("\n📏 Updated Radii (m) (estimated):")
        print(f"   Fuel outer       : {np.mean(r_updated[self.idx_fuel_outer]):.6e}")
        print(f"   Cladding inner   : {np.mean(r_updated[self.idx_clad_inner]):.6e}")
        print(f"   Cladding outer   : {np.mean(r_updated[self.idx_clad_outer]):.6e}")

        print("\n📐 Axial Geometry:")
        print(f"   Original height  : {height:.6e} m")
        print(f"   Axial disp (avg) : {avg_axial_disp:.3e} m")
        print(f"   Axial disp (min/max): {min_axial_disp:.3e} / {max_axial_disp:.3e} m")
        print(f"   Height           : {np.max(z_updated) - np.min(z_updated):.6e} m")

    def update_geometry(self, displacement):
        """
        Update effective geometry (radii and height) based on displacement field,
        without modifying the mesh coordinates.
        """

        # Reshape geometry and displacement
        u = displacement.x.array.reshape(self.mesh.geometry.x.shape)
        mesh_geometry = self.mesh.geometry.x

        # Radial: estimate r_deformed = r + u_r
        r_initial = np.sqrt(mesh_geometry[:, 0]**2 + mesh_geometry[:, 1]**2)
        u_r = np.sqrt(u[:, 0]**2 + u[:, 1]**2)
        r_deformed = r_initial + u_r

        # Update radii using max/mean in selected regions
        self.fuel_outer_radius = np.mean(r_deformed[self.idx_fuel_outer])
        self.cladding_inner_radius = np.mean(r_deformed[self.idx_clad_inner])
        self.cladding_outer_radius = np.mean(r_deformed[self.idx_clad_outer])

        # Axial: estimate z_deformed = z + u_z
        z_initial = mesh_geometry[:, 2]
        z_deformed = z_initial + u[:, 2]
        self.height = np.max(z_deformed) - np.min(z_deformed)

        print("\n✅ Updated geometric values (from displacement):")
        print(f"   Fuel outer radius       = {self.fuel_outer_radius:.6e} m")
        print(f"   Cladding inner radius   = {self.cladding_inner_radius:.6e} m")
        print(f"   Cladding outer radius   = {self.cladding_outer_radius:.6e} m")
        print(f"   Axial height            = {self.height:.6e} m")
        