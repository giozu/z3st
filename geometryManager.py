import numpy as np

class GeometryManager:
    def __init__(self, mesh, surface_tags,
             fuel_outer_label, cladding_inner_label, cladding_outer_label,
             cladding_bottom_label, cladding_top_label,
             fuel_outer_radius, cladding_inner_radius, cladding_outer_radius, height,
             dofs_fuel_outer, dofs_cladding_inner, dofs_cladding_outer,
             dofs_cladding_bottom, dofs_cladding_top):
        
        self.mesh = mesh
        self.surface_tags = surface_tags

        self.fuel_outer_label = fuel_outer_label
        self.cladding_inner_label = cladding_inner_label
        self.cladding_outer_label = cladding_outer_label
        self.cladding_bottom_label = cladding_bottom_label
        self.cladding_top_label = cladding_top_label

        self.fuel_outer_radius = fuel_outer_radius
        self.cladding_inner_radius = cladding_inner_radius
        self.cladding_outer_radius = cladding_outer_radius
        self.height = height

        self.dofs_fuel_outer = dofs_fuel_outer
        self.dofs_cladding_inner = dofs_cladding_inner
        self.dofs_cladding_outer = dofs_cladding_outer
        self.dofs_cladding_bottom = dofs_cladding_bottom
        self.dofs_cladding_top = dofs_cladding_top

        print("📏 Initial geometry loaded:")
        print(f"   Fuel outer radius     = {self.fuel_outer_radius*1e6:.3f} μm")
        print(f"   Cladding inner radius = {self.cladding_inner_radius*1e6:.3f} μm")
        print(f"   Cladding outer radius = {self.cladding_outer_radius*1e6:.3f} μm")
        print(f"   Height                = {self.height*1e6:.3f} μm")

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

    def geometry_summary(self, u_fuel, u_cladding):
        """
        Print deformed geometry summary using u_fuel and u_cladding.
        """
        print("\n📐 Updated geometric values (from displacement):")

        # Coordinate mesh
        x = self.mesh.geometry.x[:, 0]
        y = self.mesh.geometry.x[:, 1]
        z = self.mesh.geometry.x[:, 2]
        r = np.sqrt(x**2 + y**2)

        # Displacement arrays
        u_fuel_array = u_fuel.x.array.reshape(self.mesh.geometry.x.shape)
        u_clad_array = u_cladding.x.array.reshape(self.mesh.geometry.x.shape)

        # Radial displacement
        u_r_fuel = np.sqrt(u_fuel_array[:, 0]**2 + u_fuel_array[:, 1]**2)
        u_r_clad = np.sqrt(u_clad_array[:, 0]**2 + u_clad_array[:, 1]**2)

        # Axial displacement
        u_z_fuel = u_fuel_array[:, 2]
        u_z_clad = u_clad_array[:, 2]

        # Fuel outer surface
        r_fuel_outer = r[self.dofs_fuel_outer] + u_r_fuel[self.dofs_fuel_outer]
        r_fuel_outer_avg = np.mean(r_fuel_outer)

        # Cladding inner and outer surfaces
        r_clad_inner = r[self.dofs_cladding_inner] + u_r_clad[self.dofs_cladding_inner]
        r_clad_outer = r[self.dofs_cladding_outer] + u_r_clad[self.dofs_cladding_outer]

        r_clad_inner_avg = np.mean(r_clad_inner)
        r_clad_outer_avg = np.mean(r_clad_outer)

        # Apparent gap
        gap_apparent = (r_clad_inner_avg - r_fuel_outer_avg) * 1e6  # μm

        # Axial deformation (top-bottom)
        z_bot = z[self.dofs_cladding_bottom] + u_z_clad[self.dofs_cladding_bottom]
        z_top = z[self.dofs_cladding_top] + u_z_clad[self.dofs_cladding_top]

        dz = np.mean(z_top) - np.mean(z_bot)

        print(f"   Fuel outer radius       = {r_fuel_outer_avg*1e6:.3f} μm")
        print(f"   Cladding inner radius   = {r_clad_inner_avg*1e6:.3f} μm")
        print(f"   Cladding outer radius   = {r_clad_outer_avg*1e6:.3f} μm")
        print(f"   Axial height            = {dz*1e6:.3f} μm")
        print(f"   Apparent gap thickness  = {gap_apparent:.3f} μm")

    def update_geometry(self, u_fuel, u_cladding):
        """
        Update geometry based on current displacements.
        """
        # Coordinate mesh
        x = self.mesh.geometry.x[:, 0]
        y = self.mesh.geometry.x[:, 1]
        z = self.mesh.geometry.x[:, 2]
        r = np.sqrt(x**2 + y**2)

        # Displacement arrays
        u_fuel_array = u_fuel.x.array.reshape(self.mesh.geometry.x.shape)
        u_clad_array = u_cladding.x.array.reshape(self.mesh.geometry.x.shape)

        u_r_fuel = np.sqrt(u_fuel_array[:, 0]**2 + u_fuel_array[:, 1]**2)
        u_r_clad = np.sqrt(u_clad_array[:, 0]**2 + u_clad_array[:, 1]**2)
        u_z_clad = u_clad_array[:, 2]

        # Fuel outer radius
        r_fuel_outer = r[self.dofs_fuel_outer] + u_r_fuel[self.dofs_fuel_outer]
        self.fuel_outer_radius = np.mean(r_fuel_outer)

        # Cladding inner/outer radius
        r_clad_inner = r[self.dofs_cladding_inner] + u_r_clad[self.dofs_cladding_inner]
        r_clad_outer = r[self.dofs_cladding_outer] + u_r_clad[self.dofs_cladding_outer]

        self.cladding_inner_radius = np.mean(r_clad_inner)
        self.cladding_outer_radius = np.mean(r_clad_outer)

        # Axial deformation (top-bottom)
        z_bot = z[self.dofs_cladding_bottom] + u_z_clad[self.dofs_cladding_bottom]
        z_top = z[self.dofs_cladding_top] + u_z_clad[self.dofs_cladding_top]

        dz = np.mean(z_top) - np.mean(z_bot)
        self.height = dz

    def compute_deformed_geometry(self, u_fuel, u_cladding):
        """
        Compute the deformed fuel outer radius, cladding inner radius, and cladding outer radius,
        and the apparent fuel-cladding gap (mean, min, max).

        Parameters:
        ----------
        u_fuel : dolfinx.fem.Function
            Fuel displacement field
        u_cladding : dolfinx.fem.Function
            Cladding displacement field

        Returns:
        -------
        gap_mean, gap_min, gap_max,
        fuel_outer_mean, fuel_outer_min, fuel_outer_max,
        clad_inner_mean, clad_inner_min, clad_inner_max,
        clad_outer_mean, clad_outer_min, clad_outer_max
        """

        # Coordinates
        x = self.mesh.geometry.x[:, 0]
        y = self.mesh.geometry.x[:, 1]
        r = np.sqrt(x**2 + y**2)

        # Displacements
        u_fuel_array = u_fuel.x.array.reshape(self.mesh.geometry.x.shape)
        u_clad_array = u_cladding.x.array.reshape(self.mesh.geometry.x.shape)

        u_r_fuel = np.sqrt(u_fuel_array[:, 0]**2 + u_fuel_array[:, 1]**2)
        u_r_clad = np.sqrt(u_clad_array[:, 0]**2 + u_clad_array[:, 1]**2)

        # Masks for surfaces
        mask_fuel_outer = self.dofs_fuel_outer
        mask_clad_inner = self.dofs_cladding_inner
        mask_clad_outer = self.dofs_cladding_outer

        # Deformed radii
        r_fuel_outer = r[mask_fuel_outer] + u_r_fuel[mask_fuel_outer]
        r_clad_inner = r[mask_clad_inner] + u_r_clad[mask_clad_inner]
        r_clad_outer = r[mask_clad_outer] + u_r_clad[mask_clad_outer]

        # Compute apparent gap
        gap = r_clad_inner - r_fuel_outer

        # Statistics
        gap_mean = np.mean(gap)
        gap_min = np.min(gap)
        gap_max = np.max(gap)

        fuel_outer_mean = np.mean(r_fuel_outer)
        fuel_outer_min = np.min(r_fuel_outer)
        fuel_outer_max = np.max(r_fuel_outer)

        clad_inner_mean = np.mean(r_clad_inner)
        clad_inner_min = np.min(r_clad_inner)
        clad_inner_max = np.max(r_clad_inner)

        clad_outer_mean = np.mean(r_clad_outer)
        clad_outer_min = np.min(r_clad_outer)
        clad_outer_max = np.max(r_clad_outer)

        return (gap_mean, gap_min, gap_max,
                fuel_outer_mean, fuel_outer_min, fuel_outer_max,
                clad_inner_mean, clad_inner_min, clad_inner_max,
                clad_outer_mean, clad_outer_min, clad_outer_max)