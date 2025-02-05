"""

Zest: A FEniCSx module for fuel performance modeling
=====================================================
Author: Giovanni Zullo
Date: 16/01/2025

Description:
------------
output.py contains functions to visualize the results of the ZEST code.

"""

# Import necessary modules
import matplotlib.pyplot
import pyvista
import dolfinx
import matplotlib
import ufl
import numpy as np

def plot_results(displacement, problem):
    """
    Plot the displacement field of the problem
    """

    print(f"\nFunction plot_results()")

    pyvista.start_xvfb() # Start virtual framebuffer

    # Create plotter and pyvista grid
    p = pyvista.Plotter(off_screen=True)  # Enable off-screen rendering
    topology, cell_types, geometry = dolfinx.plot.vtk_mesh(problem.V_m)
    grid = pyvista.UnstructuredGrid(topology, cell_types, geometry)

    # Attach vector values to grid and warp grid by vector
    grid["u"] = displacement.x.array.reshape((geometry.shape[0], 3))

    p.add_mesh(grid, style="wireframe", color="k") # style = "points"
    warped = grid.warp_by_vector("u")
    p.add_mesh(warped, show_edges=True, cmap=matplotlib.colormaps.get_cmap("gist_rainbow").resampled(25), lighting=True)
    p.show_axes()

    figure_as_array = p.screenshot("displacement.png")
    return figure_as_array

def plot_temperature(temperature, problem):
    """
    Plot the temperature field of the problem
    """

    print(f"\nFunction plot_temperature()")

    pyvista.start_xvfb() # Start virtual framebuffer

    # Create plotter and pyvista grid
    topology, cell_types, geometry = dolfinx.plot.vtk_mesh(problem.V_t)
    grid = pyvista.UnstructuredGrid(topology, cell_types, geometry)

    # Attach vector values to grid and warp grid by vector
    grid.point_data["T"] = temperature.x.array.real
    grid.set_active_scalars("T")

    p = pyvista.Plotter(off_screen=True) # Enable off-screen rendering
    
    p.add_mesh(grid, show_edges=True)
    p.view_isometric()
    p.show_axes()

    figure_as_array = p.screenshot("temperature.png")
    return figure_as_array

def plot_stress_components(problem):
    """
    How to visualize the stress and strain fields in ParaView
    In UFL form, the stress and strain fields are defined as expressions (UFL expression objects) representing the symbolic form of the tensor.
    UFL is indeed designed to handle symbolic expressions for finite element methods, and it doesn't compute the actual values until it is evaluated on the mesh.
    """

    print(f"\nFunction plot_stress_components()")

    def return_tensor_component(name, component, V = dolfinx.fem.functionspace(problem.mesh, ("DG", 0)), mesh = problem.mesh):
        expr = dolfinx.fem.Expression(component, V.element.interpolation_points())
        field = dolfinx.fem.Function(V)
        return field.interpolate(expr)

    def plot_tensor_component(name, component, V = dolfinx.fem.functionspace(problem.mesh, ("DG", 0)), mesh = problem.mesh):
        """
        Plot a component of a tensor field
        """

        pyvista.start_xvfb()  # Start virtual framebuffer

        # Create plotter and pyvista grid
        p = pyvista.Plotter(off_screen=True)  # Enable off-screen rendering
        topology, cell_types, geometry = dolfinx.plot.vtk_mesh(problem.V_m)
        grid = pyvista.UnstructuredGrid(topology, cell_types, geometry)

        # Attach vector values to grid and warp grid by vector
        displacement = problem.displacement
        grid["u"] = displacement.x.array.reshape((geometry.shape[0], 3))

        p.add_mesh(grid, style="wireframe", color="k")  # style = "points"

        warped = grid.warp_by_vector("u")
        expr = dolfinx.fem.Expression(component, V.element.interpolation_points())
        field = dolfinx.fem.Function(V)
        field.interpolate(expr)
        
        warped.cell_data[name] = field.x.petsc_vec
        warped.set_active_scalars(name)

        p.add_mesh(warped, show_edges=True, cmap=matplotlib.colormaps.get_cmap("gist_rainbow").resampled(25), lighting=True)
        p.show_axes()

        figure_as_array = p.screenshot(f"{name}.png")

        return figure_as_array

    stress = problem.stress_field

    sigma_xx = stress[0, 0]
    sigma_xy = stress[0, 1]
    sigma_xz = stress[0, 2]

    sigma_yx = stress[1, 0]
    sigma_yy = stress[1, 1]
    sigma_yz = stress[1, 2]

    sigma_zx = stress[2, 0]
    sigma_zy = stress[2, 1]
    sigma_zz = stress[2, 2]

    strain = problem.strain_field

    strain_xx = strain[0, 0]
    strain_xy = strain[0, 1]
    strain_xz = strain[0, 2]

    strain_yx = strain[1, 0]
    strain_yy = strain[1, 1]
    strain_yz = strain[1, 2]

    strain_zx = strain[2, 0]
    strain_zy = strain[2, 1]
    strain_zz = strain[2, 2]

    print("sigma_xx ufl_shape:", sigma_xx.ufl_shape)
    print("strain_xx ufl_shape:", strain_xx.ufl_shape)

    plot_tensor_component("sigma_zz", sigma_zz)
    plot_tensor_component("strain_zz", strain_zz)

def plot_sigma_rr_at_z(problem, z_value=1.0):
    """
    Plot sigma_rr as a function of the radius r at a specific z value.
    """
    print(f"\nFunction plot_sigma_rr_at_z()")

    # Coordinate system
    x = ufl.SpatialCoordinate(problem.mesh)[0]
    y = ufl.SpatialCoordinate(problem.mesh)[1]
    z = ufl.SpatialCoordinate(problem.mesh)[2]

    # Stress components
    stress = problem.stress_field
    sigma_xx = stress[0, 0]
    sigma_xy = stress[0, 1]
    sigma_yy = stress[1, 1]

    # Radial stress component
    sigma_rr = (sigma_xx * x**2 + 2 * x * y * sigma_xy + sigma_yy * y**2) / (x**2 + y**2)

    # Function space for sigma_rr
    V_scalar = dolfinx.fem.functionspace(problem.mesh, ("DG", 0))
    sigma_rr_expr = dolfinx.fem.Expression(sigma_rr, V_scalar.element.interpolation_points())
    sigma_rr_func = dolfinx.fem.Function(V_scalar)
    sigma_rr_func.interpolate(sigma_rr_expr)

    # Extract values at z = z_value
    dofs_z = dolfinx.fem.locate_dofs_geometrical(V_scalar, lambda x: np.isclose(x[2], z_value))
    sigma_rr_values = sigma_rr_func.x.array[dofs_z]

    # Compute radius values
    coordinates = problem.mesh.geometry.x[dofs_z]
    radius_values = np.sqrt(coordinates[:, 0]**2 + coordinates[:, 1]**2)

    # Plot sigma_rr as a function of radius
    import matplotlib.pyplot as plt

    plt.figure()
    plt.plot(radius_values, sigma_rr_values, 'o')
    plt.xlabel("Radius (r)")
    plt.ylabel("sigma_rr")
    plt.title(f"Radial stress component sigma_rr at z = {z_value}")
    plt.savefig(f"sigma_rr_z_{z_value}.png")
    plt.show()

def plot_sigma_xx_at_z(problem, z_value=1.0):
    """
    Plot sigma_xx as a function of the radius r at a specific z value.
    """
    print(f"\nFunction plot_sigma_xx_at_z()")

    # Coordinate system
    x = ufl.SpatialCoordinate(problem.mesh)[0]
    y = ufl.SpatialCoordinate(problem.mesh)[1]
    z = ufl.SpatialCoordinate(problem.mesh)[2]

    # Stress components
    stress = problem.stress_field
    sigma_xx = stress[0, 0]

    # Function space for sigma_xx
    V_scalar = dolfinx.fem.functionspace(problem.mesh, ("DG", 0))
    sigma_xx_expr = dolfinx.fem.Expression(sigma_xx, V_scalar.element.interpolation_points())
    sigma_xx_func = dolfinx.fem.Function(V_scalar)
    sigma_xx_func.interpolate(sigma_xx_expr)

    # Extract values at z = z_value
    dofs_z = dolfinx.fem.locate_dofs_geometrical(V_scalar, lambda x: np.isclose(x[2], z_value))
    sigma_xx_values = sigma_xx_func.x.array[dofs_z]

    # Compute radius values
    coordinates = problem.mesh.geometry.x[dofs_z]
    radius_values = np.sqrt(coordinates[:, 0]**2 + coordinates[:, 1]**2)
    
    # Plot sigma_xx as a function of radius
    import matplotlib.pyplot as plt

    plt.figure()
    plt.plot(radius_values, sigma_xx_values, 'o')
    plt.xlabel("Radius (r)")
    plt.ylabel("sigma_xx")
    plt.title(f"Stress component sigma_xx at z = {z_value}")
    plt.savefig(f"sigma_xx_z_{z_value}.png")
    plt.show()

def print_x_at_z_zero(problem):
    """
    Print the values of x at z = 0.
    """
    z_value = 0.0
    coordinates = problem.mesh.geometry.x
    X_coord2 = Xh0.x.array.reshape(len(mesh.geometry.x[:,1]),2)

    x = coordinates[:, 0]
    y = coordinates[:, 1]
    z = coordinates[:, 2]

    print(x)
    print(y)
    print(z)

    # x_at_z_zero = coordinates[np.isclose(coordinates[:, 2], z_value)]
    # print(f"Values of x at z = {z_value}:")
    # print(x_at_z_zero)

def to_paraview(problem):
    """
    Write the displacement, stress, and strain fields to XDMF files for visualization in ParaView
    """

    print(f"\nFunction to_paraview()")

    mesh = problem.mesh
    gdim = mesh.geometry.dim
    mesh_nodes = len(problem.mesh.geometry.x)

    # Define the cell type as a string (e.g., "triangle" or "tetrahedron")
    # cell_type = "tetrahedron"  # Update this based on your mesh type
    # cell_type = "hexahedron"
    # cell_type = mesh.ufl_cell().cellname()

    output_file = "displacement.xdmf"
    with dolfinx.io.XDMFFile(mesh.comm, output_file, "w") as xdmf:
        xdmf.write_mesh(mesh)
        problem.displacement.name = "displacement"
        xdmf.write_function(problem.displacement)
    print(f"Displacement field written to {output_file}")

    stress_space = dolfinx.fem.functionspace(mesh, ("Lagrange", 1, (gdim, gdim)))
    stress_expression = dolfinx.fem.Expression(problem.stress_field, stress_space.element.interpolation_points())
    stress_values = dolfinx.fem.Function(stress_space)
    stress_values.interpolate(stress_expression)

    output_file = "stress_tensor.xdmf"
    with dolfinx.io.XDMFFile(mesh.comm, output_file, "w") as xdmf:
        xdmf.write_mesh(mesh)
        stress_values.name = "stress"
        xdmf.write_function(stress_values)
    print(f"Stress tensor written to {output_file}")

    strain_space = dolfinx.fem.functionspace(mesh, ("Lagrange", 1, (gdim, gdim)))
    strain_expression = dolfinx.fem.Expression(problem.strain_field, strain_space.element.interpolation_points())
    strain_values = dolfinx.fem.Function(strain_space)
    strain_values.interpolate(strain_expression)

    output_file = "strain_tensor.xdmf"
    with dolfinx.io.XDMFFile(mesh.comm, output_file, "w") as xdmf:
        xdmf.write_mesh(mesh)
        strain_values.name = "strain"
        xdmf.write_function(strain_values)
    print(f"Strain tensor written to {output_file}")
