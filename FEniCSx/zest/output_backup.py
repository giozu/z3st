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
import matplotlib
import matplotlib.pyplot
import matplotlib.pyplot as plt
import pyvista
import dolfinx
from dolfinx.mesh import compute_midpoints
import ufl
import numpy as np
import os

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

        figure_as_array = p.screenshot(f"output/{name}.png")

        return figure_as_array

    stress = problem.stress

    sigma_xx = stress[0, 0]
    sigma_xy = stress[0, 1]
    sigma_xz = stress[0, 2]

    sigma_yx = stress[1, 0]
    sigma_yy = stress[1, 1]
    sigma_yz = stress[1, 2]

    sigma_zx = stress[2, 0]
    sigma_zy = stress[2, 1]
    sigma_zz = stress[2, 2]

    strain = problem.strain

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

def plot_sigma_rr_at_z(problem, z_target=0.005, tol=0.0001):
    """
    Compute and plot the radial stress σ_rr as a function of radius at a specific height z_target.
    
    Parameters:
    - problem: The finite element problem containing the mesh and stress field.
    - z_target: The height (z-value) at which to extract and plot σ_rr.
    - tol: Tolerance to select points near the target z-value.
    """

    # Compute cell centers (for DG(0) function space)
    cell_centers = compute_midpoints(problem.mesh, problem.mesh.topology.dim, np.arange(problem.mesh.topology.index_map(problem.mesh.topology.dim).size_local))

    # Extract cell center coordinates
    x = cell_centers[:, 0]
    y = cell_centers[:, 1]
    z = cell_centers[:, 2]  # Heights
    r = np.sqrt(x**2 + y**2)  # Compute the radius

    unique_z_values = np.unique(cell_centers[:, 2])
    print("Unique z-values in the mesh:", unique_z_values)

    # Select only the points close to the chosen z-value
    mask = np.abs(z - z_target) < tol
    num_points = np.sum(mask)  # Count the number of selected points

    # Debug print
    print(f"Checking z = {z_target:.6f}: Found {num_points} points.")

    if num_points == 0:
        print(f"No points found at z = {z_target} within tolerance {tol}. Try increasing `tol`.")
        return

    # Define the function space for stress tensor (DG, 0)
    V_tensor = dolfinx.fem.functionspace(problem.mesh, ("DG", 0, (3, 3)))

    # Interpolate the stress tensor
    stress_expr = dolfinx.fem.Expression(problem.stress, V_tensor.element.interpolation_points())
    stress_func = dolfinx.fem.Function(V_tensor)
    stress_func.interpolate(stress_expr)

    # Extract the stress tensor values
    stress_values = stress_func.x.array.reshape(-1, 3, 3)  # Reshape for a 3x3 tensor

    # Compute the radial stress component σ_rr for selected points
    sigma_xx = stress_values[mask, 0, 0]  # σ_xx
    sigma_xy = stress_values[mask, 0, 1]  # σ_xy
    sigma_yy = stress_values[mask, 1, 1]  # σ_yy
    sigma_rr = (sigma_xx * x[mask]**2 + 2 * x[mask] * y[mask] * sigma_xy + sigma_yy * y[mask]**2) / (x[mask]**2 + y[mask]**2)

    # Plot the radial stress distribution
    plt.figure(figsize=(8, 6))
    plt.scatter(r[mask], sigma_rr, label=f"σ_rr at z={z_target:.4f}m", color="b")
    plt.xlabel("Radius [m]")
    plt.ylabel("Radial Stress σ_rr [Pa]")
    plt.legend()
    plt.grid()
    plt.title(f"Radial Stress σ_rr at z={z_target:.4f}m")
    plt.show()

def plot_sigma_rr_at_heights(problem, z_values=[0.005, 0.01, 0.015], tol=0.0001):
    """
    Compute and plot the radial stress σ_rr as a function of radius for multiple heights (z).
    
    Parameters:
    - problem: The finite element problem containing the mesh and stress field.
    - z_values: List of heights at which to extract and plot σ_rr.
    - tol: Tolerance to select points near the target z-values.
    """

    # Compute cell centers (for DG(0) function space)
    cell_centers = compute_midpoints(problem.mesh, problem.mesh.topology.dim, np.arange(problem.mesh.topology.index_map(problem.mesh.topology.dim).size_local))

    # Define the function space for stress tensor (DG, 0)
    V_tensor = dolfinx.fem.functionspace(problem.mesh, ("DG", 0, (3, 3)))

    # Interpolate the stress tensor
    stress_expr = dolfinx.fem.Expression(problem.stress, V_tensor.element.interpolation_points())
    stress_func = dolfinx.fem.Function(V_tensor)
    stress_func.interpolate(stress_expr)

    # Extract the stress tensor values
    stress_values = stress_func.x.array.reshape(-1, 3, 3)  # Reshape for a 3x3 tensor

    # Extract cell center coordinates
    x = cell_centers[:, 0]
    y = cell_centers[:, 1]
    z = cell_centers[:, 2]  # Heights
    r = np.sqrt(x**2 + y**2)  # Compute the radius

    # Create figure
    plt.figure(figsize=(8, 6))

    for z_target in z_values:
        # Select only the points close to the chosen z-value
        mask = np.abs(z - z_target) < tol
        num_points = np.sum(mask)  # Count the number of selected points

        print(f"Checking z = {z_target:.6f}: Found {num_points} points.")

        if num_points == 0:
            print(f"No points found at z = {z_target} within tolerance {tol}.")
            continue

        # Compute the radial stress component σ_rr
        sigma_xx = stress_values[mask, 0, 0]  # σ_xx
        sigma_xy = stress_values[mask, 0, 1]  # σ_xy
        sigma_yy = stress_values[mask, 1, 1]  # σ_yy
        sigma_rr = (sigma_xx * x[mask]**2 + 2 * x[mask] * y[mask] * sigma_xy + sigma_yy * y[mask]**2) / (x[mask]**2 + y[mask]**2)

        # Plot the profile for this height
        plt.scatter(r[mask], sigma_rr, label=f"σ_rr at z={z_target:.4f}m", alpha=0.7)

    # Final plot settings
    plt.xlabel("Radius [m]")
    plt.ylabel("Radial Stress σ_rr [Pa]")
    plt.legend()
    plt.grid()
    plt.title("Radial Stress σ_rr at Different Heights")
    plt.show()

def to_paraview(problem):
    """
    Export the displacement, stress, strain, and temperature fields to XDMF files for visualization in ParaView.

    Parameters:
    problem : object
        An object that contains the FEM problem with the fields:
        - displacement
        - stress
        - strain
        - temperature

    The resulting files will be:
    - 'displacement.xdmf' for the displacement field
    - 'stress_tensor.xdmf' for the stress field
    - 'strain_tensor.xdmf' for the strain field
    - 'temperature.xdmf' for the temperature field
    """
    
    # Create the output directory if it doesn't exist
    output_dir = "output"
    os.makedirs(output_dir, exist_ok=True)

    print(f"\nFunction to_paraview()")

    mesh = problem.mesh
    gdim = mesh.geometry.dim

    # Write displacement field to XDMF file
    # displacement_space = dolfinx.fem.functionspace(mesh, ("Lagrange", 1, (gdim,)))
    # displacement_expression = dolfinx.fem.Expression(problem.displacement, displacement_space.element.interpolation_points())
    # displacement_values = dolfinx.fem.Function(displacement_space)
    # displacement_values.interpolate(displacement_expression)
    output_file = os.path.join(output_dir, "displacement.xdmf")
    with dolfinx.io.XDMFFile(mesh.comm, output_file, "w") as xdmf:
        xdmf.write_mesh(mesh)
        xdmf.write_function(problem.displacement)
    print(f"Displacement field written to {output_file}")

    # Write stress tensor field to XDMF file
    stress_space = dolfinx.fem.functionspace(mesh, ("Lagrange", 1, (gdim, gdim)))
    stress_expression = dolfinx.fem.Expression(problem.stress, stress_space.element.interpolation_points())
    stress_values = dolfinx.fem.Function(stress_space)
    stress_values.interpolate(stress_expression)

    output_file = os.path.join(output_dir, "stress.xdmf")
    with dolfinx.io.XDMFFile(mesh.comm, output_file, "w") as xdmf:
        xdmf.write_mesh(mesh)
        stress_values.name = "stress"
        xdmf.write_function(stress_values)
    print(f"Stress tensor written to {output_file}")

    # Write strain tensor field to XDMF file
    strain_space = dolfinx.fem.functionspace(mesh, ("Lagrange", 1, (gdim, gdim)))
    strain_expression = dolfinx.fem.Expression(problem.strain, strain_space.element.interpolation_points())
    strain_values = dolfinx.fem.Function(strain_space)
    strain_values.interpolate(strain_expression)

    output_file = os.path.join(output_dir, "strain.xdmf")
    with dolfinx.io.XDMFFile(mesh.comm, output_file, "w") as xdmf:
        xdmf.write_mesh(mesh)
        strain_values.name = "strain"
        xdmf.write_function(strain_values)
    print(f"Strain tensor written to {output_file}")

    # Write temperature field to XDMF file
    output_file = os.path.join(output_dir, "temperature.xdmf")
    with dolfinx.io.XDMFFile(problem.mesh.comm, output_file, "w") as file:
        file.write_mesh(problem.mesh)
        file.write_function(problem.temperature)
    print(f"Temperature field written to {output_file}")

    # Complete output file (in progress)
    # Write results field to XDMF file
    # output_file = os.path.join(output_dir, "paraview_results.xdmf")
    # with dolfinx.io.XDMFFile(mesh.comm, output_file, "w") as xdmf:
    #     xdmf.write_mesh(mesh)
    #     problem.displacement.name = "displacement"
    #     problem.temperature.name = "temperature"
    #     stress_values.name = "stress"
    #     strain_values.name = "strain"
    #     xdmf.write_function(problem.displacement)
    #     xdmf.write_function(problem.temperature)
    #     xdmf.write_function(stress_values)
    #     xdmf.write_function(strain_values)

    # print(f"Displacement field written to {output_file}")
    # print(f"Temperature field written to {output_file}")
    # print(f"Stress tensor written to {output_file}")
    # print(f"Strain tensor written to {output_file}")
