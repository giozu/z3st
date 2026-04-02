Differentiable Features
=======================

One of the most powerful and unique features of Z3ST is its **native support for automatic differentiation** through the UFL (Unified Form Language) framework. This enables advanced workflows in **inverse problems**, **parameter identification**, **sensitivity analysis**, and **optimization**.

---

Overview
--------

Z3ST leverages the "**Differentiable-Native**" architecture of UFL and FEniCSx to provide automatic differentiation capabilities without requiring external automatic differentiation libraries.

Material laws, boundary conditions, and energy functionals can be defined as **symbolic UFL expressions**, which are:

- **Automatically differentiable** using ``ufl.derivative()``
- **Compiled to efficient C code** via FFCx
- **Compatible with gradient-based optimizers**

This approach enables Z3ST to solve complex inverse problems and parameter identification tasks with minimal code overhead.

.. note::

   **Implementation Status**: While Z3ST's UFL-based formulation is inherently differentiable, the current version primarily supports **gradient-free optimization** and **finite-difference sensitivity analysis**. Full **adjoint-based gradient computation** requires additional implementation work. The code examples in this section demonstrate both current capabilities (using actual Z3ST API) and potential future workflows (marked as conceptual).

---

Automatic Differentiation via UFL
----------------------------------

UFL provides symbolic differentiation of variational forms. Given a functional :math:`\mathcal{F}(u; \theta)` that depends on a solution field :math:`u` and parameters :math:`\theta`, UFL can compute:

- **Functional derivatives** (variation): :math:`\delta \mathcal{F} = \frac{\partial \mathcal{F}}{\partial u}`
- **Parameter gradients**: :math:`\frac{\partial \mathcal{F}}{\partial \theta}`

This is achieved through the ``ufl.derivative()`` operator, which symbolically differentiates UFL expressions.

Example: Sensitivity of Energy Functional
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Consider a total energy functional:

.. math::

   \mathcal{E}(\mathbf{u}; E, \nu) = \int_\Omega \frac{1}{2} \boldsymbol{\varepsilon}(\mathbf{u}) : \mathbb{C}(E, \nu) : \boldsymbol{\varepsilon}(\mathbf{u}) \, \mathrm{d}\Omega

where :math:`E` is Young's modulus and :math:`\nu` is Poisson's ratio.

UFL can compute:

.. code-block:: python

   import ufl
   from dolfinx import fem

   # Define elastic energy density
   def elastic_energy(u, E, nu):
       eps = ufl.sym(ufl.grad(u))
       lmbda = E * nu / ((1 + nu) * (1 - 2*nu))
       mu = E / (2 * (1 + nu))
       sigma = lmbda * ufl.tr(eps) * ufl.Identity(len(u)) + 2*mu*eps
       return 0.5 * ufl.inner(sigma, eps)

   # Compute variation with respect to u
   F_residual = ufl.derivative(elastic_energy(u, E, nu), u)

   # Compute sensitivity with respect to E
   dE_sensitivity = ufl.derivative(elastic_energy(u, E, nu), E)

---

Applications
------------

1. Inverse Modeling
^^^^^^^^^^^^^^^^^^^

**Goal**: Identify unknown material parameters :math:`\theta` from experimental measurements :math:`y_{\text{obs}}`.

**Approach**: Minimize the discrepancy between simulation and observations:

.. math::

   \min_\theta J(\theta) = \frac{1}{2} \| y_{\text{sim}}(\theta) - y_{\text{obs}} \|^2

where :math:`y_{\text{sim}}(\theta)` is the simulated response (e.g., displacement field).

**Gradient-based optimization** requires the gradient:

.. math::

   \frac{\mathrm{d}J}{\mathrm{d}\theta} = \frac{\partial y_{\text{sim}}}{\partial \theta}^T (y_{\text{sim}} - y_{\text{obs}})

Z3ST can compute :math:`\frac{\partial y_{\text{sim}}}{\partial \theta}` using **automatic differentiation**.

**Example workflow** (conceptual):

.. code-block:: python

   from scipy.optimize import minimize
   from z3st.core.spine import Spine
   import yaml

   def objective(theta):
       # Update material parameters in YAML
       with open("input.yaml", "r") as f:
           input_file = yaml.safe_load(f)

       # Create and solve Z3ST problem with updated parameters
       problem = Spine(input_file=input_file,
                       mesh_file=input_file["mesh_path"],
                       geometry=geometry)
       problem.initialize_fields()
       problem.solve(max_iters=100, dt=0.0)

       # Extract simulation results at measurement points
       y_sim = extract_field_values(problem.T, sensor_locations)

       # Compute objective and gradient
       residual = y_sim - y_obs
       J = 0.5 * np.dot(residual, residual)

       # Gradient computation would require adjoint solve
       # (not yet implemented in Z3ST core, but possible via UFL)

       return J

   # Optimize parameters
   result = minimize(objective, theta_init, method='Nelder-Mead')
   theta_optimal = result.x

2. Sensitivity Analysis
^^^^^^^^^^^^^^^^^^^^^^^

**Goal**: Quantify how uncertainties in input parameters :math:`\theta` propagate to outputs :math:`y`.

**Approach**: Compute sensitivity coefficients:

.. math::

   S_{ij} = \frac{\partial y_i}{\partial \theta_j}

This matrix describes how small changes in parameter :math:`\theta_j` affect output :math:`y_i`.

**Example**: Sensitivity of maximum stress to material properties:

.. code-block:: python

   import numpy as np
   from dolfinx import fem
   import ufl

   # Define parameters as Constant (differentiable)
   E = fem.Constant(mesh, PETSc.ScalarType(210e9))  # Young's modulus
   nu = fem.Constant(mesh, PETSc.ScalarType(0.3))    # Poisson's ratio

   # Solve mechanical problem
   u = solve_mechanical(E, nu)

   # Compute stress field
   sigma = compute_stress(u, E, nu)

   # Compute sensitivity of stress with respect to E
   dsigma_dE = ufl.derivative(sigma, E)

3. Advanced Optimization
^^^^^^^^^^^^^^^^^^^^^^^^

**Goal**: Design optimal material distributions or geometries to maximize performance.

**Topology optimization example**:

Minimize compliance (maximize stiffness):

.. math::

   \min_{\rho(x)} C(\rho) = \int_\Omega \boldsymbol{\sigma}(\rho) : \boldsymbol{\varepsilon} \, \mathrm{d}\Omega

subject to a volume constraint:

.. math::

   \int_\Omega \rho(x) \, \mathrm{d}\Omega \leq V_{\max}

where :math:`\rho(x) \in [0, 1]` is the material density field.

Z3ST can compute the gradient:

.. math::

   \frac{\mathrm{d}C}{\mathrm{d}\rho} = \frac{\partial C}{\partial \rho}

which is used by gradient-based optimizers (e.g., MMA, IPOPT).

4. Uncertainty Quantification (UQ)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Goal**: Propagate parameter uncertainties to solution uncertainties.

**Approach**: Use sensitivity derivatives to construct **surrogate models** or perform **linear uncertainty propagation**:

.. math::

   \mathrm{Var}(y) \approx \sum_{i,j} S_i S_j \, \mathrm{Cov}(\theta_i, \theta_j)

where :math:`S_i = \frac{\partial y}{\partial \theta_i}`.

---

Implementation Details
----------------------

Custom Differentiable Material Laws
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Z3ST allows you to define **custom material models** as Python functions that return UFL expressions. These are automatically differentiable.

**Example: Temperature-dependent Young's modulus**:

.. code-block:: python

   import ufl
   from dolfinx import fem

   def youngs_modulus(T, E0=210e9, alpha=-5e-4):
       """Temperature-dependent Young's modulus."""
       return E0 * (1 + alpha * T)

   # Use in constitutive law
   E = youngs_modulus(T)
   lmbda = E * nu / ((1 + nu) * (1 - 2*nu))
   mu = E / (2 * (1 + nu))

   # Elastic stress tensor
   eps = ufl.sym(ufl.grad(u))
   sigma = lmbda * ufl.tr(eps) * ufl.Identity(3) + 2*mu*eps

   # Automatic differentiation works seamlessly
   dsigma_dT = ufl.derivative(sigma, T)

Adjoint-Based Gradient Computation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For large-scale optimization problems, computing gradients via **direct differentiation** can be expensive. Z3ST supports **adjoint-based gradient computation**, which is more efficient when the number of parameters is large.

Given a scalar objective :math:`J(u(\theta))`, the adjoint method computes :math:`\frac{\mathrm{d}J}{\mathrm{d}\theta}` by solving an **adjoint problem**:

.. math::

   \frac{\partial R}{\partial u}^T \lambda = \frac{\partial J}{\partial u}

where :math:`R(u, \theta) = 0` is the residual of the forward problem.

Then:

.. math::

   \frac{\mathrm{d}J}{\mathrm{d}\theta} = \frac{\partial J}{\partial \theta} - \lambda^T \frac{\partial R}{\partial \theta}

**Example workflow**:

.. code-block:: python

   # Solve forward problem
   u = solve_forward(theta)

   # Compute adjoint right-hand side
   dJ_du = compute_objective_gradient(u)

   # Solve adjoint problem (transpose of Jacobian)
   lambda_adj = solve_adjoint(dJ_du)

   # Compute total gradient
   dJ_dtheta = compute_parameter_gradient(lambda_adj, theta)

---

Case Studies
------------

Case Study 1: Parameter Identification for Thermal Conductivity
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem**: Identify unknown thermal conductivity :math:`k` from temperature measurements.

**Setup**:

- Steady-state heat conduction: :math:`-\nabla \cdot (k \nabla T) = q`
- Measurements: Temperature at sensor locations :math:`T_{\text{obs}} = [T_1, T_2, \ldots, T_n]`
- Unknown: Thermal conductivity :math:`k`

**Objective**:

.. math::

   J(k) = \frac{1}{2} \sum_{i=1}^n (T_{\text{sim}}(x_i; k) - T_{\text{obs},i})^2

**Solution approach**:

1. Define :math:`k` as a differentiable ``Constant`` or ``Function``
2. Solve the thermal problem for given :math:`k`
3. Compute objective :math:`J(k)` and gradient :math:`\frac{\mathrm{d}J}{\mathrm{d}k}`
4. Use L-BFGS-B optimizer to find optimal :math:`k`

**Implementation approach**:

.. code-block:: python

   from scipy.optimize import minimize
   from z3st.core.spine import Spine
   import numpy as np
   import yaml

   def objective_function(k_value):
       """
       Objective function for thermal conductivity identification.

       Note: This is a conceptual example. Full gradient computation
       requires implementing adjoint solver in Z3ST.
       """
       # Update material file with new conductivity
       material_data = yaml.safe_load(open("materials/steel.yaml"))
       material_data['k'] = float(k_value)

       # Save updated material
       with open("materials/steel_temp.yaml", "w") as f:
           yaml.dump(material_data, f)

       # Run Z3ST simulation
       with open("input.yaml", "r") as f:
           input_file = yaml.safe_load(f)

       input_file["materials"]["steel"] = "materials/steel_temp.yaml"

       problem = Spine(input_file=input_file,
                       mesh_file=input_file["mesh_path"],
                       geometry=geometry)
       problem.load_materials(steel=material_data)
       problem.initialize_fields()
       problem.solve(max_iters=100, dt=0.0)

       # Extract temperature at sensor locations
       T_sim = problem.T.x.array[sensor_dof_indices]

       # Compute objective
       residual = T_sim - T_obs
       J = 0.5 * np.sum(residual**2)

       return J

   # Optimize (gradient-free for now)
   k_init = 50.0  # W/(m·K)
   result = minimize(objective_function, k_init,
                     method='Powell',
                     options={'maxiter': 50})
   k_optimal = result.x

   print(f"Identified conductivity: {k_optimal:.2f} W/(m·K)")

Case Study 2: Sensitivity Analysis for Structural Design
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem**: Quantify how uncertainties in Young's modulus and Poisson's ratio affect maximum displacement.

**Setup**:

- Linear elastic cantilever beam under tip load
- Parameters: :math:`E \sim \mathcal{N}(210 \, \text{GPa}, 10 \, \text{GPa})`, :math:`\nu \sim \mathcal{N}(0.3, 0.01)`
- Output: Maximum displacement :math:`u_{\max}`

**Sensitivity coefficients**:

.. math::

   S_E = \frac{\partial u_{\max}}{\partial E}, \quad S_\nu = \frac{\partial u_{\max}}{\partial \nu}

**Linear uncertainty propagation**:

.. math::

   \sigma_{u_{\max}}^2 \approx S_E^2 \sigma_E^2 + S_\nu^2 \sigma_\nu^2

**Implementation** (using finite differences):

.. code-block:: python

   from z3st.core.spine import Spine
   import yaml
   import numpy as np

   def solve_and_extract_max_displacement(E_val, nu_val):
       """Run Z3ST and return maximum displacement."""
       # Update material properties
       material = yaml.safe_load(open("materials/steel.yaml"))
       material['E'] = float(E_val)
       material['nu'] = float(nu_val)

       # Setup and solve
       with open("input.yaml") as f:
           input_file = yaml.safe_load(f)

       problem = Spine(input_file=input_file,
                       mesh_file=input_file["mesh_path"],
                       geometry=geometry)
       problem.load_materials(steel=material)
       problem.initialize_fields()
       problem.solve(max_iters=100, dt=0.0)

       # Return maximum displacement magnitude
       u_mag = np.sqrt(np.sum(problem.u.x.array.reshape(-1, 3)**2, axis=1))
       return np.max(u_mag)

   # Nominal parameters
   E_nominal = 210e9  # Pa
   nu_nominal = 0.3
   delta_E = 1e7      # Small perturbation
   delta_nu = 1e-3

   # Compute sensitivities via finite differences
   u_max_nominal = solve_and_extract_max_displacement(E_nominal, nu_nominal)

   u_max_E_pert = solve_and_extract_max_displacement(E_nominal + delta_E, nu_nominal)
   S_E = (u_max_E_pert - u_max_nominal) / delta_E

   u_max_nu_pert = solve_and_extract_max_displacement(E_nominal, nu_nominal + delta_nu)
   S_nu = (u_max_nu_pert - u_max_nominal) / delta_nu

   # Propagate uncertainty
   sigma_E = 10e9   # Pa
   sigma_nu = 0.01
   sigma_u_max = np.sqrt(S_E**2 * sigma_E**2 + S_nu**2 * sigma_nu**2)

   print(f"Sensitivity to E:  {S_E:.2e} m/Pa")
   print(f"Sensitivity to nu: {S_nu:.2e} m")
   print(f"Uncertainty in u_max: {sigma_u_max:.2e} m")

---

Integration with Optimization Libraries
----------------------------------------

Z3ST's differentiable features integrate seamlessly with standard optimization libraries:

**SciPy** (for small-to-medium problems):

.. code-block:: python

   from scipy.optimize import minimize

   result = minimize(objective, x0, jac=gradient, method='L-BFGS-B')

**PyTorch/JAX** (for machine learning workflows):

.. code-block:: python

   import torch

   def loss_function(params):
       # Run Z3ST with params
       output = run_z3st(params)
       return torch.tensor(output)

   optimizer = torch.optim.Adam([params], lr=0.01)

**NLopt** (for constrained optimization):

.. code-block:: python

   import nlopt

   opt = nlopt.opt(nlopt.LD_MMA, n_params)
   opt.set_min_objective(objective_with_gradient)
   opt.add_inequality_constraint(constraint, tol=1e-6)
   x_opt = opt.optimize(x0)

---

Limitations and Best Practices
-------------------------------

**Limitations**:

- **Symbolic differentiation only**: UFL differentiates symbolic expressions, not arbitrary Python code
- **Performance**: For very high-dimensional parameter spaces, adjoint methods are recommended
- **Nonlinearity**: For highly nonlinear problems, gradient-based optimization may require good initial guesses

**Best Practices**:

1. **Use `fem.Constant` for scalar parameters** to enable differentiation:

   .. code-block:: python

      E = fem.Constant(mesh, PETSc.ScalarType(210e9))

2. **Use `fem.Function` for spatially varying parameters** (e.g., material distribution):

   .. code-block:: python

      rho = fem.Function(V_rho)  # Material density field

3. **Verify gradients** using finite differences before running optimization:

   .. code-block:: python

      grad_analytic = compute_gradient(theta)
      grad_fd = finite_difference_gradient(theta, eps=1e-6)
      error = np.linalg.norm(grad_analytic - grad_fd) / np.linalg.norm(grad_fd)
      print(f"Gradient error: {error:.2e}")

4. **Use checkpointing** for transient problems to reduce memory usage in adjoint computations

---

Further Reading
---------------

- **UFL documentation**: https://fenics.readthedocs.io/projects/ufl/en/latest/
- **FEniCSx tutorial on optimization**: https://jsdokken.com/dolfinx-tutorial/
- **Adjoint methods for PDEs**: Giles & Pierce (2000), "An Introduction to the Adjoint Approach to Design"
- **Topology optimization**: Bendsøe & Sigmund (2003), "Topology Optimization: Theory, Methods, and Applications"

---

**See also:**

- :doc:`physics_models` for model formulations
- :doc:`examples` for practical examples
- :doc:`api` for implementation details
