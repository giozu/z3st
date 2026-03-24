Physics Models
==============

Z3ST provides a comprehensive set of physical models for multiphysics thermo-mechanical analysis, including:

- **Thermal conduction** with volumetric heating and gap conductance
- **Linear elasticity** and **J2 plasticity** for mechanical response
- **Phase-field fracture** (damage mechanics) for crack propagation
- **Cluster dynamics** for defect evolution in irradiated materials

All models are implemented using the FEniCSx variational framework and are designed to be modular, extensible, and coupled through a staggered solution scheme.

---

Thermal Model
-------------

The thermal model solves the **heat conduction equation** with support for:

- Volumetric heat sources (e.g., nuclear heating, internal dissipation)
- Temperature-dependent thermal conductivity
- Robin-type boundary conditions for gap conductance between subdomains
- Dirichlet and Neumann boundary conditions

Mathematical Formulation
^^^^^^^^^^^^^^^^^^^^^^^^

The governing equation for thermal diffusion in a domain :math:`\Omega` is:

.. math::

   \rho c_p \frac{\partial T}{\partial t} - \nabla \cdot (k \nabla T) = q \quad \text{in } \Omega

where:

- :math:`T` is the temperature field (K)
- :math:`\rho` is the material density (kg/m³)
- :math:`c_p` is the specific heat capacity (J/(kg·K))
- :math:`k` is the thermal conductivity (W/(m·K))
- :math:`q` is the volumetric heat source (W/m³)

**Boundary Conditions:**

- **Dirichlet:** :math:`T = T_0` on :math:`\Gamma_D`
- **Neumann:** :math:`-k \nabla T \cdot \mathbf{n} = q_0` on :math:`\Gamma_N`
- **Robin (gap conductance):** :math:`-k \nabla T \cdot \mathbf{n} = h(T - T_{\infty})` on :math:`\Gamma_R`

**Weak Formulation:**

Find :math:`T \in V_T` such that:

.. math::

   \int_\Omega \rho c_p \frac{\partial T}{\partial t} v \, \mathrm{d}\Omega
   + \int_\Omega k \nabla T \cdot \nabla v \, \mathrm{d}\Omega
   + \int_{\Gamma_R} h T v \, \mathrm{d}s
   = \int_\Omega q v \, \mathrm{d}\Omega
   + \int_{\Gamma_N} q_0 v \, \mathrm{d}s
   + \int_{\Gamma_R} h T_{\infty} v \, \mathrm{d}s

for all test functions :math:`v \in V_T`.

Implementation Notes
^^^^^^^^^^^^^^^^^^^^

- The thermal model is implemented in :class:`z3st.models.thermal_model.ThermalModel`
- Supports both **steady-state** and **transient** analysis
- Time integration uses backward Euler (implicit) scheme
- Volumetric heat sources can be spatially dependent (defined via UFL expressions)
- Material properties can be temperature-dependent through user-defined Python functions

**See also:**

- :doc:`examples` for thermal benchmark cases
- :doc:`usage` for YAML thermal model configuration

---

Mechanical Model
----------------

The mechanical model solves **linear elasticity** and **J2 plasticity** problems with support for:

- Linear elastic isotropic and anisotropic materials
- Thermo-mechanical coupling (thermal strain)
- Large deformation (finite strain) formulations (planned)
- Plasticity with isotropic hardening (J2 model)
- Custom differentiable constitutive laws

Mathematical Formulation: Linear Elasticity
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The governing equation for linear elasticity is the **equilibrium equation**:

.. math::

   \nabla \cdot \boldsymbol{\sigma} + \mathbf{f} = \mathbf{0} \quad \text{in } \Omega

where:

- :math:`\mathbf{u}` is the displacement field (m)
- :math:`\boldsymbol{\sigma}` is the Cauchy stress tensor (Pa)
- :math:`\mathbf{f}` is the body force (N/m³)

**Constitutive relation (isotropic linear elasticity):**

.. math::

   \boldsymbol{\sigma} = \lambda \, \mathrm{tr}(\boldsymbol{\varepsilon}) \mathbf{I} + 2\mu \boldsymbol{\varepsilon} - 3K\alpha_T (T - T_{\text{ref}}) \mathbf{I}

where:

- :math:`\boldsymbol{\varepsilon} = \frac{1}{2}(\nabla \mathbf{u} + \nabla \mathbf{u}^T)` is the strain tensor
- :math:`\lambda` and :math:`\mu` are the Lamé parameters
- :math:`K = \lambda + \frac{2\mu}{3}` is the bulk modulus
- :math:`\alpha_T` is the coefficient of thermal expansion (1/K)
- :math:`T_{\text{ref}}` is the reference temperature (K)

**Weak Formulation:**

Find :math:`\mathbf{u} \in V_u` such that:

.. math::

   \int_\Omega \boldsymbol{\sigma}(\mathbf{u}) : \nabla \mathbf{v} \, \mathrm{d}\Omega
   = \int_\Omega \mathbf{f} \cdot \mathbf{v} \, \mathrm{d}\Omega
   + \int_{\Gamma_N} \mathbf{t}_0 \cdot \mathbf{v} \, \mathrm{d}s

for all test functions :math:`\mathbf{v} \in V_u`.

J2 Plasticity (von Mises)
^^^^^^^^^^^^^^^^^^^^^^^^^

For elasto-plastic materials, the constitutive law follows the **J2 plasticity model** with isotropic hardening:

**Yield criterion:**

.. math::

   f(\boldsymbol{\sigma}, \varepsilon_p) = \sigma_{\text{eq}} - \sigma_y(\varepsilon_p) \leq 0

where:

- :math:`\sigma_{\text{eq}} = \sqrt{\frac{3}{2} \mathbf{s} : \mathbf{s}}` is the von Mises equivalent stress
- :math:`\mathbf{s} = \boldsymbol{\sigma} - \frac{1}{3}\mathrm{tr}(\boldsymbol{\sigma})\mathbf{I}` is the deviatoric stress
- :math:`\sigma_y(\varepsilon_p)` is the yield stress as a function of equivalent plastic strain
- :math:`\varepsilon_p` is the accumulated plastic strain

**Flow rule (associative plasticity):**

.. math::

   \dot{\boldsymbol{\varepsilon}}^p = \dot{\lambda} \frac{\partial f}{\partial \boldsymbol{\sigma}} = \dot{\lambda} \frac{3}{2} \frac{\mathbf{s}}{\sigma_{\text{eq}}}

where :math:`\dot{\lambda}` is the plastic multiplier determined by the consistency condition.

**Differentiable Constitutive Laws**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

One of the unique features of Z3ST is the ability to define **custom material laws as UFL expressions**, which are automatically differentiable. This enables:

- **Inverse modeling**: Identify material parameters from experimental data
- **Sensitivity analysis**: Compute gradients with respect to material properties
- **Advanced optimization**: Integrate with gradient-based optimizers

See :doc:`differentiable_features` for examples and applications.

Implementation Notes
^^^^^^^^^^^^^^^^^^^^

- The mechanical model is implemented in :class:`z3st.models.mechanical_model.MechanicalModel`
- Supports both **2D** (plane stress, plane strain) and **3D** geometries
- **Axisymmetric** formulation for rotationally symmetric problems
- Plasticity model implemented in :class:`z3st.models.plasticity_model`
- Linear solver options: direct (MUMPS), iterative (AMG, Hypre)

**See also:**

- :doc:`examples` for mechanical benchmark cases
- :doc:`usage` for YAML mechanical model configuration

---

Phase-Field Damage (Fracture Mechanics)
----------------------------------------

The phase-field damage model enables the simulation of **crack initiation, propagation, and branching** in brittle materials without requiring explicit crack tracking.

Z3ST implements the **AT1** and **AT2** models based on variational fracture mechanics, where the sharp crack is regularized by a continuous damage field :math:`D \in [0, 1]`.

Mathematical Formulation
^^^^^^^^^^^^^^^^^^^^^^^^

The total energy functional is:

.. math::

   \mathcal{E}(\mathbf{u}, D) = \int_\Omega \psi_e(\boldsymbol{\varepsilon}, D) \, \mathrm{d}\Omega + \mathcal{G}_c \int_\Omega \gamma(D, \nabla D) \, \mathrm{d}\Omega

where:

- :math:`\psi_e(\boldsymbol{\varepsilon}, D)` is the degraded elastic energy density
- :math:`\mathcal{G}_c` is the critical energy release rate (J/m²)
- :math:`\gamma(D, \nabla D)` is the crack surface density function

**Crack surface density (AT2 model):**

.. math::

   \gamma(D, \nabla D) = \frac{1}{2\ell_c} D^2 + \frac{\ell_c}{2} |\nabla D|^2

where :math:`\ell_c` is the characteristic length scale (m).

**Degradation function:**

.. math::

   g(D) = (1 - D)^2 + k_{\text{res}}

where :math:`k_{\text{res}}` is a small residual stiffness to avoid singularity.

**Governing equations:**

The system is solved using a **staggered scheme**:

1. **Mechanical problem** (with fixed :math:`D`):

   .. math::

      \nabla \cdot \left[ g(D) \boldsymbol{\sigma}_0(\boldsymbol{\varepsilon}) \right] = \mathbf{0}

2. **Damage problem** (with fixed :math:`\mathbf{u}`):

   .. math::

      -\ell_c^2 \nabla^2 D + D = \frac{\ell_c}{\mathcal{G}_c} \mathcal{H}

   where :math:`\mathcal{H}` is the crack driving force (elastic energy density).

**Irreversibility constraint:**

The damage field is constrained to be **monotonically increasing**:

.. math::

   D^{n+1} \geq D^n

This is enforced using a history field to track the maximum driving force.

Implementation Notes
^^^^^^^^^^^^^^^^^^^^

- The damage model is implemented in :class:`z3st.models.damage_model.DamageModel`
- Supports both **AT1** (linear crack density) and **AT2** (quadratic crack density)
- **Spectral decomposition** of strain energy for tension/compression split
- **Volumetric-deviatoric split** for alternative driving force definitions
- Irreversibility enforced via history field tracking

**See also:**

- :doc:`examples` for phase-field fracture cases
- :doc:`usage` for damage model configuration

---

Cluster Dynamics Model
-----------------------

The cluster dynamics model simulates **defect cluster evolution** in irradiated materials by solving a **1D advection-diffusion equation** in size space.

This model is particularly relevant for nuclear materials under neutron irradiation, where point defects (vacancies and interstitials) agglomerate into clusters of varying sizes.

Mathematical Formulation
^^^^^^^^^^^^^^^^^^^^^^^^

The governing equation for the cluster concentration :math:`c_n(x, t)` as a function of cluster size :math:`n` is:

.. math::

   \frac{\partial c_n}{\partial t} = -\frac{\partial J_n}{\partial n} + S_n - R_n

where:

- :math:`c_n(x, t)` is the concentration of clusters of size :math:`n` (m⁻³)
- :math:`J_n` is the flux of clusters in size space
- :math:`S_n` is the source term (production rate)
- :math:`R_n` is the sink term (recombination rate)

**Cluster flux in size space:**

.. math::

   J_n = \beta_n c_n - D_n \frac{\partial c_n}{\partial n}

where:

- :math:`\beta_n` is the bias factor (drift velocity)
- :math:`\D_n` is the diffusion coefficient in size space

**Boundary conditions:**

- **Nucleation boundary** (smallest clusters): prescribed flux or concentration
- **Large cluster limit**: absorbing or reflecting boundary

Implementation Notes
^^^^^^^^^^^^^^^^^^^^

- The cluster dynamics model is implemented in :class:`z3st.models.cluster_dynamic_model.ClusterDynamicModel`
- Uses **Discontinuous Galerkin (DG)** finite element discretization in size space
- Coupled to the thermal model via temperature-dependent reaction rates
- Mass conservation is enforced throughout the simulation

**See also:**

- :doc:`examples` for cluster dynamics test cases
- :doc:`usage` for cluster dynamics configuration

---

Gap Conductance Model
----------------------

The gap conductance model provides **thermal coupling between subdomains** through a **Robin-type boundary condition**. This is essential for multi-component systems where components are in thermal contact but not physically bonded.

Mathematical Formulation
^^^^^^^^^^^^^^^^^^^^^^^^

At the interface :math:`\Gamma_{\text{gap}}` between two subdomains, the heat flux is given by:

.. math::

   q = h_{\text{gap}} (T_1 - T_2)

where:

- :math:`h_{\text{gap}}` is the gap conductance (W/(m²·K))
- :math:`T_1, T_2` are temperatures on either side of the gap

**Robin boundary condition:**

.. math::

   -k_1 \nabla T_1 \cdot \mathbf{n} = h_{\text{gap}} (T_1 - T_2)

   -k_2 \nabla T_2 \cdot \mathbf{n} = h_{\text{gap}} (T_2 - T_1)

Implementation Notes
^^^^^^^^^^^^^^^^^^^^

- The gap model is implemented in :class:`z3st.models.gap_model.GapModel`
- Can use **fixed conductance** or **pressure-dependent** models
- Supports contact between multiple material domains
- Integrated into the thermal model weak form

**See also:**

- :doc:`usage` for gap model configuration

---

Coupled Thermo-Mechanical Analysis
-----------------------------------

Z3ST uses a **staggered coupling scheme** to solve coupled thermo-mechanical problems:

1. **Solve thermal problem** with fixed displacement field
2. **Solve mechanical problem** with updated temperature field
3. **Check convergence** of both fields
4. **Apply relaxation** and repeat until convergence

**Convergence criteria:**

.. math::

   \|T^{k+1} - T^k\| < \text{tol}_T

   \|\mathbf{u}^{k+1} - \mathbf{u}^k\| < \text{tol}_u

**Adaptive relaxation:**

The relaxation factors :math:`\omega_T` and :math:`\omega_u` can be dynamically adjusted based on convergence history to improve stability and convergence rate.

**See also:**

- :doc:`usage` for solver settings and coupling options
- :doc:`api` for solver implementation details
