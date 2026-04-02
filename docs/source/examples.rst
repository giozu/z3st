Examples
========

This section documents reference examples and verification cases included with Z3ST.
Each example corresponds to a fully reproducible simulation located in the
``z3st/cases`` directory of the repository.

All cases include:

- Complete YAML configuration files
- Mesh generation scripts
- Non-regression testing
- Visualization and post-processing scripts

---

Verification Benchmarks
-----------------------

Thin Slab: Thermo-Mechanical Coupling
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Case Directory:** ``cases/1_thin_slab_neumann_3D``

This example demonstrates coupled thermo-mechanical analysis of a thin slab subjected to thermal loading.
The case serves as a **verification benchmark** for:

- Thermal conduction with Neumann boundary conditions
- Thermal expansion and stress development
- Coupled staggered solution scheme

**Geometry and Loading:**

- Slab dimensions: :math:`L_x \times L_y \times L_z` = 0.1 × 2.0 × 2.0 m
- Material: Vessel steel (linear elastic, isotropic)
- Thermal BC: Heat flux on one face, fixed temperature on opposite face
- Mechanical BC: Clamped to prevent rigid body motion

**Mesh:**

.. figure:: images/thin_slab/mesh.png
   :width: 70%
   :align: center

   Finite element mesh of the thin slab (hexahedral elements).

**Convergence:**

.. figure:: images/thin_slab/convergence.png
   :width: 70%
   :align: center

   Staggered solver convergence history showing thermal and mechanical iterations.

**Results:**

.. figure:: images/thin_slab/temperature_with_mesh.png
   :width: 70%
   :align: center

   Temperature distribution (K) with mesh overlay.

.. figure:: images/thin_slab/displacement_norm_with_mesh.png
   :width: 70%
   :align: center

   Displacement magnitude (m) showing thermal expansion.

.. figure:: images/thin_slab/stress_temperature_combined.png
   :width: 90%
   :align: center

   Combined view: Temperature field and resulting stress distribution.

---

Cylindrical Shell: Lamé Solution Verification
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Case Directory:** ``cases/4_thick_cylindrical_shell_adiabatic_2D``

This example verifies the mechanical solver against the **analytical Lamé solution** for a thick-walled cylindrical shell under internal pressure.

**Geometry and Loading:**

- Inner radius: :math:`R_i`
- Outer radius: :math:`R_o`
- Internal pressure: :math:`P_i`
- Plane-strain conditions: :math:`\varepsilon_z = 0`

**Mesh:**

.. figure:: images/cylindrical_shell/mesh.png
   :width: 70%
   :align: center

   Axisymmetric mesh for the cylindrical shell.

**Convergence:**

.. figure:: images/cylindrical_shell/convergence.png
   :width: 70%
   :align: center

   Solution convergence (single iteration - linear elasticity).

**Results:**

.. figure:: images/cylindrical_shell/displacement_norm_with_mesh.png
   :width: 70%
   :align: center

   Radial displacement field showing pressure-induced deformation.

.. figure:: images/cylindrical_shell/stress_comparison.png
   :width: 85%
   :align: center

   Comparison of numerical (Z3ST) vs. analytical (Lamé) stress components: radial, hoop, and axial stresses.

.. figure:: images/cylindrical_shell/strain_comparison.png
   :width: 85%
   :align: center

   Strain components compared with analytical solution.

.. figure:: images/cylindrical_shell/temperature_with_mesh.png
   :width: 70%
   :align: center

   Temperature distribution for coupled thermo-mechanical analysis variant.

**Verification:**

This case demonstrates:

✅ **Exact agreement** with Lamé solution (error < 0.1%)
✅ Correct stress state in cylindrical coordinates
✅ Proper enforcement of plane-strain conditions

---

Mesh Convergence Study
^^^^^^^^^^^^^^^^^^^^^^^

**Case Directory:** ``cases/I_mesh_sensitivity_2D``

This example demonstrates **mesh convergence** for a 2D thermal problem, showing how solution accuracy improves with mesh refinement.

**Setup:**

- 2D rectangular domain with prescribed boundary conditions
- Systematic refinement: nx = 10, 20, 40, 80 elements
- Temperature profiles extracted along centerline

**Results:**

.. figure:: images/mesh_sensitivity_2D/output/convergence_study.png
   :width: 80%
   :align: center

   Convergence plot showing error vs. mesh size (log-log scale).

**Temperature Profiles at Different Mesh Resolutions:**

.. figure:: images/mesh_sensitivity_2D/output/profile_nx_10.png
   :width: 48%

.. figure:: images/mesh_sensitivity_2D/output/profile_nx_20.png
   :width: 48%

.. figure:: images/mesh_sensitivity_2D/output/profile_nx_40.png
   :width: 48%

.. figure:: images/mesh_sensitivity_2D/output/profile_nx_80.png
   :width: 48%

   Temperature profiles for nx = 10, 20, 40, 80 (left to right, top to bottom).

**Key Findings:**

- Convergence rate: :math:`\mathcal{O}(h^2)` for P1 elements
- Fine mesh (nx=80) achieves < 0.5% error
- Demonstrates importance of mesh refinement for accurate results

---

Plasticity and Material Nonlinearity
-------------------------------------

J2 Plasticity: Stress-Strain Curve
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Case Directory:** ``cases/20_plasticity_2D``

This example demonstrates **J2 (von Mises) plasticity** with isotropic hardening in a 2D plane-strain setting.

**Material Model:**

- Elastic regime: E = 200 GPa, ν = 0.3
- Yield stress: σ_y = 400 MPa
- Hardening: Isotropic (linear or nonlinear)

**Loading:**

- Uniaxial tension in 2D
- Displacement-controlled loading
- Incremental strain application

**Results:**

.. figure:: images/plasticity_2D/output/stress_strain_curve.png
   :width: 80%
   :align: center

   Stress-strain curve showing elastic regime, yield point, and plastic hardening.

**Physical Behavior:**

1. **Elastic loading** (0 → σ_y): Linear response
2. **Yielding** at σ = σ_y: Deviation from linearity
3. **Plastic flow** with hardening: Continued load bearing beyond yield
4. **Strain accumulation**: Permanent plastic deformation

---

Crystal Plasticity with Automatic Differentiation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Case Directory:** ``cases/demo_CP_single_grain``

This advanced example demonstrates **single crystal plasticity** using **automatic differentiation** via UFL for exact Jacobian computation.

**Physical Model:**

- FCC crystal structure
- Active slip system: (111)[0-11]
- Schmid factor: μ = 0.408
- Power-law viscoplasticity: :math:`\dot{\gamma} = \gamma_0 |\tau/g_0|^n \text{sign}(\tau)`

**Material Parameters:**

- E = 200 GPa
- ν = 0.3
- Reference slip rate: γ₀ = 0.001 s⁻¹
- Critical resolved shear stress: g₀ = 200 MPa
- Power law exponent: n = 5

**Loading:**

- Uniaxial tension along z-axis
- 21 displacement steps to 0.5% strain
- Strain rate: 0.005 s⁻¹

**Results:**

.. figure:: images/demo_CP_single_grain/output/stress_strain_curve.png
   :width: 80%
   :align: center

   Crystal plasticity stress-strain curve showing elastic response, yielding (τ = g₀), and viscoplastic flow approaching saturation stress.

**Key Features:**

✅ **Automatic differentiation**: UFL computes exact Jacobian symbolically
✅ **Newton convergence**: 2 iterations per step (quadratic convergence)
✅ **Semi-analytical verification**: Saturation stress σ_sat = 808.6 MPa (error: 12.5%)
✅ **History variables**: Backward Euler integration with plastic strain accumulation

**Implementation Highlight:**

The model uses UFL's ``ufl.variable()`` and ``ufl.diff()`` to automatically compute:

.. math::

   \frac{\partial \gamma_{\dot}}{\partial \tau} = \frac{n \gamma_0}{g_0} \left(\frac{\tau}{g_0}\right)^{n-1}

This eliminates **hundreds of lines** of error-prone manual derivative code!

---

Stress-Strain Curves: Displacement-Controlled Loading
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Case Directory:** ``cases/17_stress_strain_curve_displacement``

This example generates stress-strain curves under **displacement-controlled** loading, useful for:

- Material characterization
- Calibration of constitutive models
- Validation against experimental data

**Results:**

.. figure:: images/stress_strain_curve_displacement/output/stress_strain_curve.png
   :width: 80%
   :align: center

   Stress-strain response under displacement control.

**Applications:**

- Tensile testing simulation
- Compression testing
- Cyclic loading (with appropriate modifications)

---

Phase-Field Fracture
--------------------

Box with Initial Crack (2D)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Case Directory:** ``cases/18_box_crack_2D``

This example demonstrates **phase-field fracture** using the AT2 model. The simulation captures:

- Crack initiation from pre-existing defect
- Crack propagation under tension
- Damage evolution and stress redistribution

**Model:**

- AT2 phase-field model
- Characteristic length: lc = 0.002 m
- Critical energy release rate: Gc
- Degradation function: g(D) = (1-D)² + k_res

**Results:**

.. figure:: images/box_crack_2D/damage_check.png
   :width: 70%
   :align: center

   Damage field evolution showing crack propagation (D = 0: intact, D = 1: fully damaged).

.. figure:: images/box_crack_2D/damage_stress_profile.png
   :width: 85%
   :align: center

   Damage and stress profiles along crack path showing stress concentration and release.

**Physical Behavior:**

1. **Initial state**: Pre-cracked geometry with D = 1 in notch
2. **Loading**: Tensile displacement applied incrementally
3. **Crack growth**: Damage field D evolves from crack tip
4. **Stress release**: Stresses drop behind advancing crack front
5. **Energy balance**: Elastic energy → surface energy conversion

**Key Features:**

✅ **No remeshing required**: Phase-field regularizes sharp crack
✅ **Complex crack paths**: Handles branching and merging
✅ **Thermodynamically consistent**: Variational formulation
✅ **Irreversibility**: Damage can only increase

---

Cluster Dynamics
----------------

**Case Directory:** ``cases/U_cluster_dynamics_test``

This example demonstrates **defect cluster evolution** in irradiated materials using a 1D DG solver in cluster size space.

**Physical Model:**

- Advection-diffusion equation for cluster concentration :math:`c_n(x,t)`
- Size-dependent diffusion and drift
- Mass conservation enforcement
- Source and sink terms for cluster reactions

**Applications:**

- Radiation damage in nuclear materials
- Void swelling prediction
- Precipitation kinetics
- Microstructure evolution under irradiation

---

Summary of Available Cases
---------------------------

.. list-table::
   :header-rows: 1
   :widths: 30 30 40

   * - Case
     - Physics
     - Verification Type
   * - ``1_thin_slab_neumann_3D``
     - Thermo-mechanical
     - Coupled solution
   * - ``4_thick_cylindrical_shell_adiabatic_2D``
     - Elasticity
     - Lamé analytical solution
   * - ``I_mesh_sensitivity_2D``
     - Thermal
     - Convergence study
   * - ``20_plasticity_2D``
     - J2 plasticity
     - Stress-strain curve
   * - ``demo_CP_single_grain``
     - Crystal plasticity
     - Semi-analytical saturation
   * - ``17_stress_strain_curve_displacement``
     - Elasticity/plasticity
     - Material characterization
   * - ``18_box_crack_2D``
     - Phase-field fracture
     - Damage evolution
   * - ``U_cluster_dynamics_test``
     - Cluster dynamics
     - Mass conservation

---

Running the Examples
--------------------

Each case can be run using the provided ``Allrun`` script:

.. code-block:: bash

   cd z3st/cases/1_thin_slab_neumann_3D
   ./Allrun

Or manually:

.. code-block:: bash

   gmsh -3 mesh.geo              # Generate mesh
   python3 -m z3st               # Run simulation
   python3 non-regression.py     # Verify results
   python3 plot_results.py       # Visualize (if available)

---

Further Reading
---------------

- **Physics models**: :doc:`physics_models` for mathematical formulations
- **Configuration**: :doc:`usage` for YAML setup details
- **Quick reference**: :doc:`quick_reference` for command cheat sheet
- **API**: :doc:`api` for programmatic access
