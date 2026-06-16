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

Linear thermo-mechanics
-----------------------

Thin Slab: Thermo-Mechanical Coupling
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. **Case Directory:** ``cases/verification/thermal/thin_slab_neumann_3D``

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

Cylindrical Shell: Lamé Solution Verification
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. **Case Directory:** ``cases/verification/thermal/thick_cylindrical_shell_adiabatic_2D``

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

Agreement with Lamé solution
Correct stress state in cylindrical coordinates
Proper enforcement of plane-strain conditions

Plasticity and Material Nonlinearity
-------------------------------------

J2 Plasticity: Stress-Strain Curve
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Case Directory:** ``cases/verification/plasticity/j2_hardening_2D``

This example demonstrates **J2 (von Mises) plasticity** with isotropic hardening in a 2D plane-strain setting.

**Material Model:**

- Elastic regime: :math:`E = 200\,\mathrm{GPa}`, :math:`\nu = 0.3`
- Yield stress: :math:`\sigma_y = 400\,\mathrm{MPa}`
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

1. **Elastic loading** (:math:`0 \to \sigma_y`): Linear response
2. **Yielding** at :math:`\sigma = \sigma_y`: Deviation from linearity
3. **Plastic flow** with hardening: Continued load bearing beyond yield
4. **Strain accumulation**: Permanent plastic deformation

Crystal Plasticity with Automatic Differentiation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Case Directory:** ``cases/verification/plasticity/crystal_single_grain``

This advanced example demonstrates **single crystal plasticity** using **automatic differentiation** via UFL for exact Jacobian computation.

**Physical Model:**

- FCC crystal structure
- Active slip system: (111)[0-11]
- Schmid factor: :math:`\mu = 0.408`
- Power-law viscoplasticity: :math:`\dot{\gamma} = \gamma_0 |\tau/g_0|^n \text{sign}(\tau)`

**Material Parameters:**

- :math:`E = 200\,\mathrm{GPa}`
- :math:`\nu = 0.3`
- Reference slip rate: :math:`\dot{\gamma}_0 = 0.001\,\mathrm{s}^{-1}`
- Critical resolved shear stress: :math:`g_0 = 200\,\mathrm{MPa}`
- Power law exponent: :math:`n = 5`

**Loading:**

- Uniaxial tension along z-axis
- 21 displacement steps to 0.5% strain
- Strain rate: :math:`0.005\,\mathrm{s}^{-1}`

**Results:**

.. figure:: images/demo_CP_single_grain/output/stress_strain_curve.png
   :width: 80%
   :align: center

   Crystal plasticity stress-strain curve showing elastic response, yielding (:math:`\tau = g_0`), and viscoplastic flow approaching saturation stress.

**Key Features:**

**Automatic differentiation**: UFL computes exact Jacobian symbolically
**Newton convergence**: 2 iterations per step (quadratic convergence)
**Semi-analytical verification**: Saturation stress :math:`\sigma_{sat} = 808.6\,\mathrm{MPa}` (error: 12.5%)
**History variables**: Backward Euler integration with plastic strain accumulation

**Implementation Highlight:**

The model uses UFL's ``ufl.variable()`` and ``ufl.diff()`` to automatically compute:

.. math::

   \frac{\partial \dot{\gamma}}{\partial \tau} = \frac{n \gamma_0}{g_0} \left(\frac{\tau}{g_0}\right)^{n-1}

This eliminates **hundreds of lines** of error-prone manual derivative code!

Stress-Strain Curves: Displacement-Controlled Loading
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Case Directory:** ``cases/verification/plasticity/stress_strain_displacement``

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

Phase-Field Fracture
--------------------

Box with Initial Crack (2D)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Case Directory:** ``cases/regression/box_crack_2D``

This example demonstrates **phase-field fracture** using the AT2 model. The simulation captures:

- Crack initiation from pre-existing defect
- Crack propagation under tension
- Damage evolution and stress redistribution

**Model:**

- AT2 phase-field model
- Characteristic length: :math:`\ell_c = 0.002\,\mathrm{m}`
- Critical energy release rate: :math:`G_c`
- Degradation function: :math:`g(D) = (1-D)^2 + k_{res}`

**Results:**

.. figure:: images/box_crack_2D/damage_check.png
   :width: 70%
   :align: center

   Damage field evolution showing crack propagation (:math:`D = 0`: intact, :math:`D = 1`: fully damaged).

.. figure:: images/box_crack_2D/damage_stress_profile.png
   :width: 85%
   :align: center

   Damage and stress profiles along crack path showing stress concentration and release.

**Physical Behavior:**

1. **Initial state**: Pre-cracked geometry with :math:`D = 1` in notch
2. **Loading**: Tensile displacement applied incrementally
3. **Crack growth**: Damage field D evolves from crack tip
4. **Stress release**: Stresses drop behind advancing crack front
5. **Energy balance**: Elastic energy → surface energy conversion

**Key Features:**

**No remeshing required**: Phase-field regularizes sharp crack
**Complex crack paths**: Handles branching and merging
**Thermodynamically consistent**: Variational formulation
**Irreversibility**: Damage can only increase

Single-Edge-Notched Shear Test (Miehe 2010 benchmark)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Case Directory:** ``cases/benchmarks/sen_shear``

A standard phase-field-fracture benchmark: a square specimen with a horizontal
pre-notch is loaded in shear. The curved crack path that develops reproduces the
result of Miehe et al., *Comput. Methods Appl. Mech. Engrg.* 199 (2010)
2765--2778 -- the foundational phase-field paper.

**Setup:**

- Coupled mechanics + phase-field damage (steel)
- Horizontal pre-notch on the left edge
- Imposed horizontal displacement :math:`u_x = 15\,\mu\mathrm{m}`

**Results:**

.. figure:: images/sen_shear/SENS_ux.png
   :width: 48%

.. figure:: images/sen_shear/SENS_damage_final.png
   :width: 48%

   Imposed shear displacement (left) and the resulting curved crack path (right).

.. figure:: images/sen_shear/SENS_vonMises.png
   :width: 55%
   :align: center

   Von Mises stress concentration at the propagating crack tip.

UO\ :sub:`2` Pellet Thermal-Shock Cracking (McClenny 2022)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Case Directory:** ``cases/benchmarks/pellet_quench_2D_xy``

The fully coupled showpiece, in which thermal, mechanical, and phase-field damage
all meet. A 2D plane-strain transverse cross-section of a UO\ :sub:`2` pellet has
a cold-contact wedge on its rim; the tensile hoop-stress ring drives discrete
radial cracks. It reproduces McClenny et al., *J. Nucl. Mater.* 565 (2022).

**Model:**

- AT1 crack density + Amor split + Ambati hybrid constraint
- Fully coupled :math:`T \to \boldsymbol{\varepsilon}_{el} \to D`
- 2D plane-strain half-disc with mirror symmetry

**Results -- fields:**

.. figure:: images/full_cylinder_cracking/temperature_field.png
   :width: 48%

.. figure:: images/full_cylinder_cracking/stress_hoop_field.png
   :width: 48%

   Cold-contact wedge cools the rim (left); the tensile hoop-stress ring it sets
   up drives the cracking (right).

**Results -- simulation against experiment:**

.. figure:: images/full_cylinder_cracking/damage_field.png
   :width: 48%

.. figure:: images/full_cylinder_cracking/UO2_damage_sample.png
   :width: 48%

   Simulated damage with discrete radial cracks (left) against a cross-section of
   a real cracked UO\ :sub:`2` pellet (right).

.. figure:: images/full_cylinder_cracking/thermal_shock_results.png
   :width: 90%
   :align: center

   Quantitative verification: radial temperature profile, temperature history at
   the contact rim, and damage penetration, against McClenny Fig. 7b.

Numerical analysis
------------------

Mesh Convergence Study
^^^^^^^^^^^^^^^^^^^^^^^

**Case Directory:** ``cases/studies/mesh_sensitivity_2D``

This example demonstrates **mesh convergence** for a 2D thermal problem, showing how solution accuracy improves with mesh refinement.

**Setup:**

- 2D rectangular domain with prescribed boundary conditions
- Systematic refinement: :math:`n_x = 10, 20, 40, 80` elements
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
