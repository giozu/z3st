Physics Models
==============

Z3ST provides a set of physical models for multiphysics thermo-mechanical
analysis:

- **Thermal conduction** with volumetric heating, temperature-dependent
  properties, and gap conductance between bodies
- **Mechanical equilibrium** with five constitutive routes: isotropic
  (Lamé) and anisotropic (Voigt) linear elasticity, Neo-Hookean
  hyperelasticity, J2 plasticity, and a custom Python route
- **J2 and crystal plasticity** for the inelastic response
- **Phase-field fracture** (AT1 / AT2) for crack initiation and propagation
- **Gap conductance** for heat transfer across the interface between
  separate bodies
- **Cluster dynamics** for defect-cluster evolution in size space

All models are written as FEniCSx (UFL) variational forms. A central feature
is that the constitutive physics is expressed symbolically and differentiated
**automatically**: the first Piola--Kirchhoff stress is obtained as
:math:`\boldsymbol{P} = \partial \psi / \partial \boldsymbol{F}` and the
consistent Newton tangent as the UFL derivative of the residual, so no
constitutive Jacobian is hand-coded. The models are modular and coupled
through a staggered solution scheme with adaptive relaxation.

.. contents:: On this page
   :local:
   :depth: 1

----

Thermal Model
-------------

The thermal model solves the heat-conduction equation

.. math::

   \rho C_p \frac{\partial T}{\partial t} - \nabla \cdot (k \nabla T) = q''',

in steady-state or transient mode (selected by configuration), with backward
Euler time integration in the transient case.

- :math:`T` -- temperature (K)
- :math:`\rho` -- density (kg/m³), :math:`C_p` -- specific heat (J/(kg·K))
- :math:`k` -- thermal conductivity (W/(m·K)), constant or a user-supplied
  Python function of temperature
- :math:`q'''` -- volumetric heat source (W/m³)

**Boundary conditions** are applied on a disjoint partition
:math:`\partial\Omega = \partial\Omega_T \cup \partial\Omega_q \cup \partial\Omega_R`:

- Dirichlet: :math:`T = T_d`
- Neumann: :math:`-k \nabla T \cdot \boldsymbol{n} = q_n`
- Robin: :math:`-k \nabla T \cdot \boldsymbol{n} = h\,(T - T_{ext})`, in either a
  convective variant (user film coefficient :math:`h`) or a gap-coupled
  variant in which :math:`h` comes from the :ref:`gap-conductance model
  <gap-conductance>`.

**Volumetric source** ``q'''`` supports a fissile mode
(:math:`q''' = \mathrm{LHR}/A`, linear heat rate divided by the pellet area),
an exponential :math:`\gamma`-heating decay :math:`q''' = q_0 \exp(-\mu_\gamma r)`
in rectangular geometry, a modified-Bessel decay
:math:`K_0(\mu_\gamma r)/K_0(\mu_\gamma R_i)` in cylindrical geometry, and a
spherical decay.

**Weak form.** Find :math:`T \in V_T` such that for all :math:`v \in V_T`

.. math::

   \int_\Omega \rho C_p \frac{\partial T}{\partial t} v \, \mathrm{d}\Omega
   + \int_\Omega k \nabla T \cdot \nabla v \, \mathrm{d}\Omega
   + \int_{\Gamma_R} h\, T v \, \mathrm{d}s
   = \int_\Omega q''' v \, \mathrm{d}\Omega
   + \int_{\Gamma_N} q_n v \, \mathrm{d}s
   + \int_{\Gamma_R} h\, T_{ext}\, v \, \mathrm{d}s .

Implemented in :class:`z3st.models.thermal_model.ThermalModel`.

.. figure:: images/thin_slab/stress_temperature_combined.png
   :width: 80%
   :align: center

   Coupled thermo-mechanical thin slab: temperature and thermal stress through
   the thickness, numerical (markers) against the analytical solution (lines).

----

Mechanical Model
----------------

Quasi-static mechanical equilibrium is the linear-momentum balance

.. math::

   -\nabla \cdot \boldsymbol{\sigma} = \boldsymbol{f},

with body force :math:`\boldsymbol{f}` (e.g. gravity along the regime's vertical
axis) and the small-strain tensor
:math:`\boldsymbol{\varepsilon} = \tfrac{1}{2}(\nabla \boldsymbol{u} + \nabla \boldsymbol{u}^\top)`.
The out-of-plane components follow the active regime: retained (plane strain),
reconstructed from the in-plane solution (plane stress), solved as part of
:math:`\boldsymbol{u}` (3D), or expressed in cylindrical components
(axisymmetric).

**Weak form.** Find :math:`\boldsymbol{u} \in V_u` such that for all
:math:`\boldsymbol{v} \in V_u`

.. math::

   \int_\Omega \boldsymbol{\sigma}(\boldsymbol{u}) : \nabla \boldsymbol{v} \, \mathrm{d}\Omega
   = \int_\Omega \boldsymbol{f} \cdot \boldsymbol{v} \, \mathrm{d}\Omega
   + \int_{\Gamma_N} \boldsymbol{t}_0 \cdot \boldsymbol{v} \, \mathrm{d}s .

Five constitutive routes are selected by the ``constitutive`` field of each
material card.

Isotropic linear elasticity (Lamé)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The default route is

.. math::

   \boldsymbol{\sigma} = \mathbb{C} : \big( \boldsymbol{\varepsilon} - \alpha (T - T_{ref}) \boldsymbol{I} \big),

with :math:`\mathbb{C}` built from the Lamé parameters :math:`\lambda, G`. The
thermal-stress contribution
:math:`\boldsymbol{\sigma}_{th} = -(3\lambda + 2G)\,\alpha (T - T_{ref})\,\boldsymbol{I}`
is moved to the right-hand side of the momentum balance.

.. note::

   **Consistent damage degradation of the thermal stress.** When the phase-field
   damage block is active, *both* the elastic and the thermal-stress
   contributions are weighted by the same degradation function :math:`g(D)`, so
   the effective stress driving equilibrium is

   .. math::

      \boldsymbol{\sigma}_{eff} = g(D)\, \mathbb{C} : \big( \boldsymbol{\varepsilon} - \alpha (T - T_{ref}) \boldsymbol{I} \big).

   A fully damaged cell (:math:`D \to 1`, :math:`g \to K`) then recovers the
   traction-free crack-face limit for both driving forces. Degrading only the
   elastic term while leaving :math:`\boldsymbol{\sigma}_{th}` at full magnitude
   would force cracked cells to carry an unopposed thermal body force against a
   stiffness reduced to :math:`g(D)\to K=10^{-6}`, producing strains that grow as
   :math:`1/K` and a catastrophic loss of energy conservation. The consistent
   :math:`g(D)` weighting of both contributions is therefore essential to any
   thermo-mechanical-fracture problem in which a real crack opens.

Anisotropic elasticity (Voigt)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A user-supplied :math:`6\times 6` elasticity matrix in Voigt form,
:math:`\boldsymbol{\sigma} = \mathbb{C}_{\text{Voigt}} : \boldsymbol{\varepsilon}`,
makes anisotropic elasticity available without changing the assembly code.

Neo-Hookean hyperelasticity
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. note::

   This route is **implemented and verified** (case
   ``1_thin_slab_non_linear``), not planned.

Compressible Neo-Hookean hyperelasticity uses the strain-energy density

.. math::

   \psi(\boldsymbol{F}) = \frac{\mu}{2}(I_C - 3) - \mu \ln J + \frac{\lambda}{2} (\ln J)^2,

with :math:`\boldsymbol{C} = \boldsymbol{F}^\top \boldsymbol{F}`,
:math:`I_C = \mathrm{tr}\,\boldsymbol{C}`, :math:`J = \det \boldsymbol{F}`, and
:math:`\boldsymbol{F} = \boldsymbol{I} + \nabla \boldsymbol{u}`. The first
Piola--Kirchhoff and Cauchy stresses are obtained by **automatic
differentiation** of the energy,

.. math::

   \boldsymbol{P} = \frac{\partial \psi}{\partial \boldsymbol{F}}, \qquad
   \boldsymbol{\sigma} = J^{-1} \boldsymbol{P} \boldsymbol{F}^\top ,

so a single line, ``P = ufl.diff(psi, F)``, replaces a hand-derived stress
tensor; the consistent tangent is the UFL derivative of the residual. The
problem is solved with the PETSc SNES Newton driver and a configurable line
search.

Custom constitutive route
^^^^^^^^^^^^^^^^^^^^^^^^^^^

A user-supplied stress function loaded from an arbitrary Python module gives an
explicit extension point for arbitrary constitutive laws; this is the route used
by the crystal-plasticity demonstration below.

**Boundary conditions.** Vector and per-component Dirichlet conditions
(:math:`u_x = u_x^d`, ...), :math:`\mathrm{Clamp}_i` and :math:`\mathrm{Slip}_i`
variants that constrain a single component to zero, and scalar Neumann tractions
along the facet normal. Any value may be supplied as a list of length
``n_steps`` for step-wise loading.

Implemented in :class:`z3st.models.mechanical_model.MechanicalModel`; supports
2D (plane stress / plane strain), 3D, and axisymmetric regimes.

.. figure:: images/cylindrical_shell/stress_comparison.png
   :width: 75%
   :align: center

   Linear elasticity: radial and hoop stress in a thick cylindrical shell,
   numerical against the analytical Lamé solution.

J2 plasticity (von Mises)
^^^^^^^^^^^^^^^^^^^^^^^^^^

A small-strain J2 model with linear isotropic hardening, integrated by a
return-mapping update at quadrature points. The yield function is

.. math::

   f(\boldsymbol{\sigma}, p) = \sigma_{eq} - \sigma_y(p), \qquad
   \sigma_y(p) = \sigma_0 + H p,

with the von Mises equivalent stress
:math:`\sigma_{eq} = \sqrt{\tfrac{3}{2}\,\boldsymbol{s}:\boldsymbol{s}}`,
deviatoric stress :math:`\boldsymbol{s}`, hardening modulus :math:`H`, and
cumulative plastic strain :math:`p`. An elastic predictor is computed each
increment; if the yield surface is violated, a radial-return correction is
applied,

.. math::

   \Delta p = \frac{\langle f_{trial} \rangle_+}{3G + H}, \qquad
   \boldsymbol{\sigma} = \boldsymbol{\sigma}_{trial} - 3G \Delta p\, \boldsymbol{n}, \qquad
   \boldsymbol{n} = \frac{\boldsymbol{s}_{trial}}{\sigma_{eq,\,trial}} .

The plastic-strain tensor and :math:`p` live on quadrature spaces and are
updated at the end of each converged staggered iteration. Implemented in
:class:`z3st.models.plasticity_model`.

.. figure:: images/plasticity_2D/output/stress_strain_curve.png
   :width: 65%
   :align: center

   J2 plasticity: stress-strain response (numerical vs analytical, plane strain),
   recovering the elastic slope and the linear hardening branch.

Crystal plasticity (single grain)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``plasticity.mode: custom`` hook replaces the J2 update with a user-supplied
routine. The demonstration case ``demo_CP_single_grain`` implements
rate-dependent single-crystal plasticity on one 3D grain. The stress follows the
elastic law on the plastic-corrected strain,

.. math::

   \boldsymbol{\sigma} = \mathbb{C} : (\boldsymbol{\varepsilon} - \boldsymbol{\varepsilon}^p),

with the plastic strain accumulated over the active slip systems through the
Schmid tensor and a power-law slip rate,

.. math::

   \dot{\boldsymbol{\varepsilon}}^p = \sum_s \dot{\gamma}^s\, \boldsymbol{P}^s, \qquad
   \boldsymbol{P}^s = \tfrac{1}{2}\big(\boldsymbol{m}^s \otimes \boldsymbol{n}^s + \boldsymbol{n}^s \otimes \boldsymbol{m}^s\big),

.. math::

   \dot{\gamma}^s = \dot{\gamma}_0 \left| \frac{\tau^s}{g_0} \right|^{n} \operatorname{sign}(\tau^s), \qquad
   \tau^s = \boldsymbol{\sigma} : \boldsymbol{P}^s ,

where :math:`\boldsymbol{m}^s` and :math:`\boldsymbol{n}^s` are the slip
direction and slip-plane normal, :math:`\tau^s` is the resolved shear stress,
:math:`g_0` is the slip resistance (CRSS), :math:`\dot{\gamma}_0` is the
reference slip rate, and :math:`n` is the power-law exponent. Time integration is
backward Euler,
:math:`\boldsymbol{\varepsilon}^p_{n+1} = \boldsymbol{\varepsilon}^p_{n} + \Delta t\, \dot{\boldsymbol{\varepsilon}}^p_{n+1}`,
with the plastic strain stored as a history variable. The demo uses an FCC
system :math:`(111)[01\bar{1}]` (Schmid factor :math:`\approx 0.408`) and obtains
the **exact Jacobian by automatic differentiation** -- precisely the case that is
painful to differentiate by hand. The result is verified against saturation-stress
theory.

.. figure:: images/demo_CP_single_grain/output/stress_strain_curve.png
   :width: 65%
   :align: center

   Crystal plasticity (single grain): the Z3ST response saturates towards the
   analytical saturation stress :math:`\sigma_{sat}`.

----

Phase-Field Damage (Fracture Mechanics)
---------------------------------------

Brittle fracture is modelled with the variational phase-field approach, in which
the sharp crack is regularised by a continuous damage field
:math:`D \in [0, 1]`. The total energy is

.. math::

   E_{tot}(\boldsymbol{u}, D) = \int_\Omega g(D)\, \psi^+(\boldsymbol{\varepsilon}) \, \mathrm{d}\Omega
   + \int_\Omega \psi^-(\boldsymbol{\varepsilon}) \, \mathrm{d}\Omega
   + G_c \int_\Omega \gamma(D, \nabla D) \, \mathrm{d}\Omega,

with the quadratic degradation function :math:`g(D) = (1-D)^2 + K` and
regularisation :math:`K = 10^{-6}`. Two crack-density functionals are supported,

.. math::

   \gamma_{AT2} = \frac{1}{2}\left( \frac{D^2}{\ell_c} + \ell_c\, |\nabla D|^2 \right), \qquad
   \gamma_{AT1} = \frac{3}{8}\left( \frac{D}{\ell_c} + \ell_c\, |\nabla D|^2 \right),

where :math:`\ell_c` is the regularisation length. AT1 has an analytical elastic
threshold below which damage cannot grow; AT2 has none.

**Energy splits.** Two splits define the crack-driving force :math:`\psi^+`: the
Miehe spectral split (positive principal strains drive cracking) and the Amor
volumetric/deviatoric split,

.. math::

   \psi^+ = \tfrac{1}{2}\lambda \langle \mathrm{tr}\,\boldsymbol{\varepsilon} \rangle_+^2
   + G\, \mathrm{dev}(\boldsymbol{\varepsilon}) : \mathrm{dev}(\boldsymbol{\varepsilon}),

which is cheaper and well suited to AT1.

**Staggered governing equations.** With fixed :math:`D`, solve mechanics with the
degraded stress; with fixed :math:`\boldsymbol{u}`, solve the damage sub-problem
driven by the history field :math:`H`. The history is stored cell-by-cell on the
:math:`DG_0` space: non-dimensional for AT2 (:math:`H = (2\ell_c/G_c)\,\psi^+`)
and as the physical positive elastic energy density for AT1.

**Irreversibility** is enforced point-wise after every staggered iteration,

.. math::

   D^{n+1} = \min\!\big(1, \max(D^{n+1}, D^{n})\big).

**Hybrid constraint (optional).** The Ambati-Gerasimov-De Lorenzis hybrid
constraint sets the :math:`\psi^+` contribution to zero in cells where
:math:`\psi^- > \psi^+` (predominantly compression), suppressing spurious crack
growth while preserving the energy that resists fracture.

**Thermo-mechanical coupling.** When thermal, mechanical, and damage blocks are
all active, :math:`\psi^+` is evaluated on the **elastic** strain
:math:`\boldsymbol{\varepsilon}_{el} = \boldsymbol{\varepsilon} - \alpha(T - T_{ref})\boldsymbol{I}`,
so that unconstrained uniform thermal expansion produces zero releasable energy
and zero crack-driving force. In ``regime: 2d`` plane strain, the z-component of
the thermal eigenstrain is suppressed in the damage driver only (the geometrically
blocked z-expansion would otherwise drive damage everywhere through the deviatoric
channel); the mechanical equilibrium still uses the full plane-strain thermal
stress.

**Pre-crack seeding.** Damage Dirichlet conditions :math:`D = D_d` fix healthy
regions; imposing :math:`D = 1` on a thin internal slit seeds a pre-crack and lets
the staggered scheme propagate damage from a known position rather than relying on
spontaneous nucleation.

Implemented in :class:`z3st.models.damage_model.DamageModel`.

.. figure:: images/sen_shear/SENS_damage_final.png
   :width: 55%
   :align: center

   Single-edge-notched shear test (steel): the curved crack path reproduces the
   benchmark of Miehe et al., *Comput. Methods Appl. Mech. Engrg.* 199 (2010).

----

.. _gap-conductance:

Gap Conductance Model
---------------------

Heat transfer between two paired bodies separated by a small gap is a Robin
coupling with an effective film coefficient :math:`h_{gap}`, taken either from a
fixed user value or from a gas-conduction correlation,

.. math::

   k_{gas} = c \cdot 10^{-4}\, T_{gap}^{0.79}, \qquad
   h_{gap} = \frac{k_{gas}}{d_{gap}},

where :math:`T_{gap} = \tfrac{1}{2}(T_{inner} + T_{outer})` is the mean surface
temperature, :math:`c` a user prefactor, and :math:`d_{gap}` the mean
centroid-to-centroid distance between the two paired facet groups, computed once
per simulation with a SciPy cKDTree query. The gap pair is declared through the
``pair`` field of a Robin boundary condition in ``boundary_conditions.yaml``; no
specialised contact element is required. This is the capability that makes Z3ST
**multi-body** -- for example, heat transfer from a fuel pellet to its cladding.
Implemented in :class:`z3st.models.gap_model.GapModel`.

----

Cluster Dynamics Model
----------------------

A standalone one-dimensional advection-diffusion solver in cluster-size space
:math:`n` for defect-cluster size distributions :math:`c(n, t)`,

.. math::

   \frac{\partial c}{\partial t} + v\, \frac{\partial c}{\partial n} - D\, \frac{\partial^2 c}{\partial n^2} = 0,

with user-selectable initial conditions (constant on a labelled region, or a
Gaussian of prescribed mean, width, and amplitude). The discretisation uses a
:math:`DG_1` space with an upwind interior-facet flux for the advective term and
a symmetric interior-penalty Galerkin (SIPG) treatment of the diffusive term;
time integration is implicit Euler. The total cluster mass
:math:`\int c\, n \, \mathrm{d}n` is rescaled to its initial value at every step
to enforce conservation, and the local Péclet number is logged for diagnostics.

.. note::

   Cluster dynamics is a newer, exploratory capability -- a path towards a
   genuine micro-to-continuum link -- and is not yet the subject of a published
   verification study.

Implemented in :class:`z3st.models.cluster_dynamic_model.ClusterDynamicModel`.

----

Coupled Thermo-Mechanical Analysis
----------------------------------

The default coupling between thermal, mechanical, damage, and cluster physics is
a **staggered scheme**. Within each time step the active blocks are solved in
sequence -- thermal, then mechanical with the updated temperature, then damage
with the updated displacement (after refreshing the history field), then cluster
-- and each new field is under-relaxed,

.. math::

   \phi \leftarrow \alpha\, \phi^{new} + (1 - \alpha)\, \phi^{old},

with separate factors :math:`\alpha_T, \alpha_u, \alpha_D`. The factors are
**adapted automatically**: an exponential moving average of the relative residual
is tracked, and each :math:`\alpha` is grown when convergence is monotone or
shrunk on divergence, within hard bounds
:math:`\alpha \in [\alpha_{min}, \alpha_{max}]`. The staggered loop converges when
the relative or absolute change in every active field drops below tolerance
simultaneously,

.. math::

   \|T^{k+1} - T^k\| < \mathrm{tol}_T, \qquad
   \|\boldsymbol{u}^{k+1} - \boldsymbol{u}^k\| < \mathrm{tol}_u, \qquad
   \|D^{k+1} - D^k\| < \mathrm{tol}_D .

**Solvers.** Each block is a separate variational problem solved through PETSc: a
direct LU solver (MUMPS), or CG/GMRES with smoothed-aggregation AMG (PETSc GAMG)
or HYPRE BoomerAMG. Symmetric positive-definite blocks (thermal, damage) use
conjugate gradients; the mechanical block uses GMRES. Non-linear mechanical
problems (hyperelasticity, custom non-linear route) use PETSc SNES with a
Newton line search.

----

Application: thermal-shock cracking of a UO\ :sub:`2` pellet
------------------------------------------------------------

The full coupled set -- thermal, mechanical, and phase-field damage -- meets in
the UO\ :sub:`2` thermal-shock case (``14_full_cylinder_cracking_2D_xy``), a 2D
plane-strain transverse cross-section reproducer of McClenny et al.,
*J. Nucl. Mater.* 565 (2022). A cold-contact wedge cools the rim of a hot disc;
the tensile hoop-stress ring it sets up drives discrete radial cracks (AT1 +
Amor + hybrid, fully coupled :math:`T \to \boldsymbol{\varepsilon}_{el} \to D`).

.. figure:: images/full_cylinder_cracking/temperature_field.png
   :width: 49%

.. figure:: images/full_cylinder_cracking/stress_hoop_field.png
   :width: 49%

   Cold-contact wedge cools the rim (left); the tensile hoop-stress ring it sets
   up drives the cracking (right).

.. figure:: images/full_cylinder_cracking/damage_field.png
   :width: 49%

.. figure:: images/full_cylinder_cracking/UO2_damage_sample.png
   :width: 49%

   Simulated damage with discrete radial cracks at the rim (left) against a
   cross-section of a real cracked UO\ :sub:`2` pellet (right).

.. figure:: images/full_cylinder_cracking/thermal_shock_results.png
   :width: 90%
   :align: center

   Quantitative verification: radial temperature profile, temperature history at
   the contact rim, and damage penetration, against McClenny Fig. 7b.

----

**See also**

- :doc:`examples` -- the full verification-case catalogue
- :doc:`differentiable_features` -- automatic differentiation in Z3ST
- :doc:`usage` -- YAML configuration of every model
- :doc:`api` -- implementation reference
