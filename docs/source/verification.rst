Verification and testing
========================

This page lists what the Z3ST cases check, against which reference, and with
which error. All case paths are relative to ``z3st/cases/``.

Regression and verification
---------------------------

Two different checks are run on a case, and they answer different questions.

**Regression against a gold.** A case that carries
``output/non-regression_gold.json`` compares the metrics of the current run
with the values stored in that file. The comparison detects a change in the
results between two versions of the code. It does not establish that either
version is correct: a gold blessed from a wrong result protects the wrong
result. An example is recorded in the software paper. The rod case
``regression/pwr_rod_2D`` generated 223.7 W instead of the nominal 200 W,
because the radial power shape was normalised by the mean of its nodal values
rather than by its integral. Its regression test passed throughout, because the
burnup it checked was the nodal mean, which that normalisation preserves. The
defect was found by a comparison with OFFBEAT and is corrected in version
0.4.0.

**Verification against a reference.** A case whose metric has an analytical
solution, or an independent numerical solution of the same problem, compares
the run with that reference within a stated tolerance. This is the check that
says whether the result is right, to the accuracy of the reference and of the
discretisation.

Each ``non-regression.py`` writes ``output/non-regression.json`` with two
verdicts: ``summary`` (against the reference, where the case has one) and
``regression`` (against the gold). A metric with no closed form is written with
``rel_error = 0``, so it is protected by the gold comparison only.

Verification cases
------------------

The cases below have an analytical or independent reference. The values are
those of Table 6 of the software paper and agree with the
``output/non-regression_gold.json`` of each case.

.. list-table::
   :header-rows: 1
   :widths: 38 6 18 22 10

   * - Case
     - Dim.
     - Physics
     - Reference
     - Error
   * - ``verification/mechanics/lame_plane_strain_2D``
     - ax.
     - elasticity
     - Lamé, plane strain
     - 2.7e-5
   * - ``verification/mechanics/lame_gps_3D``
     - 3D
     - elasticity
     - Lamé, generalised plane strain
     - 2.2e-1
   * - ``verification/mechanics/mariotte_thin_shell``
     - ax.
     - elasticity
     - Mariotte membrane limit
     - 2.1e-2
   * - ``verification/thermal/spherical_shell``
     - 3D
     - gamma heating, thermo-elasticity
     - spherical shell, closed form
     - 6.0e-3
   * - ``verification/thermal/thick_cylindrical_shell_non_adiabatic_2D``
     - ax.
     - thermo-elasticity
     - hollow cylinder, analytical
     - 2.8e-2
   * - ``verification/mechanics/thermal_gradient_3D``
     - 3D
     - thermo-elasticity
     - imposed gradient, analytical
     - 3.5e-2
   * - ``verification/plasticity/j2_hardening_2D``
     - 2D
     - J2 plasticity
     - linear hardening, closed form
     - 2.6e-2
   * - ``verification/fuel/creep``
     - ax.
     - Norton creep
     - closed form
     - 1.2e-14
   * - ``verification/fuel/creep_irradiation``
     - ax.
     - irradiation creep
     - closed form
     - 8.9e-14
   * - ``verification/fuel/creep_shrink_fit_2D``
     - ax.
     - contact, creep relaxation
     - radial J2 solution
     - 1.9e-2
   * - ``verification/fuel/shrink_fit_disk_3d``
     - 3D
     - contact
     - Lamé interference pressure
     - 1.2e-2
   * - ``verification/thermal/coaxial_gap_3D``
     - 3D
     - gap conductance
     - analytical
     - 9.0e-6
   * - ``verification/fuel/burnup``
     - ax.
     - burnup accumulation
     - closed form
     - 1.1e-7

"ax." is the axisymmetric regime. The error is the relative difference against
the reference, taken as the largest over the stress components where several
are checked. The exceptions are:

- ``coaxial_gap_3D``: the maximum temperature error (``Linf_error_T``).
- ``creep`` and ``creep_irradiation``: the creep strain.
- ``burnup``: the mean burnup.
- ``creep_shrink_fit_2D`` and ``shrink_fit_disk_3d``: the contact pressure. For
  ``creep_shrink_fit_2D`` the value is the largest of the three relaxation
  times, at 600 days.
- ``j2_hardening_2D``: the largest :math:`\sigma_{xx}` error over the 21 load
  steps.
- ``lame_gps_3D``: the error is dominated by the axial stress of the
  generalised plane-strain constraint. The radial and hoop errors are
  1.1e-2 and 2.6e-3.
- ``mariotte_thin_shell``: the axial stress, the larger of the axial and hoop
  errors (the hoop error is 1.2e-2). The radial stress vanishes in the
  membrane limit and its relative error (0.58 in the gold) is not meaningful.

Mesh convergence
----------------

**Two-dimensional thermo-elastic slab** (``studies/mesh_sensitivity_2D``). A
slab 0.1 m wide and 1 m high with bilinear (Q1) elements, Dirichlet
temperatures of 490 and 480 K on its two faces and a volumetric gamma source.
``mesh_sensitivity.py`` regenerates the mesh with 10, 20, 40 and 80 cells across
the width, :math:`h = 10^{-2}` to :math:`1.25\times10^{-3}` m, with 40 cells
over the height in every run. The temperature and the stress
:math:`\sigma_{yy}` on the horizontal cut at mid-height are compared with their
analytical profiles in a discrete :math:`L_2` norm, normalised by the imposed
temperature difference and by the peak stress. The errors are written to
``convergence_data.txt``:

.. list-table::
   :header-rows: 1

   * - :math:`h` (m)
     - error in :math:`T`
     - error in :math:`\sigma_{yy}`
   * - 1.0e-2
     - 6.1e-3
     - 1.05e-2
   * - 5.0e-3
     - 1.55e-3
     - 2.43e-3
   * - 2.5e-3
     - 3.9e-4
     - 5.5e-4
   * - 1.25e-3
     - 9.9e-5
     - 2.7e-4

The temperature error falls by ratios of 3.93, 3.96 and 3.96, so the observed
order is two. The stress, compared at the cell centres, falls by 4.32, 4.39
and 2.08. It converges at second order down to 5.5e-4. On the finest mesh the
remaining error lies in the vertical direction, which the study does not
refine.

**Three-dimensional spherical shell** (``verification/thermal/spherical_shell``).
A steel shell of radii 2.0 and 2.5 m heated by gamma radiation entering through
the inner surface, with surfaces held at 494 and 491 K, the inner surface free
and the outer one clamped. ``exact.py`` evaluates the closed-form temperature,
displacement and stresses from the case input files. One octant is meshed with
structured hexahedra. ``mesh_convergence.py`` runs three meshes and writes
``convergence_data.txt``, and ``plots.py`` draws the profiles and the
convergence plot. Relative :math:`L_2` errors:

.. list-table::
   :header-rows: 1

   * - cells
     - :math:`T`
     - :math:`\sigma_{rr}`
     - :math:`\sigma_{\theta\theta}`
     - :math:`u_r`
   * - 480
     - 1.97e-3
     - 1.63e-2
     - 2.16e-2
     - 1.82e-3
   * - 3840
     - 5.10e-4
     - 4.25e-3
     - 5.98e-3
     - 4.14e-4
   * - 30720
     - 1.29e-4
     - 1.13e-3
     - 1.61e-3
     - 9.52e-5

Each refinement divides the temperature error by 3.87 and 3.94, the stress
errors by 3.61 to 3.83 and the displacement error by 4.41 and 4.35, so every
quantity converges at second order. The gold is written on the 3840-cell mesh.
There is no model difference between Z3ST and this reference, and the error
goes to zero with the mesh.

**Three-dimensional contact** (``verification/fuel/shrink_fit_disk_3d``). A
quarter of a solid disc of radius 4.1 mm inside a ring of radii 4.13 and
4.75 mm, 1 mm thick, with an initial gap of 30 µm, closed by uniform heating
of the disc and loaded through penalty contact. The reference is the
plane-stress Lamé interference pressure. At the final step the computed
pressure is 74.19 MPa against the analytical 74.97 MPa, a difference of 1.2 %
in the gold (the paper quotes 1.0 % at the final step and 1.2 % at the first
closed step, the largest). The software paper reports three meshes, with
21606, 49131 and 102228 displacement degrees of freedom, and errors of 1.39,
1.16 and 1.10 %, approaching a floor near 1.05 %. The floor is a model
difference: the penalty contact admits a penetration of order
:math:`p/k_\mathrm{pen}`, and the reference is a plane-stress solution while the
computed body has a symmetry plane and a free surface. The case is therefore
accurate to about 1 % against its reference, and is not claimed to converge to
it. The three-mesh series is not stored in the repository.

**Porosity migration** (``verification/fuel/porosity_migration``). The software
paper reports a radial mesh refinement from 100 to 200 and 400 elements, which
leaves the centre temperature at 2949.1, 2949.0 and 2949.0 K and the void
radius at 0.1600, 0.1600 and 0.1625 :math:`R_o`. The refined meshes are not
stored in the repository.

The phase-field case ``benchmarks/damage/pellet_quench_2D_xy`` resolves the
regularisation length with four elements and is not refined further.

Independent reference for creep relaxation
------------------------------------------

``verification/fuel/creep_shrink_fit_2D`` is an axisymmetric UO\ :sub:`2`
pellet shrink-fitted into a Zircaloy-4 cladding, whose contact pressure relaxes
by Norton creep of the cladding. The reference is ``reference_1d.py``, a radial
solution of the same joint (elastic solid shaft, J2 hub in generalised plane
strain, penalty spring in series, explicit creep update converged in time),
which shares no code with Z3ST. The case asserts four values:

.. list-table::
   :header-rows: 1

   * - t (days)
     - Z3ST (MPa)
     - ``reference_1d.py`` (MPa)
     - deviation
   * - 0 (elastic)
     - 24.703
     - 24.699
     - 1.7e-4
   * - 600
     - 7.256
     - 7.120
     - 1.9 %
   * - 1240
     - 5.140
     - 5.063
     - 1.5 %
   * - 2500
     - 3.648
     - 3.605
     - 1.2 %

The remaining deviation is an error of the time step: doubling ``n_steps``
takes the 2500-day value from 3.648 to 3.633 MPa. The closed form of Esposito
et al. (*Int. J. Pressure Vessels and Piping* 185 (2020) 104126) falls about
3 % below the radial reference, because it uses the Tresca criterion in the
hub. It is plotted by ``plots.py`` and not asserted. Two corrections are
applied to that closed form in the case, and are documented in the case
``README.md``:

- Eq. (4) of the paper carries the Poisson terms of shaft and hub with swapped
  signs. They cancel for equal materials, as in the paper's own validation,
  and raise the elastic factor by 3 % for UO\ :sub:`2` and Zircaloy.
  ``case_params.elastic_factor`` carries the Lamé signs.
- The hub factor :math:`k_2` is taken from Eq. (20) of the paper.

Guard cases
-----------

``verification/thermal/two_heated_materials_2D`` guards a defect of the volumetric
heat source. The source is nodal, and a node on the interface between two
heated materials belongs to both. Before version 0.4.0 the source of each
material was added at those nodes, so the interface nodes carried twice the
source. The case splits a gamma-heated cylindrical wall into two identical
materials and compares its temperature, node by node, with the same wall meshed
as one material (``single/``), with a tolerance of 1e-6.

The suites
----------

**Local suite.** ``z3st/cases/non-regression_local.sh`` discovers its cases:
every directory under ``z3st/cases/`` that contains both an ``Allrun`` and an
``output/non-regression_gold.json`` is a member, ``sandbox/`` is never
scanned, and the cases listed in ``suite_exclude.txt`` are skipped. A case
fails if ``Allrun`` exits non-zero, if ``output/non-regression.json`` is
missing, or if its ``summary`` or ``regression`` verdict is not ``PASS``.

.. code-block:: bash

   cd z3st/cases
   ./non-regression_local.sh --list    # discovered cases and exclusions
   ./non-regression_local.sh           # run all
   ./non-regression_local.sh verification/fuel/creep   # run named cases

68 case directories carry a gold. Five are listed in ``suite_exclude.txt``
(``benchmarks/damage/sen_shear``, ``benchmarks/damage/sen_tension``,
``benchmarks/damage/plate_thermal_shock_2D``,
``verification/fuel/creep_law_discovery`` and ``regression/pwr_rod_2D``), so
the local suite has 63 cases, 17 of them three-dimensional. The excluded cases
are run by hand. ``regression/pwr_rod_2D`` carries a gold and is run by hand
outside the suite, because the contact path it covers is also covered by
``creep_shrink_fit_2D``.

**Continuous integration.** ``.github/workflows/ci.yml`` runs on every push to
``main``, ``develop`` and ``development/*`` and on every pull request, inside
the ``dolfinx/dolfinx:v0.11.0`` container. It calls
``z3st/cases/non-regression_github.sh``, which runs the 21 cases listed in
``z3st/cases/cases_ci.txt``. A case enters that list if it has a stable gold and
exercises a physics path no other listed case reaches. Every case in both
suites runs serially, so neither suite tests parallel execution.

**Blessing a gold.** A gold is created or updated by copying the verdict file
of a run:

.. code-block:: bash

   cp output/non-regression.json output/non-regression_gold.json

Do this only after the run has been reviewed against its reference and, for a
changed gold, after the cause of the change is understood. Blessing a gold
makes the case a member of the local suite.

**Static checks.** ``python -m z3st.utils.audit_checks`` runs checks on the
tree that need no simulation (``--list`` describes each one). Among them it
reports physics models that no gold-carrying case reaches, golds that contain
NaN or inf, and case paths cited in the documentation that do not exist.
