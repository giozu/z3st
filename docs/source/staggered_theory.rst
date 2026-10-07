.. _staggered-theory:

The Staggered Solver: What It Computes
======================================

Z3ST solves temperature, displacement, phase field, cluster distribution and
porosity one block at a time within each time step and repeats the sequence until every block passes
its staggered tolerance (:ref:`coupled-scheme`). This page states what the
converged result is, which theoretical guarantees apply to the formulation
implemented, and what has been measured in their place.

.. contents:: On this page
   :local:
   :depth: 1


The fixed point
---------------

Write the discrete equations of one time step as block residuals
:math:`R_T(T, \boldsymbol u, d) = 0`, :math:`R_u(T, \boldsymbol u, d) = 0` and
:math:`R_d(T, \boldsymbol u, d) = 0`, the last followed by the projection
:math:`d \leftarrow \min(1, \max(d, d^n))` and with :math:`\mathcal H` updated
from the current displacement. One staggered iteration is a non-linear block
Gauss-Seidel sweep,

.. math::

   \tilde T^{k} : R_T(\tilde T^{k}, \boldsymbol u^{k-1}, d^{k-1}) = 0, \qquad
   \tilde{\boldsymbol u}^{k} : R_u(T^{k}, \tilde{\boldsymbol u}^{k}, d^{k-1}) = 0, \qquad
   \tilde d^{k} : R_d(T^{k}, \boldsymbol u^{k}, \tilde d^{k}) = 0,

each followed by the relaxation
:math:`X^{k} = \omega_X \tilde X^{k} + (1 - \omega_X) X^{k-1}`. The loop stops
when :math:`\|\tilde X^k - X^{k-1}\|_2/\|\tilde X^k\|_2` is below
``stag_tol`` for every field in the same iteration (or the absolute norm, with
``convergence: norm``).

At an exact fixed point, :math:`\tilde X^k = X^{k-1}` for every block, so every
block residual vanishes at the same state. That state solves the coupled
discrete system that a monolithic solver applied to the same equations would
solve. Z3ST has no monolithic route, so this equivalence is a property of the
equations and has not been checked numerically in the code. At a finite
tolerance the loop stops at a state where the increments are small, and the
residuals are not evaluated. The distance from the fixed point is therefore
measured by varying the tolerance (below).

Two consequences of the implementation:

- The convergence measure uses the unrelaxed solve output :math:`\tilde X^k`,
  so it does not shrink when the relaxation factor is small.
- A step that reaches ``max_iters`` is accepted with a warning unless adaptive
  time stepping is on, in which case it is bisected. An accepted
  non-converged step is not a fixed point.


Energy-based guarantees and the hybrid formulation
--------------------------------------------------

For the variational phase-field model, in which the elastic energy is
:math:`\int g(d)\,\psi(\boldsymbol\varepsilon)\,\mathrm{d}x` with no split and
the phase-field equation is the stationarity condition of the same energy
functional, the energy is convex in :math:`\boldsymbol u` at fixed :math:`d`
and convex in :math:`d` at fixed :math:`\boldsymbol u`. The staggered scheme
is then the alternate minimisation of Bourdin, Francfort and Marigo
[BourdinFrancfortMarigo2000]_: each block lowers the energy, and the iterates
converge monotonically to a stationary point [FarrellMaurini2017]_. For that
formulation the a-posteriori truncation :math:`d \leftarrow \max(d, d^n)` has
been analysed by Almi [Almi2020]_.

The default route of Z3ST is the hybrid formulation of Ambati et al.
[Ambati2015]_. The stress is degraded as
:math:`g(d)\,\mathbb C:\boldsymbol\varepsilon_{el}`, while the phase field is
driven by :math:`\psi^+` through a history field, with a split selected by
``damage.split`` and the constraint that suppresses growth where
:math:`\psi^- > \psi^+`. These equations are not the stationarity conditions
of one energy functional [Ambati2015]_. For the hybrid formulation, therefore:

- no energy-descent property of the staggered iteration holds, and
  convergence of the iteration is not guaranteed by the argument above,
- the elastic and fracture energies written to ``energies.txt`` are
  diagnostics, and their sum is not a conserved quantity.

The thermal block adds a further departure from the variational setting. The
temperature enters the mechanics through the eigenstrain and the mechanics does
not enter the heat equation, apart from the gap conductance and the contact
pressure. The relaxation also changes the iteration map, so even for a
variational model a relaxed sweep is not an exact block minimisation.

The relations used to convert between :math:`G_c` and :math:`\sigma_c`
(AT1: :math:`G_c = \tfrac83 \ell\sigma_c^2/E`, AT2:
:math:`G_c = \tfrac{256}{27}\ell\sigma_c^2/E`) are the critical stresses of the
homogeneous one-dimensional solution of the variational AT1 and AT2 models
[Pham2011]_ [Tanne2018]_. Under the hybrid formulation with a split, the
multiaxial strength is set by the split as well.


What is measured instead
------------------------

The paper accompanying Z3ST version 0.4 reports the following measurements.

**Tolerance sweeps.** The same tolerance applied to every field was varied from
:math:`10^{-3}` to :math:`10^{-6}`:

==================  =========  ================  ==============  =================
Case                Tolerance  Iterations/step   Error vs Lamé   Peak fuel T (K)
==================  =========  ================  ==============  =================
Shrink fit, 2D      1e-3       11.8              0.48 %          --
Shrink fit, 2D      1e-6       15.9              0.48 %          --
Shrink fit, 3D      1e-3       15.4              0.47 %          --
Shrink fit, 3D      1e-6       18.7              0.47 %          --
Rod, no creep       1e-3       12.7              --              1076.41
Rod, no creep       1e-6       24.1              --              1076.30
==================  =========  ================  ==============  =================

The error is against the Lamé interference pressure with the ring compliance
at the ring's inner radius. The error of the
contact pressure against the analytical Lamé interference
pressure does not change over the three decades, so the coupling error is
below 0.005 percentage points of that error. In the rod case, with the
cladding creep off, a 200-day horizon and adaptive time stepping off, the peak
fuel temperature spans 0.11 K and the hot gap 3 nm. Each decade of tolerance
costs on average 1.4, 1.1 and 3.8 more staggered iterations per step for the
three cases.

**Thermal shock.** In ``cases/benchmarks/damage/pellet_quench_2D_xy`` the
fracture energy is non-decreasing at every one of the 100 steps, which is the
check of irreversibility that this quantity supports. Over the transient the
elastic energy falls by 9.09 J and the fracture energy rises by 4.90 J. The
two terms are not a closed balance, because the thermal eigenstrain does work
on the body and the history field does not derive from an energy functional.
The case takes 3265 staggered iterations over its 100 steps, with a median of
25 per step and up to 159 during crack growth.


Relaxation and acceleration
---------------------------

- The **adaptive** relaxation (``relax_adaptive``) multiplies each factor by
  ``relax_growth`` when the residual is below its exponential moving average
  and by ``relax_shrink`` otherwise. It is a heuristic and is not Aitken's
  method.
- **Aitken's** :math:`\Delta^2` method (``relax_aitken``) acts on the
  displacement only, with a separate Aitken option for porosity
  (``porosity.aitken``).
- Both are clamped to [``relax_min``, ``relax_max``], default [0.05, 1.0].
  With the default upper bound no factor exceeds one, so the scheme never
  over-relaxes. The over-relaxed alternate minimisation of Farrell and
  Maurini [FarrellMaurini2017]_ is therefore not what Z3ST does by default.

Gerasimov and De Lorenzis [Gerasimov2016]_ show that the results of staggered
schemes for phase-field fracture depend on the staggered tolerance. The
tolerance sweeps above report that dependence for Z3ST.


Stability of the computed state
-------------------------------

A converged staggered iteration delivers a state at which the block equations
hold. For softening models such a state may be unstable, and first-order
solvers, staggered or monolithic, can converge to an unstable branch
[LeonBaldelliMaurini2021]_. The stability criterion of the variational theory
is the positivity of the second variation of the energy on the cone of
admissible directions [Pham2011]_. Z3ST does not evaluate it, and for the
hybrid formulation there is no energy whose second variation could be
evaluated.


References
----------

.. [Ambati2015] M. Ambati, T. Gerasimov, L. De Lorenzis, *A review on
   phase-field models of brittle fracture and a new fast hybrid formulation*,
   Comput. Mech. 55 (2015) 383--405. `doi:10.1007/s00466-014-1109-y
   <https://doi.org/10.1007/s00466-014-1109-y>`_

.. [BourdinFrancfortMarigo2000] B. Bourdin, G. A. Francfort, J.-J. Marigo,
   *Numerical experiments in revisited brittle fracture*, J. Mech. Phys. Solids
   48 (2000) 797--826. `doi:10.1016/S0022-5096(99)00028-9
   <https://doi.org/10.1016/S0022-5096(99)00028-9>`_

.. [FarrellMaurini2017] P. E. Farrell, C. Maurini, *Linear and nonlinear solvers
   for variational phase-field models of brittle fracture*, Int. J. Numer.
   Methods Eng. 109 (2017) 648--667. `doi:10.1002/nme.5300
   <https://doi.org/10.1002/nme.5300>`_

.. [Almi2020] S. Almi, *Irreversibility and alternate minimization in phase field
   fracture: a viscosity approach*, Z. Angew. Math. Phys. (ZAMP) 71 (2020) 128.
   `doi:10.1007/s00033-020-01357-x
   <https://doi.org/10.1007/s00033-020-01357-x>`_

.. [Pham2011] K. Pham, H. Amor, J.-J. Marigo, C. Maurini, *Gradient damage models
   and their use to approximate brittle fracture*, Int. J. Damage Mech. 20 (2011)
   618--652. `doi:10.1177/1056789510386852
   <https://doi.org/10.1177/1056789510386852>`_

.. [Tanne2018] E. Tanné, T. Li, B. Bourdin, J.-J. Marigo, C. Maurini, *Crack
   nucleation in variational phase-field models of brittle fracture*, J. Mech.
   Phys. Solids 110 (2018) 80--99. `doi:10.1016/j.jmps.2017.09.006
   <https://doi.org/10.1016/j.jmps.2017.09.006>`_

.. [Gerasimov2016] T. Gerasimov, L. De Lorenzis, *A line search assisted
   monolithic approach for phase-field computing of brittle fracture*, Comput.
   Methods Appl. Mech. Eng. 312 (2016) 276--303.
   `doi:10.1016/j.cma.2015.12.017 <https://doi.org/10.1016/j.cma.2015.12.017>`_

.. [LeonBaldelliMaurini2021] A. A. León Baldelli, C. Maurini, *Numerical
   bifurcation and stability analysis of variational gradient-damage models for
   phase-field fracture*, J. Mech. Phys. Solids 152 (2021) 104424.
   `doi:10.1016/j.jmps.2021.104424 <https://doi.org/10.1016/j.jmps.2021.104424>`_
