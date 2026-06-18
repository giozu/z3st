.. _staggered-theory:

Mathematical Foundations of the Staggered Solver
================================================

Z3ST solves coupled thermo-mechanical-damage problems with a **staggered**
scheme: within every time step it solves the thermal, mechanical, and
phase-field damage subproblems in sequence and iterates them, with adaptive
relaxation, until convergence (see :ref:`coupled-scheme`). A recurring and
entirely legitimate question is whether this is *mathematically correct*: is
solving three weak forms in turn, to convergence, the same as solving a single
coupled problem :math:`\boldsymbol{F}(\boldsymbol{u}, T, D) = \boldsymbol{0}`?
And shouldn't one instead solve *one global minimization problem*?

This page answers both questions. The short version is:

- **Yes** -- iterated to convergence, the staggered scheme solves exactly the
  same coupled first-order system that a monolithic solver would, and converges
  to the *same* critical point.
- **No** -- the correct solution concept is *not* a global minimization. For
  softening gradient-damage models, global minimization is neither attainable
  nor physically desirable; the physically correct answer comes from *local*
  minimization, which is precisely what the staggered (alternate-minimization)
  scheme realizes.
- The genuine criterion of rigor is *not* "monolithic versus staggered" but the
  **second-order stability** of the computed state -- a condition that a
  monolithic solver does not check either.

The argument below follows the variational theory of gradient-damage models
developed by Marigo, Maurini, Pham and co-workers
[MarigoMauriniPham2016]_ [PhamMarigo2010]_, the numerical analysis of its
solvers [FarrellMaurini2017]_, and the incremental variational formulation of
the *thermo*-mechanical-damage coupling [Kamagate2025]_.

.. contents:: On this page
   :local:
   :depth: 1


The variational structure: energy and three principles
-------------------------------------------------------

Z3ST's phase-field fracture model (AT1/AT2) is the gradient-damage
regularization of Griffith fracture introduced by [BourdinFrancfortMarigo2000]_.
The total energy of the body, at a fixed load, is

.. math::

   \mathcal{E}(\boldsymbol{u}, D) =
   \int_\Omega \tfrac{1}{2}\, g(D)\, \mathbb{A}_0\,
   \boldsymbol{\varepsilon}(\boldsymbol{u}) : \boldsymbol{\varepsilon}(\boldsymbol{u})
   \; \mathrm{d}x
   \;+\; \frac{G_c}{c_w} \int_\Omega
   \left( \frac{w(D)}{\ell_c} + \ell_c\, \nabla D \cdot \nabla D \right)
   \mathrm{d}x ,

where :math:`D \in [0,1]` is the damage field, :math:`g(D)` the
stiffness degradation, :math:`w(D)` the dissipation function, :math:`G_c`
the fracture toughness, and :math:`\ell_c` the internal length. The first integral
is the stored elastic energy; the second is the energy dissipated to create the
smeared crack.

The decisive point is that the quasi-static evolution is **not** governed by
"minimize :math:`\mathcal{E}`". Following the energetic theory of rate-independent
processes [MarigoMauriniPham2016]_ [PhamMarigo2010]_, the state
:math:`(\boldsymbol{u}_t, D_t)` must satisfy, at every load :math:`t`, three
principles:

#. **Irreversibility.** Damage cannot heal: :math:`\dot{D} \ge 0`. This turns
   the damage subproblem into a *variational inequality* (a unilateral, bound-
   constrained problem), not a plain equation.

#. **Stability.** The state must be a *local* minimum of the energy along every
   admissible direction: for all admissible :math:`(\boldsymbol{v}, \beta)` with
   :math:`\beta \ge D_t` there exists :math:`\bar{h} > 0` such that

   .. math::

      \mathcal{E}\big(\boldsymbol{u}_t + h(\boldsymbol{v} - \boldsymbol{u}_t),\,
      D_t + h(\beta - D_t)\big) \ge \mathcal{E}(\boldsymbol{u}_t, D_t)
      \qquad \forall\, h \in [0, \bar{h}] .

#. **Energy balance.** The energy released equals the work of external forces
   plus the energy dissipated -- no spurious creation or destruction of energy.

The first-order (necessary) optimality conditions of principles 1--2 are the
weak forms Z3ST actually solves:

.. math::

   \mathcal{E}_{\boldsymbol{u}}(\boldsymbol{u}, D; \boldsymbol{v}) = 0
   \quad \forall \boldsymbol{v},
   \qquad\qquad
   \mathcal{E}_{D}(\boldsymbol{u}, D; \beta - D) \ge 0
   \quad \forall \beta \ge D ,

i.e. mechanical equilibrium and the (irreversibility-constrained) damage
criterion. Together with the heat equation these form the coupled system
:math:`\boldsymbol{F}(\boldsymbol{u}, T, D) = \boldsymbol{0}` (with the damage
equation read as a variational inequality).

.. note::

   Already here the framing "solve one minimization problem" is too narrow: the
   correct object is an *evolution* satisfying irreversibility, stability, and
   energy balance. Minimization enters only as the *stationarity* (first-order)
   part of stability.


Separate convexity and alternate minimization
----------------------------------------------

The energy :math:`\mathcal{E}(\boldsymbol{u}, D)` is **not jointly convex**
in :math:`(\boldsymbol{u}, D)` -- the product term
:math:`g(D)\,\boldsymbol{\varepsilon}(\boldsymbol{u}):\boldsymbol{\varepsilon}(\boldsymbol{u})`
couples the two fields non-convexly. This is the root of every numerical
difficulty in phase-field fracture, and it is also why a *single* global
minimization is ill-posed (many local minimizers).

However -- and this is the key structural fact -- :math:`\mathcal{E}` **is convex
separately** in each field when the other is frozen [FarrellMaurini2017]_:

- with :math:`D` fixed, :math:`\boldsymbol{u} \mapsto \mathcal{E}` is a
  standard (convex) linear-elasticity problem;
- with :math:`\boldsymbol{u}` fixed, :math:`D \mapsto \mathcal{E}` is a
  convex, bound-constrained Helmholtz-type problem.

The **alternate-minimization** algorithm of [BourdinFrancfortMarigo2000]_ -- which
is exactly the staggered scheme -- exploits this: alternately minimize over one
field with the other fixed, and iterate. Because each subproblem has a *unique*
solution that *lowers* the energy, the algorithm produces a sequence of states
with monotonically decreasing energy, and therefore **converges to a critical
point** of the coupled problem [FarrellMaurini2017]_. In the words of Farrell and
Maurini, *"at each iteration before convergence, the optimization subproblem has a
unique solution with a lower energy, and thus, the algorithm converges
monotonically to a stationary point."*

This is the precise sense in which the staggered scheme is *well-posed*: each
block is a convex problem with a unique solution, and the outer iteration is a
non-linear block Gauss--Seidel iteration whose fixed point is a solution of the
full coupled system.

Beyond convergence to *a* critical point, the limit of the alternate-minimization
evolution has been characterized rigorously. As the time step is refined, the
time-discrete staggered solutions converge to a parametrized evolution that
satisfies Griffith's criterion and is thermodynamically consistent with damage
irreversibility [Knees2017]_ [Almi2019]_. Crucially, enforcing irreversibility by
the *a posteriori* pointwise truncation that Z3ST uses,
:math:`D^{n+1} = \min(1, \max(D^{n+1}, D^{n}))`, is exactly the variant analysed
in [Almi2020]_ and shown to converge to the correct evolution. The staggered
scheme is therefore not a mere numerical expedient: its converged result is a
provably sound approximation of the underlying brittle-fracture evolution.


Equivalence with a monolithic solve, at convergence
----------------------------------------------------

Let the staggered loop run until the change in *every* field, between two
successive sweeps, drops below tolerance simultaneously:

.. math::

   \| T^{k+1} - T^{k} \| < \mathrm{tol}_T, \qquad
   \| \boldsymbol{u}^{k+1} - \boldsymbol{u}^{k} \| < \mathrm{tol}_u, \qquad
   \| D^{k+1} - D^{k} \| < \mathrm{tol}_D .

At that point an additional sweep no longer changes the state: the iteration has
reached a **fixed point of the staggered map**, which is by definition a solution
of the coupled first-order system
:math:`\boldsymbol{F}(\boldsymbol{u}, T, D) = \boldsymbol{0}`. A monolithic
Newton solver applied to the same system would converge to the *same* critical
point (within the same basin of attraction). The two approaches are therefore
**equivalent at convergence** -- they are two algorithms for the same discrete
solution, not two different solutions. This has been confirmed numerically: a
monolithic quasi-Newton (BFGS) solve of the coupled system *"yields identical
results to the AM/staggered solver"* [Wu2020]_.

Two honest qualifications:

- **Iterate to convergence, not once.** The equivalence holds only if the inner
  staggered loop is iterated to convergence at each time step. A *single*
  alternate-minimization pass per step (as in some explicit-in-time schemes)
  introduces a splitting error and is *not* equivalent to solving
  :math:`\boldsymbol{F} = \boldsymbol{0}` [FarrellMaurini2017]_. Z3ST iterates the
  loop to a tolerance (`stag_tol`) and reports whether the step converged; this
  is what makes the equivalence hold.

- **Increment versus residual.** Z3ST's stopping test is on the *increment* of
  the fields. At the fixed point increment and coupled residual both vanish, so
  the criteria agree; for a fully rigorous statement one would additionally
  monitor the residual of the coupled system. Reporting that residual is the
  single cheapest way to make the equivalence argument airtight.

The robustness advantage of the staggered scheme is also a direct consequence of
the non-convexity: the *plain* monolithic Newton Jacobian is *indefinite*, so
Newton's method "does not converge unless extremely small continuation steps are
taken" [FarrellMaurini2017]_. Alternate minimization is unconditionally robust
precisely because it only ever solves convex subproblems. Two complementary
families of remedies exist in the literature, and both *confirm* rather than
contradict the staggered result:

- **Accelerate the staggered scheme.** Over-relaxation (a non-linear
  successive-over-relaxation reading of Gauss--Seidel) [FarrellMaurini2017]_,
  Anderson acceleration [Storvik2020]_, or composing alternate minimization with
  Newton once inside its basin of attraction [FarrellMaurini2017]_ cut the
  iteration count by large factors while reaching the *same* converged solution.
  Z3ST's adaptive (Aitken-type) relaxation is exactly such an over-relaxed
  alternate minimization.
- **Make the monolithic solve robust.** Replacing plain Newton with a quasi-Newton
  (BFGS) scheme [Wu2020]_ [Kristensen2020]_, an exact line search [Heinzmann2025]_,
  or a fracture-energy arc-length method with under-relaxation [Bharali2022]_
  restores convergence on the indefinite system -- and again reaches the same
  critical point as the staggered scheme.

The practical takeaway is that staggered-versus-monolithic is a question of
*efficiency and robustness*, not of *which solution is correct*: at convergence
the two coincide.


Why not global minimization
----------------------------

It is tempting to argue that the "mathematically pure" formulation is a single
*global* minimization of :math:`\mathcal{E}`. For softening gradient-damage models
this is incorrect on two levels.

**Global minimization is not attainable and is physically wrong for nucleation.**
Tanné, Li, Bourdin, Marigo and Maurini [Tanne2018]_ show that, for a bar of size
:math:`L`, *global* minimization nucleates a crack only at the Griffith-like load,
with critical stress

.. math::

   \sigma_c^{\text{global}} \sim \sqrt{\frac{2 G_c E}{L}} ,

which diverges as :math:`L \to 0` and produces a spurious structural size effect --
the very pathology of Griffith/LEFM that the model was meant to cure. **Local**
minimization at a *fixed* internal length :math:`\ell_c` instead nucleates at the
material strength,

.. math::

   \sigma_c^{\text{local}} \sim \sqrt{\frac{G_c E}{\ell_c}} ,

an intrinsic critical stress independent of structural size. Alternate
minimization realizes this local minimization. In other words, the staggered
scheme does not merely *approximate* the "right" global answer -- for crack
*nucleation* it gives the *physically correct* answer that global minimization
would miss.

This is why, in Z3ST and in the gradient-damage literature, :math:`\ell_c` is a
**material parameter** linked to the tensile strength
(:math:`\sigma_c \sim \sqrt{G_c E / \ell_c}`) [Tanne2018]_ [PhamMarigo2010]_, not a
purely numerical regularization to be driven to zero.

**With temperature, the coupled problem is a saddle point, not a minimum.** When
the heat equation is added, the incremental variational principle for the
thermo-mechanical-damage problem is convex in :math:`(\boldsymbol{u}, D)` but
*concave* in the temperature :math:`T` [Kamagate2025]_. The correct statement is
therefore a min-max,

.. math::

   (\boldsymbol{u}, T, D)
   = \arg\,\inf_{\boldsymbol{u}, D}\, \sup_{T}\;
   I_n(\boldsymbol{u}, T, D) ,

so even *ideally* there is no single global minimization to solve. Notably,
Kamagaté et al. [Kamagate2025]_ derive this rigorous variational framework for
exactly the problem Z3ST targets (gradient damage + thermoelasticity + heat
conduction, applied to thermal-shock plate cracking) and implement it -- in
FEniCS -- with a **semi-staggered** algorithm iterated to convergence, identical
in spirit to Z3ST's solver.


The real criterion of rigor: second-order stability
----------------------------------------------------

Solving :math:`\boldsymbol{F}(\boldsymbol{u}, T, D) = \boldsymbol{0}` -- whether by
the staggered scheme or by a monolithic solver -- delivers a *first-order*
critical point. For softening models this is **necessary but not sufficient**:
because the energy is non-convex, a critical point may be *unstable*, and standard
first-order algorithms can converge to such an unstable branch
[LeonBaldelliMaurini2021]_.

The rigorous criterion is the **second-order stability** condition (principle 2
above), i.e. positivity of the reduced Hessian of the energy restricted to the
cone of admissible (irreversibility-respecting) directions -- the criterion first
formulated for gradient-damage models by [Pham2011]_. León Baldelli and
Maurini [LeonBaldelliMaurini2021]_ formulate this as an eigenvalue problem on the
reduced Hessian and use it to *"filter out unstable solutions provided by standard
first-order minimization algorithms."*

The important consequence for Z3ST users is that **this gap is identical for the
staggered and the monolithic scheme** -- neither tests stability by default, so
switching solver does not close it.

What is more surprising, a *more robust* solver can be *worse* at exposing the
problem. The robust monolithic schemes of the previous section owe their
robustness to a **positive-definite approximation of the Hessian**: quasi-Newton
(BFGS) updates, by construction, build a positive-definite operator, and that is
precisely what lets them converge on the indefinite phase-field system where plain
Newton stalls. But the indefiniteness they smooth over is exactly the signature of
loss of stability. Lacking the *exact* second variation, a quasi-Newton iteration
can therefore converge smoothly onto an **unstable** equilibrium branch and report
success, with no indication that the computed state is not a physical minimum
[Terzi2025]_. Robustness of the *solve* and detection of *instability* are thus in
tension: the trick that buys convergence also hides the bifurcation.

.. warning::

   Convergence of the non-linear solver -- staggered or monolithic, Newton or
   quasi-Newton -- certifies only that a *first-order* critical point has been
   found. It says nothing about whether that state is *stable* (a local minimum).
   A quasi-Newton monolithic solver may even mask the instability, because its
   positive-definite Hessian surrogate cannot see the negative eigenvalue that
   marks the loss of stability. The only reliable diagnosis is an explicit
   **second-order check**: the sign of the smallest eigenvalue of the reduced
   Hessian, restricted to the irreversibility cone [LeonBaldelliMaurini2021]_
   [Terzi2025]_.

The rigorous upgrade path is therefore *not* "rewrite the solver as monolithic"
but "add a second-order stability check" -- a constrained-eigenvalue problem
solved once on the converged first-order state, cheap relative to the time
stepping itself. This is the actual research frontier and a natural future
addition to Z3ST's verification suite; it would let Z3ST not only *find* an
evolution but *certify* that the branch it follows is the physically stable one.


Summary: why the Z3ST staggered solver is sound
-----------------------------------------------

.. note::

   - The phase-field energy is **separately convex**; the staggered scheme is the
     associated **alternate minimization**, in which every subproblem is convex
     and uniquely solvable, and the energy decreases monotonically to a critical
     point [BourdinFrancfortMarigo2000]_ [FarrellMaurini2017]_.
   - Iterated to convergence, the staggered scheme solves the **same coupled
     system** :math:`\boldsymbol{F}(\boldsymbol{u}, T, D) = \boldsymbol{0}` as a
     monolithic solver, and reaches the **same critical point**.
   - It is **more robust** than a monolithic Newton solve, whose Jacobian is
     indefinite; Z3ST's adaptive relaxation is an over-relaxed alternate
     minimization in the sense of [FarrellMaurini2017]_.
   - It realizes **local** minimization, which for softening models is the
     *physically correct* concept (correct nucleation at the material strength);
     **global** minimization is neither attainable nor desirable [Tanne2018]_.
   - With temperature the coupled problem is a **saddle point**, and the rigorous
     variational reference for precisely this thermo-mechanical-damage problem is
     solved -- in FEniCS -- with a **staggered** algorithm too [Kamagate2025]_.
   - The genuine criterion beyond :math:`\boldsymbol{F} = \boldsymbol{0}` is
     **second-order stability** [LeonBaldelliMaurini2021]_, a check that is
     independent of the staggered-versus-monolithic choice.

For the concrete implementation -- block ordering, relaxation, solvers, and
convergence test -- see :ref:`coupled-scheme`.


References
----------

.. [BourdinFrancfortMarigo2000] B. Bourdin, G. A. Francfort, J.-J. Marigo,
   *Numerical experiments in revisited brittle fracture*, J. Mech. Phys. Solids
   48 (2000) 797--826. `doi:10.1016/S0022-5096(99)00028-9
   <https://doi.org/10.1016/S0022-5096(99)00028-9>`_

.. [PhamMarigo2010] K. Pham, J.-J. Marigo, *Approche variationnelle de
   l'endommagement: I. Les concepts fondamentaux / II. Les modèles à gradient*,
   C. R. Mécanique 338 (2010). `doi:10.1016/j.crme.2010.02.001
   <https://doi.org/10.1016/j.crme.2010.02.001>`_

.. [MarigoMauriniPham2016] J.-J. Marigo, C. Maurini, K. Pham, *An overview of
   the modelling of fracture by gradient damage models*, Meccanica 51 (2016)
   3107--3128. `doi:10.1007/s11012-016-0538-4
   <https://doi.org/10.1007/s11012-016-0538-4>`_

.. [FarrellMaurini2017] P. E. Farrell, C. Maurini, *Linear and nonlinear solvers
   for variational phase-field models of brittle fracture*, Int. J. Numer.
   Methods Eng. 109 (2017) 648--667. `doi:10.1002/nme.5300
   <https://doi.org/10.1002/nme.5300>`_

.. [Tanne2018] E. Tanné, T. Li, B. Bourdin, J.-J. Marigo, C. Maurini, *Crack
   nucleation in variational phase-field models of brittle fracture*, J. Mech.
   Phys. Solids 110 (2018) 80--99. `doi:10.1016/j.jmps.2017.09.006
   <https://doi.org/10.1016/j.jmps.2017.09.006>`_

.. [LeonBaldelliMaurini2021] A. A. León Baldelli, C. Maurini, *Numerical
   bifurcation and stability analysis of variational gradient-damage models for
   phase-field fracture*, J. Mech. Phys. Solids 152 (2021) 104424.
   `doi:10.1016/j.jmps.2021.104424 <https://doi.org/10.1016/j.jmps.2021.104424>`_

.. [Kamagate2025] B. Kamagaté, L. Cheng, R. Abdelmoula, E. Danho, D. Kondo,
   *An incremental variational method to the coupling between gradient damage,
   thermoelasticity and heat conduction*, C. R. Mécanique 353 (2025) 1063--1084.
   `doi:10.5802/crmeca.325 <https://doi.org/10.5802/crmeca.325>`_

.. [Pham2011] K. Pham, H. Amor, J.-J. Marigo, C. Maurini, *Gradient damage models
   and their use to approximate brittle fracture*, Int. J. Damage Mech. 20 (2011)
   618--652. `doi:10.1177/1056789510386852
   <https://doi.org/10.1177/1056789510386852>`_

.. [Knees2017] D. Knees, M. Negri, *Convergence of alternate minimization schemes
   for phase-field fracture and damage*, Math. Models Methods Appl. Sci. (M3AS)
   27 (2017) 1743--1794. `doi:10.1142/S0218202517500312
   <https://doi.org/10.1142/S0218202517500312>`_

.. [Almi2019] S. Almi, S. Belz, M. Negri, *Convergence of discrete and continuous
   unilateral flows for Ambrosio--Tortorelli energies and application to
   mechanics*, ESAIM: Math. Model. Numer. Anal. (M2AN) 53 (2019) 659--699.
   `doi:10.1051/m2an/2018057 <https://doi.org/10.1051/m2an/2018057>`_

.. [Almi2020] S. Almi, *Irreversibility and alternate minimization in phase field
   fracture: a viscosity approach*, Z. Angew. Math. Phys. (ZAMP) 71 (2020) 128.
   `doi:10.1007/s00033-020-01357-x
   <https://doi.org/10.1007/s00033-020-01357-x>`_

.. [Wu2020] J.-Y. Wu, Y. Huang, V. P. Nguyen, *On the BFGS monolithic algorithm
   for the unified phase field damage theory*, Comput. Methods Appl. Mech. Eng.
   360 (2020) 112704. `doi:10.1016/j.cma.2019.112704
   <https://doi.org/10.1016/j.cma.2019.112704>`_

.. [Kristensen2020] P. K. Kristensen, E. Martínez-Pañeda, *Phase field fracture
   modelling using quasi-Newton methods and a new adaptive step scheme*, Theor.
   Appl. Fract. Mech. 107 (2020) 102446. `doi:10.1016/j.tafmec.2019.102446
   <https://doi.org/10.1016/j.tafmec.2019.102446>`_

.. [Bharali2022] R. Bharali, S. Goswami, C. Anitescu, T. Rabczuk, *A robust
   monolithic solver for phase-field fracture integrated with fracture energy
   based arc-length method and under-relaxation*, Comput. Methods Appl. Mech. Eng.
   394 (2022) 114927. `doi:10.1016/j.cma.2022.114927
   <https://doi.org/10.1016/j.cma.2022.114927>`_

.. [Storvik2020] E. Storvik, J. W. Both, J. M. Sargado, J. M. Nordbotten, F. A.
   Radu, *An accelerated staggered scheme for variational phase-field models of
   brittle fracture*, Comput. Methods Appl. Mech. Eng. 381 (2021) 113822.
   `doi:10.1016/j.cma.2021.113822 <https://doi.org/10.1016/j.cma.2021.113822>`_

.. [Heinzmann2025] J. Heinzmann, F. Vicentini, P. Carrara, L. De Lorenzis,
   *Iterative convergence in phase-field brittle fracture computations: exact line
   search is all you need*, arXiv:2511.23064 (2025).
   `arXiv:2511.23064 <https://arxiv.org/abs/2511.23064>`_

.. [Terzi2025] M. M. Terzi, O. U. Salman, D. Faurie, A. A. León Baldelli,
   *Navigating local minima and bifurcations in brittle thin film systems with
   irreversible damage*, arXiv:2409.04307 (2025).
   `arXiv:2409.04307 <https://arxiv.org/abs/2409.04307>`_
