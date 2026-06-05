# Speaker script — Z3ST oral (FEniCS 2026, 10-minute slot)

Format: 10-minute oral, Paris, 17–19 June 2026. Budget ~9 content slides at
~60–65 s each, leaving a few seconds of slack. Backup slides are not presented
unless a question pulls them up. Aim to *finish by 9:30* so you are never cut off.

Timings are cumulative targets (mm:ss). If you are behind at slide 6, drop the
gap-conductance sentence and the mesh-sensitivity aside — they are the most
compressible.

---

## 1 · Title — 0:00 → 0:20
"Good morning. I'm Giovanni Zullo from Politecnico di Milano, and I'll show you
Z3ST — a small FEniCSx framework for coupled thermo-mechanical analysis where the
constitutive physics is written as energies and differentiated automatically."

## 2 · The problem — 0:20 → 1:40
"The components I care about — nuclear fuel, thermal-shock cracking — fail where
heat, deformation and fracture interact. The production codes, MOOSE and BISON,
are validated but heavy: changing a weak form is a real undertaking. FEniCSx is
the opposite — expressive forms, symbolic differentiation built in — but there is
no turnkey coupled thermo-mechanical-fracture driver. Z3ST fills that gap: a
configuration-driven solver on top of FEniCSx, where new physics means a new
energy, not a hand-coded tangent."

## 3 · Design at a glance — 1:40 → 2:50
"The architecture is one driver, `Spine`, that composes the physics as mixins —
thermal, mechanical, damage, plasticity, cluster dynamics, gap. Everything else is
data: geometry, boundary conditions, materials and solver settings live in a
single `input.yaml`. The code is material-agnostic — it doesn't care which solid
you give it — and a material property can be a constant or a Python function of
state, like a temperature-dependent conductivity k of T. That's a gentle on-ramp:
a student can add a new material with a few lines of Python. The same model runs in
1D, 2D, 3D, axisymmetric or plane-stress — you pick a regime, you don't rewrite the
model. You launch a case with `python -m z3st`."

## 4 · The core idea (AD) — 2:50 → 4:10  ← the slide that matters
"This is the heart of it. A constitutive model is entered as a strain-energy
density. Here, neo-Hookean. The first Piola-Kirchhoff stress is *one line* —
`P = ufl.diff(psi, F)` — the symbolic derivative of the energy with respect to the
deformation gradient. The residual is the usual inner product, and the Newton
tangent comes from `ufl.derivative` — again automatic, no analytic Jacobian. So
elasticity, hyperelasticity, J2 and crystal plasticity, phase-field — all share
this one path: write the energy, get the stress and the tangent for free. We even
have a single-crystal slip-system viscoplasticity demo whose exact Jacobian comes
straight from AD — the one case nobody wants to differentiate by hand."

## 5 · Coupled thermo-mechanics — 4:10 → 5:25
"The coupling is staggered. Each load step alternates: a thermal solve with
temperature-dependent conductivity and a heat source; a mechanical solve that
subtracts the thermal eigenstrain so stress is driven by the elastic strain; and a
damage solve driven by the positive elastic energy, kept irreversible. The trick
that makes it robust is adaptive relaxation — per-field factors that grow when the
residual shrinks and back off when it oscillates. And crucially, damage degrades
both the elastic and the thermal-stress contributions through a consistent g(D) —
without that, thermal-shock fracture is intractable."

## 6 · Coupled multi-physics, across bodies — 5:25 → 6:40
"At heart this is a multi-physics framework: thermal, mechanical and phase-field
damage, genuinely two-way coupled — temperature drives the thermal strain and
stress, and damage degrades both the elastic and the thermal-stress terms. The
second point is that it works across separate bodies: gap conductance carries heat
across a gap between two regions — think a fuel pellet and its cladding — as a
thermal interface condition. That multi-body heat transfer is essential for real
fuel-performance modelling. And one newer, exploratory capability: cluster
dynamics, which evolves a defect population in cluster-size space — it is
implemented but not yet investigated, and it is our path towards a genuine
micro-to-continuum link."

## 7 · Verified and reproducible — 6:40 → 7:50
"None of this is useful unless it's trustworthy, and verification is something I
care about a lot. Every solver and every physics — thermal, mechanical, plasticity,
fracture — is checked against an analytical solution. That's about fifty cases, each
a numerical-versus-analytical comparison, and they re-run on every single commit, so
a regression can't sneak in. The teaching cases make the point nicely: a 1D bar on a
line mesh and the same bar in full 3D give the identical displacement — regime
selection, not re-meshing. And it's built for experimentation: you can hot-reload
tolerances and relaxation factors mid-run, stream force-displacement curves per step,
and inspect the matrix directly."

## 8 · Hero case — 7:50 → 9:00
"Here's everything together: thermal-shock cracking of a UO2 pellet. Top row: a
cold-contact wedge cools the rim, and the tensile hoop-stress ring it sets up drives
the cracking. Bottom row is the payoff — on the left, the simulated damage: a set of
discrete radial cracks at the rim; on the right, a cross-section of a real UO2 pellet
showing exactly that pattern. So this isn't just a pretty field — the model
reproduces what you see in a real sample, and it matches McClenny's Figure 8. 2D
plane-strain half-disc, AT1, Amor split, hybrid, fully coupled temperature to
elastic strain to damage — all from one `input.yaml`. (Backup has the quantitative
match to McClenny 7b.)"

## 9 · Open + demo invite — 9:00 → 9:45
"Z3ST is open under Apache-2.0 on GitHub, documented, archived on Zenodo with a DOI.
If you'd like to see it run, come to the software demonstration at the poster
session — I'll take a coupled case end to end on my laptop and tune it while it
solves. Thank you — I'm happy to take questions."

---

## Anticipated questions (have the answer, not the slide)
- *Monolithic vs staggered?* Staggered with adaptive relaxation; monolithic
  phase-field Newton is on the roadmap (backup slide 1 has the loop).
- *How does it compare to BISON / MOOSE?* Not a replacement — a lightweight,
  hackable research framework; energy-first weak forms, fast to extend.
- *Performance / parallelism?* PETSc + MPI under dolfinx; HYPRE BoomerAMG on the
  larger cases. Not the focus of the talk.
- *Which elements?* Lagrange P1 for displacement/temperature; DG1 for cluster
  dynamics; quadrature elements for plasticity history.
- *Validation vs verification?* The suite is verification (vs analytical /
  reference); validation against experiment is application-specific.
- *Show me plasticity / what about crystal plasticity?* Backup "plasticity verified" —
  J2 vs analytical 2D plane strain, and the single-crystal case saturating to the
  analytical sigma_sat. Same AD path; the slip-system Jacobian is automatic.
- *Is the phase-field validated against a standard benchmark?* Backup "single-edge-
  notched shear" — the curved crack reproduces Miehe et al. (CMAME 199, 2010) at
  u_x = 15 microns.
- *Is the coupling actually verified, not just claimed?* Backup "coupled thermo-
  mechanics" — temperature and thermal stress vs analytical across the slab.

## Backup slides (do not present; flip to on a question)
1. staggered solver loop · 2. constitutive-model table · 3. plasticity verified
(J2 + crystal) · 4. coupled thermo-mechanics verified · 5. SEN shear, Miehe (2010) ·
6. hero case verified vs McClenny (2022), Fig. 7b.

## Delivery notes
- The AD slide (4) and the hero case (8) are the two the audience remembers. Slow
  down on those; speed through 3 and 7 if you must.
- Do not read the code aloud line by line — point at `ufl.diff(psi, F)` and say it
  once.
- Practised cold, this runs ~9:30. That is the target; the committee is strict on
  the 10-minute format.
