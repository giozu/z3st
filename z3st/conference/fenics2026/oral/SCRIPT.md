# Speaker script — Z3ST oral (FEniCS 2026, 10-minute slot)

Format: 10-minute oral, Paris, 17–19 June 2026. Budget 10 content slides; the
deck is tighter since the PCMI slide (6b) was added, so move briskly. Backup
slides are not presented unless a question pulls them up. Aim to *finish by 9:45*
so you are never cut off.

Timings are cumulative targets (mm:ss). The most compressible cuts, in order:
the cluster-dynamics line on slide 6, the teaching-cases aside on slide 7, and
the gap-conductance detail on slide 6 (PCMI / slide 6b already carries the
across-bodies story). Drop these first if you are behind at slide 6b.

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

## 6 · Coupled multi-physics, across bodies — 5:25 → 6:10
"At heart this is a multi-physics framework: thermal, mechanical and phase-field
damage, genuinely two-way coupled — temperature drives the thermal strain and
stress, and damage degrades both the elastic and the thermal-stress terms. And it
works across separate bodies — gap conductance carries heat across a gap between two
regions, a fuel pellet and its cladding. *(Compressible: cluster dynamics — a
newer, exploratory capability evolving a defect population in cluster-size space,
our path to a micro-to-continuum link.)*"

## 6b · PCMI: pellet–clad contact, verified — 6:10 → 7:05  ← new
"And here is the *mechanical* side of 'across bodies' — pellet–cladding
interaction. As the pellet heats it expands, closes the gap, and contacts the
cladding. Contact is a penalty: a pressure proportional to the penetration,
applied as equal and opposite tractions on the two faces. Nothing is prescribed —
the pressure *emerges*, and you can watch it on the left: the gap closes, contact
switches on, and the cladding is pushed outward — that is load transfer. And it
feeds back thermally: contact raises the gap conductance, so the fuel cools the
moment it touches. The point for this room is on the right — it's *verified*:
against the analytical Lamé interference-fit pressure, to three and a half percent,
with the stress state confirmed plane-stress. And the penalty tangent? The same AD
path — `ufl.derivative` — no hand-coded contact Jacobian."

## 7 · Verified and reproducible — 7:05 → 8:00
"None of this is useful unless it's trustworthy, and verification is something I
care about a lot. Every solver and every physics — thermal, mechanical, plasticity,
fracture — is checked against an analytical solution. That's about fifty cases, each
a numerical-versus-analytical comparison, and they re-run on every single commit, so
a regression can't sneak in. *(Compressible: the teaching cases make the point — a
1D bar on a line mesh and the same bar in full 3D give identical displacement,
regime selection not re-meshing.)* And it's built for experimentation: hot-reload
tolerances and relaxation mid-run, stream force-displacement per step, inspect the
matrix directly."

## 8 · Hero case — 8:00 → 9:05
"Here's everything together: thermal-shock cracking of a UO2 pellet. Top row: a
cold-contact wedge cools the rim, and the tensile hoop-stress ring it sets up drives
the cracking. Bottom row is the payoff — on the left, the simulated damage: a set of
discrete radial cracks at the rim; on the right, a cross-section of a real UO2 pellet
showing exactly that pattern. So this isn't just a pretty field — the model
reproduces what you see in a real sample, and it matches McClenny's Figure 8. 2D
plane-strain half-disc, AT1, Amor split, hybrid, fully coupled temperature to
elastic strain to damage — all from one `input.yaml`. (Backup has the quantitative
match to McClenny 7b.)"

## 9 · Open + demo invite — 9:05 → 9:45
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
- *Which contact method / is the pressure trustworthy?* Explicit penalty —
  pressure proportional to penetration, equal-and-opposite tractions, the tangent
  from AD (no hand-coded contact Jacobian). Verified against the analytical Lamé
  interference-fit to 3.5%, with the stress state confirmed plane-stress.
  Uniform-pressure today; per-facet (to resolve axial ridging) is the next step.
  No mesh-cutting / XFEM needed: the two bodies are separate meshes, so the
  contact interface is a mesh boundary, not an intra-element discontinuity.
- *Show me plasticity / what about crystal plasticity?* Backup "plasticity verified" —
  J2 vs analytical 2D plane strain, and the single-crystal case saturating to the
  analytical sigma_sat. Same AD path; the slip-system Jacobian is automatic.
- *Is the phase-field validated against a standard benchmark?* Backup "single-edge-
  notched shear" — the curved crack reproduces Miehe et al. (CMAME 199, 2010) at
  u_x = 15 microns.
- *Is the coupling actually verified, not just claimed?* Backup "coupled thermo-
  mechanics" — temperature and thermal stress vs analytical across the slab.
- *Can it do fuel performance — burnup, swelling, gap closure over an irradiation?*
  Backup "burnup-driven PCMI" — burnup accumulates as a per-material state field,
  swelling is an eigenstrain callable that reads it ("fuel is a material" — no
  solver change), and a 3-year rod history at 25 kW/m closes the 65 µm gap at
  ~30 MWd/kgU with ~39 MPa emergent contact pressure; the mean burnup is verified
  against the closed form to 1e-7. The peak fuel temperature *drops* on contact —
  the thermal feedback, live in the figure.

## Backup slides (do not present; flip to on a question)
1. staggered solver loop · 2. constitutive-model table · 3. plasticity verified
(J2 + crystal) · 4. coupled thermo-mechanics verified · 5. burnup-driven PCMI
(3-year rod history) · 6. SEN shear, Miehe (2010) ·
7. hero case verified vs McClenny (2022), Fig. 7b.

## Delivery notes
- The AD slide (4), PCMI (6b) and the hero case (8) are what the audience
  remembers. Slow down on those; speed through 3 and 7 if you must.
- On 6b: land "*emergent* pressure" and "verified to 3.5% against Lamé" — the two
  phrases that make it stick. Point once at the gap-closure curve, once at the
  verification curve; don't narrate both panels.
- Do not read the code aloud line by line — point at `ufl.diff(psi, F)` and say it
  once; same for the contact-tangent callback to AD on 6b.
- Practised cold, with 6b added this runs ~9:45. The compressible asides (cluster
  line on 6, teaching aside on 7) buy ~20 s of slack if the committee clock is
  tight on the 10-minute format.
