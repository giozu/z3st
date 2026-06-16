# Speaker script — Z3ST oral (FEniCS 2026, 10-minute slot)

Format: 10-minute oral, Paris, 17–19 June 2026. Budget 10 content slides; the
deck is tighter since the PCMI slide (6) was added, so move briskly. Backup
slides are not presented unless a question pulls them up. Aim to *finish by 9:45*
so you are never cut off.

Timings are cumulative targets (mm:ss). The most compressible cuts, in order:
the cluster-dynamics line on slide 5, the teaching-cases aside on slide 7, and
the gap-conductance detail on slide 5 (PCMI / slide 6 already carries the
across-bodies story). Drop these first if you are behind at slide 6.

---

## Title slide — 0:00 → 0:20  (no footer number)
"Good morning. I'm Giovanni Zullo, a post-doc at Politecnico di Milano. I'm a nuclear
engineer — I work on modelling and simulation of nuclear materials: fuel pellets,
cladding, advanced materials for extreme environments, and the behaviour of fission
products. Today I'll show you Z3ST — an open FEniCSx framework for coupled
thermo-mechanical analysis and fracture, where the constitutive physics is written as
energies and differentiated automatically."

## 1 · The problem — 0:20 → 1:40
"In my field — fuel performance — the established codes are mostly proprietary or
restricted: TRANSURANUS, FALCON, BISON. If you want something open and
finite-element, you're essentially left with MOOSE — powerful, but a heavy,
monolithic beast where adding or changing a weak form is a real undertaking.
FEniCSx is the opposite: expressive forms, automatic differentiation built in, and a
great open community — but there's no turnkey coupled thermo-mechanical-fracture
driver for fuel. Z3ST is that missing layer: an open, FEniCSx-native framework where
the whole coupled apparatus is already assembled — and verified on about forty-five
analytic cases — so a new problem is a YAML file, and new physics is a new energy,
not a hand-coded tangent."

## 2 · Design at a glance — 1:40 → 2:50
"The architecture is one driver, `Spine`, that composes the physics as mixins —
thermal, mechanical, damage, plasticity, cluster dynamics, gap. Everything else is
data: geometry, boundary conditions, materials and solver settings are all YAML,
gathered by a top-level `input.yaml`. The code is material-agnostic — it doesn't care which solid
you give it — and a material property can be a constant or a Python function of
state, like a temperature-dependent conductivity k of T. That's a gentle on-ramp:
a student can add a new material with a few lines of Python. The same model runs in
1D, 2D, 3D, axisymmetric or plane-stress — you pick a regime, you don't rewrite the
model. You launch a case with `python -m z3st`."

## 3 · The core idea (AD) — 2:50 → 4:10  ← the slide that matters
"This is the engine — and yes, you all know UFL differentiates forms. The point is
what that *buys* you. A constitutive model is a strain-energy density — here,
neo-Hookean. The stress is *one line*, `P = ufl.diff(psi, F)`; the Newton tangent
is another `ufl.derivative`. No hand-coded Jacobian, ever. And because *every*
constitutive law — elasticity, J2 and crystal plasticity, phase-field — lives
on this one path, I could assemble the entire coupled apparatus once and reuse it.
*That's* the real point: adding physics is writing an energy, and the coupling is
already done for you — you don't re-derive it, you don't re-wire it. Even the
slip-system crystal-plasticity Jacobian, the one nobody wants by hand, comes
straight from AD."

## 4 · Coupled thermo-mechanics — 4:10 → 5:25
"The coupling is staggered. Each load step alternates: a thermal solve with
temperature-dependent conductivity and a heat source; a mechanical solve that
subtracts the thermal eigenstrain so stress is driven by the elastic strain; and a
damage solve driven by the positive elastic energy, kept irreversible. The trick
that makes it robust is adaptive relaxation — per-field factors that grow when the
residual shrinks and back off when it oscillates. And if a step still stalls — too
many physics fighting at once — the solver can roll back and halve the time step on
its own, pushing through the stiff interval instead of giving up. Crucially, the
damage degradation g(D) acts on the full stress — the thermal-stress part included,
since the thermal strain sits inside the elastic strain — without that, thermal-shock
fracture is intractable."

## 5 · Coupled multi-physics, across bodies — 5:25 → 6:10
"At heart this is a multi-physics framework: thermal, mechanical and phase-field
damage, genuinely two-way coupled — temperature drives the thermal strain and
stress, and damage degrades the stress, the thermal-stress part included. And it
works across separate bodies — gap conductance carries heat across a gap between two
regions, a fuel pellet and its cladding. *(Compressible: cluster dynamics — a
newer, exploratory capability evolving a defect population in cluster-size space,
our path to a micro-to-continuum link.)*"

## 6 · PCMI: a fuel rod over a full irradiation — 6:10 → 7:05
"And here is the *mechanical* side of 'across bodies', taken all the way — a fuel
rod followed over an entire irradiation, eighteen hundred days at twenty kilowatts
per metre. Watch the left panel: as burnup builds, the pellet swells and the gap
closes; the moment it touches, a contact pressure *emerges* — nothing is
prescribed, it is a penalty proportional to the penetration. But look what it does
— it does not climb without bound, it *saturates*, because the cladding creeps
under irradiation and relaxes the contact stress about as fast as swelling builds
it. That self-limiting plateau is a real fuel-performance signature, and it falls
straight out of the coupled solve. The right panel is the thermal feedback: the
instant contact engages, the gap conductance jumps and the peak fuel temperature
*drops* — from eleven-twenty-six down to nine-seventy-six kelvin. And every
ingredient here — Fink conductivity, swelling, creep, gap, contact — is configured
in YAML; the contact is an explicit penalty for now — pressure proportional to
penetration, driven to a fixed point by the staggered loop. It is *verified* against
the analytical Lamé interference-fit pressure, to three and a half percent — I have
that in the backup if you want it."

## 7 · Verified and reproducible — 7:05 → 8:00
"None of this is useful unless it's trustworthy, and verification is something I
care about a lot. Every solver and every physics — thermal, mechanical, plasticity,
fracture — is checked against an analytical solution. That's about forty-five cases,
each a numerical-versus-analytical comparison; a fast subset re-runs in CI on every
push, the full suite locally, so a regression can't sneak in. *(Compressible: the
teaching cases make the point — a 1D bar and the same bar in full 3D give identical
displacement; you pick a regime, you don't rewrite the model.)* And it's built for
experimentation: hot-reload
tolerances and relaxation mid-run, stream force-displacement per step, inspect the
matrix directly."

## 8 · Hero case — 8:00 → 9:05
"Here's everything together: thermal-shock cracking of a UO2 pellet. Top row: a
cold-contact wedge cools the rim, and the tensile hoop-stress ring it sets up drives
the cracking. Bottom row is the payoff — on the left, the simulated damage: a set of
discrete radial cracks at the rim; on the right, a cross-section of a real UO2 pellet
showing exactly that pattern. So this isn't just a pretty field — the model
reproduces what you see in a real sample, the thermal-shock cracking McClenny et al.
observed experimentally. 2D plane-strain half-disc, AT1, Amor split, hybrid, fully
coupled temperature to elastic strain to damage — all from one `input.yaml`. (Backup
has the quantitative match to McClenny et al. 2022.)"

## 9 · Open + future work + demo invite — 9:05 → 9:45
"So, to close: Z3ST is the open, FEniCSx-native layer where the coupled apparatus is
already assembled and verified — your problem is a YAML file. It's Apache-2.0 on
GitHub, documented, archived on Zenodo with a DOI — and already in real use here at
PoliMi: I'm building grant proposals on it, with irradiated materials, coated
claddings and advanced fuels in view for the next Euratom calls; three students are
developing the code with me; and it's used as a support tool in our didactic
activities. Where it goes next: sharper contact through dolfinx-mpc — per-facet pressure for
axial ridging — coupling to the fuel codes SCIANTIX and Mérope for fission gas and
microstructure, a monolithic phase-field Newton for spontaneous nucleation, and
pushing the differentiable foundation all the way to constitutive-law discovery
through the coupled adjoint. And if you'd like to see it run today, come to the
software demonstration at the poster session — I'll take a coupled case end to end
on my laptop and tune it while it solves. Thank you — I'm happy to take questions."

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
  pressure proportional to penetration, equal-and-opposite tractions, applied as a
  fixed-point load: the pressure is frozen from the previous iterate and the staggered
  loop drives it to convergence (explicit for now — not yet a consistent Newton tangent
  for contact). Verified against the analytical Lamé interference-fit to 3.5% (backup
  slide), with the stress state confirmed plane-stress.
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
- *The abstract mentions symbolic regression / automated discovery of
  constitutive laws — where is that?* Preliminary, and I show the first step at
  the demo table: the time integrator is differentiable, so forward-mode AD
  through the implicit solver gives exact parameter sensitivities, and
  gradient-based least squares recovers Norton creep parameters from noisy
  relaxation data (exponent to ~2% in ten Gauss-Newton iterations, ~1 s).
  That's parametric identification; EUCLID-style sparse regression over a
  library of candidate constitutive models (Flaschel et al. 2022) is the roadmap — the
  differentiable foundation is the part you saw in this talk. Come to the demo
  for the live run.
- *Can it do fuel performance — burnup, swelling, gap closure over an irradiation?*
  That is slide 6 — burnup accumulates as a per-material state field,
  swelling is an eigenstrain callable that reads it (handled like any other
  material property — no solver change), and an 1800-day rod history at 20 kW/m closes the 65 µm gap at
  ~21 MWd/kgU; the contact pressure emerges and then saturates near 25 MPa as
  cladding irradiation creep relaxes it; the mean burnup is verified against the
  closed form to 2e-8. The peak fuel temperature *drops* from 1126 K to ~976 K on
  contact — the thermal feedback, live in the figure. The end-of-life stress state
  (compressive hoop at the fuel surface) is in the backup.

## Backup slides (do not present; flip to on a question)
1. staggered solver loop · 2. constitutive-model table · 3. plasticity verified
(J2 + crystal) · 4. coupled thermo-mechanics verified · 5. automated creep-law
discovery (EUCLID-style) · 6. penalty contact verified vs Lamé (the old PCMI
slide) · 7. integral-rod end-of-life stress state · 8. SEN shear, Miehe (2010) ·
9. hero case verified vs McClenny (2022).

## Delivery notes
- The AD slide (3), PCMI (6) and the hero case (8) are what the audience
  remembers. Slow down on those; speed through 2 and 7 if you must.
- On slide 6: lead with the gap closing and the *emergent* pressure, then land the two
  payoffs — the pressure *creep-saturates* (cladding creep), and the peak fuel
  temperature *drops* on contact. Point at the gap-closure/pressure panel, then the
  temperature curve. Lamé verification ("3.5%") is in the backup if asked.
- Do not read the code aloud line by line — point at `ufl.diff(psi, F)` and say it
  once. On slide 6, be clear the contact is an explicit penalty for now — don't claim
  a consistent / AD tangent for it (the AD tangent is for the constitutive stress).
- Practised cold, with slide 6 added this runs ~9:45. The compressible asides (cluster
  line on 5, teaching aside on 7) buy ~20 s of slack if the committee clock is
  tight on the 10-minute format.
