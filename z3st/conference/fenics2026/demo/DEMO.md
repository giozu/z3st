# Z3ST live software demonstration — run-sheet

FEniCS 2026, software-demonstration slot (runs during the poster session, on your
own laptop, to small groups). Paris, 17–19 June 2026.

This is the document you rehearse from. It is built to *win a prize* — see the
strategy section first, then the segmented run-sheet, then the recovery plan and
the pre-flight checklist.

---

## 0 · Prize strategy (read this once, internalise it)

The FEniCS 2026 awards (confirmed by the final-info email and the website):
- **Best Presentation by a PhD Candidate** — $500.
- **Best Presentation by a Postdoctoral Researcher** — $500.
- **Best FEniCS 2026 Poster** — $500.
- **Nate Sime's Exceptional Visualization** — surprise prize; "most beautiful
  according to the committee's subjective sense of aesthetics", composed from
  data generated primarily with FEniCS. The case-14 animation qualifies.

Judging is on a 0–7 scale across two axes:
1. **Comprehension & content** — clear background, transparent strategy and results, articulated impact.
2. **Engagement & communication** — clarity, organisation, *enthusiasm*, thoughtful answers to questions, holding attention.

What this means for you, concretely:

| Lever | How this demo pulls it |
|---|---|
| **Visualization award** (most winnable here) | The case-14 crack-propagation animation in ParaView is the centrepiece. Smooth, coloured by damage, scrubbable. Pre-baked so it is always perfect. |
| **Engagement** | The demo is *interactive*: you run things live, and you change a parameter mid-solve (hot reload) and the audience sees the solver react. People remember the thing that moved. |
| **Comprehension** | Every segment ties back to one sentence: *"new physics is a new energy, differentiated automatically."* Repeat it. |
| **Impact** | Land the validation point: case 14 reproduces McClenny et al. (JNM 2022) UO2 thermal-shock cracking. Not a toy — a published experiment. |

Two non-obvious wins:
- **Demos cycle.** A small group leaves, a new one arrives. Your core loop must be tight (~4 min) and repeatable from a clean state. Practise the *restart*, not just the run.
- **Judges wander the poster session.** They may catch any 4 minutes of you. Every loop must contain the animation and the one-sentence pitch — never bury the best thing waiting for a "later" that a judge will not stay for.
- **The abstract is a promise.** Visitors who read it expect *crystal plasticity*
  and *automated constitutive-law discovery* — segments K and M deliver both on
  demand. Offer them unprompted to anyone who mentions the abstract.

---

## 1 · The one-sentence pitch (say it in the first 20 seconds, every loop)

> "Z3ST turns FEniCSx into a YAML-driven solver for coupled thermo-mechanics and
> fracture — you write the constitutive physics as an energy, and automatic
> differentiation gives you the stress and the tangent for free."

---

## 2 · Core loop (~4–5 min, repeatable)

Run `./run_demo.sh` and step through the segments with the Enter key, so you
control the pace while you talk. Talking points below.

### A. Hook + "it just works" — 60 s
- Run the **1D bar** (`teaching/01_1D`) live: a line mesh, `regime: 1d`.
  - *"One YAML file. Steel bar, pulled. Analytical answer u(L)=PL/E. Watch."*
  - Point at the **PASS** line: `u_xL → rel err ~1e-15 → PASS`.
- Then the **same bar in full 3D** (`teaching/01_3D`), identical displacement.
  - *"Same model, I changed one flag — `regime: 3d`. No re-meshing, no new code. That is the design: one driver, many regimes."*

### B. The core idea — automatic differentiation — 60 s
- Open `z3st/models/mechanical_model.py` at the hyperelastic block (the launcher
  prints the exact lines). Point at:
  ```python
  P   = ufl.diff(psi, F)        # 1st Piola-Kirchhoff stress = d psi / d F
  Jac = ufl.derivative(R, u)    # consistent tangent — also AD
  ```
  - *"This one line is the stress. I differentiate the energy. The Newton tangent
    is another `ufl.derivative`. I never hand-derive a Jacobian — that is why a new
    material is just a new energy."*

### C. Coupled thermo-mechanics + the hot-reload "wow" — 90 s
- Start the **coupled slab** (`verification/thermal/thin_slab_neumann_2D`): thermal + mechanical,
  staggered. It solves in a few seconds and prints the staggered iteration count.
  - *"Now it is coupled: heat conduction drives thermal strain, which drives stress.
    Staggered loop with adaptive relaxation."*
- **The wow:** the launcher pauses and invites you to edit `input.yaml` live
  (e.g. drop `relax_u` or tighten `stag_tol`) and re-run — Z3ST hot-reloads
  allow-listed parameters *without restarting*.
  - *"I can retune the solver while it runs — tolerances, relaxation — and it picks
    up the change at the next step. Great for exploring a hard nonlinear case."*

### D. The showpiece — crack propagation in ParaView — 60 s  ← the prize shot
- Open the **pre-baked case-14 animation** (`./open_paraview.sh`): UO2 pellet,
  cold-contact wedge, hoop-stress ring, a radial crack growing as you scrub time.
  - *"Fully coupled: temperature → elastic strain → phase-field damage. AT1, Amor
    split, hybrid. This reproduces McClenny 2022 — a real thermal-shock experiment."*
  - Scrub the timeline so the crack visibly advances. Rotate. Let it breathe.

### E. Close — 30 s
- *"It is open, Apache-2.0, on GitHub, archived on Zenodo with a DOI, documented.
  Your case is a YAML file. Clone it tonight."*
- Hand them the handout / point at the QR codes on the attract screen.
  - `github.com/giozu/z3st` · `giozu.github.io/z3st` · DOI `10.5281/zenodo.17748028`

### P. (optional) Multi-body: pellet–cladding contact, verified — 60 s
*Not in the default tight loop — pull it up for a fuel-interested visitor or a
judge who lingers (`./run_demo.sh P`; it opens the baked figures). It pairs with
the oral's PCMI slide.*
- *"A fuel pellet heats, expands, closes the gap, and contacts the cladding — all
  coupled. Contact is a penalty: pressure proportional to penetration. Nothing is
  prescribed — the pressure **emerges**, and the cladding is driven outward, load
  transfer. And it feeds back: contact raises the gap conductance, so the fuel
  cools the moment it touches."*
- Impact: *"and it's **verified** — 3.5% against the analytical Lamé interference
  fit, the stress state confirmed plane-stress."*
- Tie back to the one sentence: *"the penalty tangent is the same AD path —
  `ufl.derivative` — no hand-coded contact Jacobian."*
- **The burnup beat (the fuel-performance closer):** open
  `baked/pcmi_burnup_curves.png`. *"That was a power ramp. Here is the same
  physics over an **1800-day irradiation**: burnup accumulates as a material
  state field, swelling is an eigenstrain that reads it — handled like any other
  material property, no solver change — and the gap closes at ~21 MWd/kgU. Watch the right panel:
  the moment contact engages, the gap conductance jumps and the **peak fuel
  temperature drops**. That is the two-way coupling, visible in one figure. And
  the burnup itself is verified against the closed form to 2e-8."*
- Optional live (during a longer chat): run the case and watch the gap close and
  contact switch on in the streamed `[contact] … CLOSED` lines.

### K. (optional) Crystal plasticity — the demo abstract's first headline — 90 s
*Pull it up for anyone who read the abstract, a plasticity person, or a judge
(`./run_demo.sh K`; runs live in ~11 s, baked fallback `baked/cp_stress_strain.png`).*
- Run `verification/plasticity/crystal_single_grain` live: single crystal, one FCC slip system
  (111)[0-11], power-law viscoplasticity, backward Euler with history variables
  in quadrature spaces.
- *"The slip-rate derivative through the Schmid tensor is the Jacobian nobody
  wants to derive by hand. Here it comes from `ufl.diff` — exactly. That's why
  Newton converges quadratically: 2 iterations per step."*
- Point at the stress-strain curve: elastic, yield at τ = g₀, saturation at the
  **semi-analytical σ_sat to 3.4%**.
- Tie back: *"the constitutive law is ~10 lines of Python — same AD path as
  everything else."*

### M. (optional) Toward constitutive-law discovery — identification from data — 90 s
*The demo abstract's second headline (Flaschel et al. 2022 — EUCLID), shown as
the first preliminary step (`./run_demo.sh M`; runs live in ~2 s, bakes
`baked/creep_identification.png`).*
- *"Same AD idea, pointed the other way: differentiate the **solver** with
  respect to the **material parameters**, and you can learn the law from data."*
- Run `identify_creep.py` live: Norton creep relaxation (the verified case),
  51 synthetic data points with 2% noise, forward-mode AD propagated through
  every implicit backward-Euler step, Gauss-Newton least squares.
- Land the numbers: *"the creep exponent comes back to ~2%, the rate to ~1%,
  in ten iterations — about a second."*
- Be honest and forward-pointing: *"this is **parametric identification**
  today; EUCLID-style sparse-regression discovery over a library of candidate
  energies is the roadmap — the differentiable foundation you just saw is the
  hard part."* (If asked: independent implementation; the EUCLID codes are
  GPL-3.0 and nothing here derives from them.)

---

## 3 · Deep-dive cards (only if a knowledgeable person stays)

Have these ready; do not volunteer them in the core loop.
- **Multiscale:** cluster dynamics solves advection–diffusion in *cluster-size*
  space and couples back to local T and strain; phase-field length scale `lc`;
  gap conductance as an interface BC.
- **Solver:** staggered with per-field adaptive relaxation (grow on shrinking
  residual, back off on oscillation); PETSc + MPI under dolfinx; MUMPS / HYPRE.
- **Damage degrades both elastic and thermal-stress terms** via a consistent
  `g(D)` — without it the thermal-shock problem is intractable. This is the
  detail that makes case 14 actually converge.
- **Verification:** ~50 cases, each with a `non-regression.py`, run on every
  commit in CI.

---

## 4 · Recovery plan (assume something will go wrong)

| If… | Do this |
|---|---|
| A live run hangs or errors | Ctrl-C, say *"let me show you the baked result"* and open the pre-rendered PNG/animation. Never debug in front of a group. |
| dolfinx import fails | The env was not activated — `conda activate z3st11` (0.11 primary; `z3st` 0.10 is the fallback). `preflight.sh` checks this; run it before the session. |
| ParaView is slow / wrong screen | Fall back to the pre-rendered PNG sequence in `demo/baked/` — `open_paraview.sh --baked` just opens the images. |
| Projector resolution is bad | Pre-set ParaView and terminal font large (see checklist). Use the PNG sequence, not interactive, if rotation stutters. |
| No internet | Everything here runs **offline**. The only online things are the GitHub/DOI links — those live on the handout and the attract screen, not the laptop. |
| Someone asks something you do not know | *"Good question — I am not certain; here is how I would find out."* Judges score *thoughtful* answers, not bluffing. |

**Golden rule:** the baked artifacts in `demo/baked/` mean you can deliver the
entire pitch with zero live computation if you must. Live runs are upside, not a
dependency.

---

## 5 · Pre-flight checklist (the morning of, and again before your slot)

Run `./preflight.sh` first — it verifies the env, the cases, and the baked
artifacts. Then, by hand:

- [ ] `conda activate z3st11` in the demo terminal; leave it activated (`z3st` 0.10 is the fallback).
- [ ] Terminal font size large (≥ 18 pt); dark-on-light or high-contrast theme.
- [ ] Laptop on mains power; screen-sleep and notifications **off**.
- [ ] ParaView opens the baked case-14 series and the animation plays smoothly.
- [ ] `teaching/01_1D`, `teaching/01_3D`, `verification/thermal/thin_slab_neumann_2D` each run clean once (warms the dolfinx import cache too).
- [ ] `demo/baked/` contains the fallback PNGs (incl. `pcmi_curves.png`, `pcmi_verification.png` for segment P, `cp_stress_strain.png` for segment K, `creep_identification.png` for segment M).
- [ ] `verification/plasticity/crystal_single_grain` runs clean once (segment K live run, ~11 s) and `python3 identify_creep.py` converges (segment M, ~2 s) — both checked by `preflight.sh`.
- [ ] Open `attract.html` once (`./attract.sh`) and confirm the loop plays and the QR codes scan.
- [ ] `../handout/handout.pdf` printed (a small stack to hand out); QR codes scan.
- [ ] Editor open on `mechanical_model.py` (AD lines) and `verification/thermal/thin_slab_neumann_2D/input.yaml` (for the hot-reload edit) in separate tabs.
- [ ] No poster was submitted — the attract screen and the handout are the
      stand-ins; both carry the QR codes (repo, docs, DOI) and your email.
- [ ] A few business cards / a printed one-pager to hand out.
- [ ] Water. You will talk for two hours straight in loops.

---

## 6 · Files in this folder

```
demo/
├── DEMO.md            this run-sheet
├── run_demo.sh        interactive, segmented launcher (A–E core + P/K/M optional)
├── preflight.sh       environment + cases + artifacts check
├── open_paraview.sh   opens the case-14 crack animation (live or --baked)
├── paraview_case14.py pvpython script: builds the damage animation / PNG sequence
├── identify_creep.py  segment M: gradient-based creep-law identification (own AD)
├── attract.html       idle-table attract loop (crack animation + pitch + QR codes)
├── attract.sh         opens the attract loop full-screen in a browser
├── qr/                QR codes for the attract loop (repo, docs, DOI) + make_qr.sh
└── baked/             pre-rendered fallback visuals (generated before the conference)

(the printed leave-behind lives one level up in ../handout/handout.pdf)
```

## 7. The idle table — attract loop

When no one is at your table, put the laptop into the **attract loop**: it plays
the case-14 crack animation full-screen with the one-sentence pitch and QR codes.
Motion on a screen pulls the cycling poster-session crowd, and it puts your
visualization in front of every passing judge even while you are talking to
someone else.

```bash
./open_paraview.sh --render   # once, to bake the frames the loop plays
./attract.sh                  # opens attract.html; press F11 / click for full screen
```

It is pure offline HTML — no internet, no ParaView, no dependencies. If you prefer
a second screen, run the attract loop on one and your live terminal on the other.
