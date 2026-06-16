# Z3ST live demo — spoken script + screen setup

Delivered in English (Paris). Action cues in [brackets]. Core loop ≈ 4–5 min,
repeatable. Keep this open on your phone or a printed sheet beside the laptop.

---

## 0 · Screen setup — what to have open BEFORE you start

Activate the env once, in the terminal you'll demo from:

```
conda activate z3st11          # 0.11 (migrated). Fallback: conda activate z3st (0.10)
cd ~/z3st/z3st/conference/fenics2026/demo
```

Then have exactly these windows open:

| Window | What | Used in |
|---|---|---|
| **Terminal** (font ≥ 18 pt) | you run `./run_demo.sh` here and press **Enter** to step through A→E | A, C, E |
| **Editor — tab 1** | `z3st/models/mechanical_model.py`, scrolled to **~line 635–656** (the neo-Hookean `psi` and `P = ufl.diff(psi, F_def)`). The launcher also prints the exact line numbers. | B |
| **Editor — tab 2** | `z3st/cases/verification/thermal/thin_slab_neumann_2D/input.yaml` (you'll edit this live) | C |
| **ParaView** | DON'T pre-open — it launches itself when you run `./open_paraview.sh` | D |
| **Browser (attract loop)** | separate, for the idle table: `./attract.sh`. Put it on a 2nd screen if you have one. | idle |

Single screen → keep Terminal + Editor side by side (or Alt-Tab); ParaView pops up for D.
Two screens → Terminal + Editor on the laptop, attract loop on the external.

Per segment, what you actually touch:
**A** terminal · **B** editor tab 1 (point) · **C** terminal → editor tab 2 (edit+save) → terminal · **D** ParaView · **E** terminal + handout.

---

## 1 · Opening (say it EVERY new group, first 20 s)

> "Hi! Let me show you Z3ST in about four minutes. The one idea: **Z3ST turns
> FEniCSx into a YAML-driven solver for coupled thermo-mechanics and fracture —
> you write the physics as an energy, and automatic differentiation gives you the
> stress and the tangent for free.** Let me show you what that means."

---

## 2 · Core loop

### A — Hook: "it just works" (≈60 s)
[Terminal: `./run_demo.sh`, Enter on the 1D bar]
> "A steel bar in tension. One YAML file — geometry, material, boundary conditions.
> The analytical answer is u(L) = PL/E… and look —" [point at the PASS line] "— it
> matches to fifteen digits.
> Now watch: **same model, I change one flag — `regime: 1d` to `regime: 3d`.**"
> [Enter → 3D runs] "Same bar, full 3D, identical displacement. No re-meshing, no new
> code. That's the design: **one driver, many regimes.**"

### B — The core idea: automatic differentiation (≈60 s)
[Switch to Editor tab 1 — mechanical_model.py, the psi / ufl.diff lines]
> "Here's the heart of it. A constitutive model is a **strain-energy density** —
> this is neo-Hookean. And the stress?" [point at `P = ufl.diff(psi, F_def)`] "**One
> line.** I differentiate the energy with respect to the deformation gradient. The
> Newton tangent is another `ufl.derivative` — also automatic.
> I **never** hand-derive a Jacobian. A new material is just a new energy —
> elasticity, plasticity, crystal plasticity, phase-field, all the same path."

### C — Coupled physics + the hot-reload "wow" (≈90 s)
[Terminal: Enter → the coupled slab runs]
> "Now let's couple it: heat conduction drives thermal strain, which drives stress —
> staggered loop with adaptive relaxation. A few seconds, there's the iteration count."
[Switch to Editor tab 2 — input.yaml; change e.g. `relax_u` or `stag_tol`, SAVE; back to terminal]
> "And my favourite part — I can **retune the solver while it runs.** Tolerances,
> relaxation — Z3ST hot-reloads them at the next step, no restart. Great for pushing
> through a stiff nonlinear case interactively."

### D — The showpiece: crack propagation (≈60 s)  ← the prize shot
[Terminal: `./open_paraview.sh` → ParaView opens; scrub the timeline, rotate]
> "Everything together: thermal-shock cracking of a **UO₂ nuclear fuel pellet.** A
> cold contact wedge cools the rim, sets up a tensile hoop-stress ring — and watch a
> **radial crack grow** as I scrub through time." [scrub slowly]
> "Fully coupled: temperature → elastic strain → phase-field damage. AT1, Amor split,
> hybrid. And it's not a toy — **it reproduces McClenny et al. 2022, a published
> thermal-shock experiment.** Same crack pattern. All from one `input.yaml`."

### E — Close (≈30 s)
[Hand out the handout / point at the QR codes on the attract screen]
> "Open source, Apache-2.0, on GitHub, archived on Zenodo with a DOI, documented and
> verified — about fifty cases checked against analytical solutions on every commit.
> **Your problem is a YAML file — clone it tonight.** I'm Giovanni, Politecnico di
> Milano. Questions?"

---

## 3 · Restart (between groups)
> "Hi — perfect timing, let me run it from the top." [reset, start at A]

## 4 · Add-ons (only if they read the abstract / are expert)  [one line + run the segment]
- **Crystal plasticity** [`./run_demo.sh K`]: "The slip-rate Jacobian — the one nobody
  wants to derive by hand — comes straight from `ufl.diff`. Newton converges in two
  iterations."
- **Constitutive-law discovery** [`./run_demo.sh M`]: "Same AD, backwards: differentiate
  the solver with respect to the material parameters, and you **learn the law from
  data** — creep exponent recovered to ~2% in about a second."
- **PCMI / fuel** [`./run_demo.sh P`]: "Pellet swells, closes the gap, contacts the
  cladding — the contact pressure **emerges**, nothing prescribed, verified to 3.5%
  against analytical Lamé."

---

## 5 · Delivery tips
- Recurring refrain: **"new physics = a new energy, differentiated automatically."**
- Slow down on **D** (the crack) and **C** (the live hot-reload) — that's what they remember.
- If anything stumbles live: *"let me show you the baked result"* → open the baked image.
  Never debug in front of a group.
- Smile, let the crack "breathe", show enthusiasm — it's scored.
- Recovery: `./open_paraview.sh --baked` (baked loop) · `conda activate z3st` (0.10 fallback).
