# Cluster Dynamics Model in Z3ST: Theoretical Background

This document describes the mathematical formulation and numerical implementation of the **Cluster Dynamics (CD)** model inside the **Z3ST** framework. The model simulates the evolution of the size distribution of defect clusters in a material subjected to irradiation or thermal stress.

## 1. Physical Context

Cluster Dynamics is a mesoscopic modelling technique that describes the population of defect clusters (e.g. vacancy aggregates, interstitial atoms, gas bubbles, dislocation loops). Unlike molecular dynamics, the system is not resolved in physical space but in **size space** $n$, where $n$ is the number of monomers that compose a cluster. The state variable is the cluster density

$$
c(n,t) \quad \big[\text{clusters} \cdot \text{m}^{-3} \cdot \text{size}^{-1}\big],
$$

representing the number of clusters of size $n$ at time $t$.

## 2. From the Master Equation to the Continuum

### 2.1 Discrete Master Equation

The fundamental process is monomer exchange. A cluster of size $n$ can grow by capturing a monomer (condensation rate $\beta_n$) or shrink by losing one (evaporation rate $\alpha_n$). The evolution of the discrete density $c_n(t)$ obeys

$$
\frac{\partial c_n}{\partial t} = J_{n-1} - J_n,
\qquad
J_n = \beta_n c_n - \alpha_{n+1} c_{n+1}.
$$

### 2.2 Fokker–Planck Approximation

Through a Taylor/Kramers–Moyal expansion (valid when $n$ varies smoothly over large clusters) the discrete system reduces to the continuous **Fokker–Planck equation**

$$
\frac{\partial c}{\partial t} = -\frac{\partial (v\,c)}{\partial n} + \frac{\partial^2 (D\,c)}{\partial n^2}.
$$

In **Z3ST** the coefficients are assumed constant with respect to $n$ (configurable through `input.yaml`), so the governing PDE reduces to the classical **advection–diffusion equation in size space**

$$
\frac{\partial c}{\partial t} + v\,\frac{\partial c}{\partial n} = D\,\frac{\partial^2 c}{\partial n^2}.
$$

### 2.3 Physical Parameters

* **Advection (drift) velocity $v$**: net cluster growth rate, $v = \beta - \alpha$.
  Set in `input.yaml` via `cluster.advection_velocity` (default `1.0`).
* **Diffusion coefficient $D$**: stochastic size fluctuations, $D = \tfrac{1}{2}(\beta + \alpha)$.
  Set in `input.yaml` via `cluster.diffusion_coefficient` (default `0.5`).

The mesh coordinate `x[0]` plays the role of the cluster size $n$.

## 3. Initial Conditions

The module [cluster_dynamic_model.py](../../models/cluster_dynamic_model.py) supports two initial-condition types, configured under `cluster.initial_condition` in `input.yaml`:

### 3.1 Constant IC

A uniform value $c_0$ is prescribed on a named region defined in `geometry.yaml`:

$$
c(n, 0) = c_0 \qquad \text{on } \Omega_{\text{IC}}.
$$

Useful for modelling a uniformly pre-populated size band (e.g. a monodisperse embryo population).

### 3.2 Gaussian IC

An analytical Gaussian profile is interpolated across the size domain:

$$
c(n, 0) = A\,\exp\!\left[-\frac{(n - \bar n)^2}{2\sigma^2}\right],
$$

with peak position $\bar n$ (`mean`), width $\sigma$ (`std_dev`), and normalization target (`amplitude`). This form is the fundamental solution of pure diffusion in size space and the natural attractor predicted by the Zeldovich–Frenkel nucleation theory near the critical size.

The test case `U_cluster_dynamics_test/input.yaml` uses `mean = 1.0`, `std_dev = 1.0`, `amplitude = 1000.0`.

## 4. Mass Conservation

The fundamental physical invariant is the **first moment** of the distribution, representing the total number of monomers bound in clusters:

$$
C_\text{tot}(t) = \int_{n_{\min}}^{n_{\max}} n\,c(n, t)\,dn.
$$

Because a straight continuous-Galerkin discretisation of a convection-dominated PDE does not guarantee discrete conservation of this moment, **Z3ST** enforces it through two mechanisms:

1. **Initialization** — in `set_cluster_initial_conditions` the initial profile is rescaled so that $\int n\,c(n,0)\,dn$ equals the user-defined target `C_tot_target`.
2. **Per-step renormalization** — at the end of each time step in `_cluster_step` the density is rescaled:

$$
c_\text{new} \leftarrow c_\text{new} \cdot \frac{C_\text{tot,target}}{C_\text{tot,current}}.
$$

This multiplicative correction keeps the total mass constant while preserving the *shape* of the distribution produced by the PDE solve.

## 5. Numerical Implementation

The solver lives in [_cluster_step](../../core/solver.py) (method of the `Solver` class). It uses a **Discontinuous Galerkin (DG)** formulation with implicit Euler time integration.

### 5.1 Time Discretization

Backward (implicit) Euler with time step $\Delta t$:

$$
\frac{c^{n+1} - c^{n}}{\Delta t} + v\,\frac{\partial c^{n+1}}{\partial n} = D\,\frac{\partial^2 c^{n+1}}{\partial n^2}.
$$

### 5.2 DG Weak Form

Let $V_c$ be a DG finite-element space on the size-space mesh, $\phi$ a test function, $\mathbf{n}$ the facet normal, and denote by $\{\!\{\cdot\}\!\}$ and $[\![\cdot]\!]$ the facet average and jump operators.

**Time-derivative (mass) term:**

$$
\int_\Omega \frac{c^{n+1}}{\Delta t}\,\phi\;dn = \int_\Omega \frac{c^{n}}{\Delta t}\,\phi\;dn.
$$

**Advection — Upwind flux:**

$$
-\int_\Omega c^{n+1}\,v\,\partial_n\phi\;dn
+ \int_{\mathcal{F}^{\mathrm{int}}}\!\!\Big( \{\!\{c^{n+1}v\,n_0\}\!\}[\![\phi]\!] + \tfrac{1}{2}|v_n^+|\,[\![c^{n+1}]\!][\![\phi]\!] \Big)\,dS
+ \int_{\partial\Omega_{\mathrm{out}}} v_n\,c^{n+1}\phi\;ds.
$$

The boundary term is activated only where $v\cdot n > 0$ (outflow), giving a natural outflow condition at $n_{\max}$ and a no-flux inflow condition at $n_{\min}$.

**Diffusion — Symmetric Interior Penalty (SIPG):**

$$
\int_\Omega D\,\partial_n c^{n+1}\,\partial_n\phi\;dn
- \int_{\mathcal{F}^{\mathrm{int}}} D\,\{\!\{\partial_n c^{n+1}\}\!\}[\![\phi\,n_0]\!]\;dS
- \int_{\mathcal{F}^{\mathrm{int}}} D\,\{\!\{\partial_n\phi\}\!\}[\![c^{n+1}\,n_0]\!]\;dS
+ \int_{\mathcal{F}^{\mathrm{int}}} \frac{D\,\gamma}{\langle h\rangle}\,[\![c^{n+1}]\!][\![\phi]\!]\;dS.
$$

The penalty parameter is $\gamma = 10.0$ (sufficient for P1 elements) and $\langle h \rangle$ is the average facet length.

### 5.3 Linear Solver

The assembled system is solved with PETSc using GMRES preconditioned by ILU (`ksp_rtol = 1e-12`, `ksp_atol = 1e-15`, `ksp_max_it = 1000`).

## 6. Stability: Péclet Number Diagnostics

The advection-to-diffusion ratio is monitored at each step via the **cell Péclet number**

$$
\mathrm{Pe} = \frac{|v|\,h}{2D},
$$

where $h = L_\text{domain}/N_\text{cells}$. When $\mathrm{Pe} > 1$ the system is advection-dominated; the DG upwind flux provides the required stabilization without any extra user action. The value is printed at every solve for diagnostic purposes.

## 7. Boundary Conditions in Size Space

Boundary conditions in $n$-space carry a direct physical meaning:

* **Inflow boundary at $n = n_{\min}$** (no-flux by default): the small-cluster population is not artificially replenished. A nucleation source term $G(n,t)$ localised near $n = 1$ may be added as a future extension.
* **Outflow boundary at $n = n_{\max}$** (upwind outflow): clusters that grow past the domain boundary are removed from the size window. Physically, this models transfer to larger aggregates not represented on the mesh.

## 8. Code Structure

* [cluster_dynamic_model.py](../../models/cluster_dynamic_model.py) — sets up the physical parameters ($v$, $D$) and initial conditions (constant or Gaussian), and computes the initial mass target `C_tot_target`.
* [solver.py :: _cluster_step](../../core/solver.py) — assembles the DG bilinear/linear forms (upwind advection + SIPG diffusion + implicit Euler), solves the linear system with GMRES/ILU, enforces mass conservation by renormalization, and reports the Péclet number and diagnostics.
* [input.yaml](input.yaml) — selects the cluster model (`models.cluster: true`), sets $v$ and $D$, and configures the initial-condition type and parameters.

## 9. Current Limitations and Planned Extensions

* **Constant coefficients** — $v$ and $D$ are $n$-independent. Realistic CD uses $v(n), D(n) \propto n^{1/3}$ (diffusion-limited capture) or stronger size dependencies.
* **No source term** — irradiation-induced nucleation $G(n,t)$ and absorption sinks $-S(n)\,c$ (grain boundaries, dislocations) are not present. They would enter as additional RHS contributions.
* **1D in size space only** — there is no spatial transport: the size distribution at each physical location evolves independently. Coupling to a spatial reaction–diffusion field is a future direction.
* **Global mass correction** — mass conservation is enforced via a multiplicative rescaling rather than by a locally conservative flux formulation. A fully conservative LDG or flux-corrected transport scheme is the natural next step.

---
*Author: Giovanni Zullo*
*Version: 0.2.0 (2026)*
