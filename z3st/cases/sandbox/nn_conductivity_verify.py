# --.. ..- .-.. .-.. ---
# Verification case: NN thermal conductivity k(T) via dolfinx-external-operator (dolfinx 0.11).
#
# Standalone sandbox prototype (does NOT touch z3st/core). Steady 1D conduction in a
# slab, no source, Dirichlet T on the two x-faces, insulated y-faces. Conductivity is a
# neural network k = NN(T) trained to fit the linear law k(T) = a + b T, for which the
# continuum solution is closed-form. We then solve the genuinely non-linear problem with
# the external-operator Newton (Tier 2) and certify three things (design note s6):
#
#   (a) NN-T matches the ANALYTIC T(x)                         -> physics + fit are correct
#   (b) NN-T matches dolfinx NewtonSolver on the SYMBOLIC k(T) -> NN route == UFL-card route
#   (c) Taylor-remainder test on the residual: slope ~ 2       -> tangent dk/dT is correct
#
# Run (z3st env): python3 nn_conductivity_verify.py
# Writes knet.pt (trained weights) and nn_conductivity_verify.png (T profiles) here.

import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from mpi4py import MPI
from petsc4py import PETSc

import basix
import dolfinx
import ufl
from dolfinx import fem, mesh
from dolfinx.fem.petsc import (
    apply_lifting,
    assemble_matrix,
    assemble_vector,
    set_bc,
    NonlinearProblem,
)

import torch

from dolfinx_external_operator import (
    FEMExternalOperator,
    evaluate_external_operators,
    evaluate_operands,
    replace_external_operators,
)

HERE = os.path.dirname(os.path.abspath(__file__))
comm = MPI.COMM_WORLD
print(f"dolfinx {dolfinx.__version__} | torch {torch.__version__}")

# --. problem definition --..
A_COEF, B_COEF = 1.5, 0.01      # true law k(T) = a + b T  [W/m/K]
T_L, T_R = 400.0, 1200.0        # Dirichlet temperatures [K]
W = 1.0                         # slab width [m]
T_MIN, T_MAX = min(T_L, T_R), max(T_L, T_R)
T0, TSCALE = 0.5 * (T_MIN + T_MAX), 0.5 * (T_MAX - T_MIN)   # normalisation


def k_true(T):
    return A_COEF + B_COEF * T


def T_analytic(x):
    """Closed-form steady T(x) for k(T)=a+bT, no source, Dirichlet ends.
    From (a+bT) dT/dx = const: a T + b T^2/2 = C x + D."""
    C = (A_COEF * (T_R - T_L) + 0.5 * B_COEF * (T_R**2 - T_L**2)) / W
    D = A_COEF * T_L + 0.5 * B_COEF * T_L**2
    return (-A_COEF + np.sqrt(A_COEF**2 + 2.0 * B_COEF * (C * x + D))) / B_COEF


# --.. ..- .-.. .-.. ---  1. train (or load) the network  --.. ..- .-.. .-.. ---
torch.manual_seed(0)
net = torch.nn.Sequential(
    torch.nn.Linear(1, 32), torch.nn.Tanh(),
    torch.nn.Linear(32, 32), torch.nn.Tanh(),
    torch.nn.Linear(32, 1),                       # linear output; tanh hidden -> smooth k'
)
weights_path = os.path.join(HERE, "knet.pt")

if os.path.exists(weights_path):
    net.load_state_dict(torch.load(weights_path))
    print(f"loaded weights from {weights_path}")
else:
    T_train = torch.linspace(T_MIN, T_MAX, 400).unsqueeze(-1)
    k_train = A_COEF + B_COEF * T_train
    Tn = (T_train - T0) / TSCALE
    opt = torch.optim.Adam(net.parameters(), lr=1e-2)
    for step in range(4000):
        opt.zero_grad()
        loss = torch.mean((net(Tn) - k_train) ** 2)
        loss.backward()
        opt.step()
    torch.save(net.state_dict(), weights_path)
    print(f"trained net, final MSE = {loss.item():.3e}, saved -> {weights_path}")
net.eval()
net.double()   # float64 inference: needed for a clean Taylor test (float32 noise floor ~1e-5)

# fit quality
with torch.no_grad():
    Tg = torch.linspace(T_MIN, T_MAX, 200, dtype=torch.float64).unsqueeze(-1)
    kg = net((Tg - T0) / TSCALE)
    rel = (kg - k_true(Tg)).abs().max().item() / k_true(T_MAX)
print(f"NN fit: max relative error over [{T_MIN:.0f},{T_MAX:.0f}] K = {rel:.2e}")


def _k_and_dk(T_np):
    """Vectorised over Gauss points: T (ndarray) -> (k, dk/dT) ndarrays."""
    t = torch.tensor((T_np.ravel() - T0) / TSCALE, dtype=torch.float64, requires_grad=True)
    k = net(t.unsqueeze(-1)).squeeze(-1)
    (dk_dtn,) = torch.autograd.grad(k, t, torch.ones_like(k))
    k_np = k.detach().numpy().astype(T_np.dtype)
    dk_np = (dk_dtn.detach().numpy() / TSCALE).astype(T_np.dtype)   # chain rule d/dT
    return k_np, dk_np


def k_external(derivatives):
    """Dispatch by derivative multi-index: (0,) -> k, (1,) -> dk/dT."""
    if derivatives == (0,):
        return lambda T_np: _k_and_dk(T_np)[0]
    if derivatives == (1,):
        return lambda T_np: _k_and_dk(T_np)[1]
    raise NotImplementedError(f"derivative {derivatives} not implemented")


# --.. ..- .-.. .-.. ---  2. mesh, spaces, BCs  --.. ..- .-.. .-.. ---
domain = mesh.create_rectangle(comm, [[0.0, 0.0], [W, 0.1]], [40, 4], mesh.CellType.triangle)
tdim = domain.topology.dim
V = fem.functionspace(domain, ("Lagrange", 1))
v = ufl.TestFunction(V)
dT = ufl.TrialFunction(V)


def _linear_guess(x):
    return T_L + (T_R - T_L) * x[0] / W


dofs_L = fem.locate_dofs_geometrical(V, lambda x: np.isclose(x[0], 0.0))
dofs_R = fem.locate_dofs_geometrical(V, lambda x: np.isclose(x[0], W))
bcs = [
    fem.dirichletbc(PETSc.ScalarType(T_L), dofs_L, V),
    fem.dirichletbc(PETSc.ScalarType(T_R), dofs_R, V),
]

deg = 2
Qe = basix.ufl.quadrature_element(domain.basix_cell(), value_shape=(), degree=deg)
Q = fem.functionspace(domain, Qe)
dx = ufl.Measure("dx", domain=domain, metadata={"quadrature_degree": deg})


# --.. ..- .-.. .-.. ---  3. NN route: external-operator Newton  --.. ..- .-.. .-.. ---
T = fem.Function(V, name="T_nn")
T.interpolate(_linear_guess)

k = FEMExternalOperator(T, function_space=Q, external_function=k_external)
F = k * ufl.inner(ufl.grad(T), ufl.grad(v)) * dx
J = ufl.algorithms.expand_derivatives(ufl.derivative(F, T, dT))

F_replaced, F_ops = replace_external_operators(F)
J_replaced, J_ops = replace_external_operators(J)
F_form = fem.form(F_replaced)
J_form = fem.form(J_replaced)


def update_operators():
    ev = evaluate_operands(F_ops)
    evaluate_external_operators(F_ops, ev)
    evaluate_external_operators(J_ops, ev)


dT_sol = fem.Function(V)
ksp = PETSc.KSP().create(comm)
ksp.setType("preonly")
ksp.getPC().setType("lu")

print("\nNN-route Newton:")
for it in range(25):
    update_operators()
    Amat = assemble_matrix(J_form, bcs=bcs)
    Amat.assemble()
    bvec = assemble_vector(F_form)
    apply_lifting(bvec, [J_form], [bcs], [T.x.petsc_vec], -1.0)
    bvec.ghostUpdate(addv=PETSc.InsertMode.ADD, mode=PETSc.ScatterMode.REVERSE)
    set_bc(bvec, bcs, T.x.petsc_vec, -1.0)

    ksp.setOperators(Amat)
    ksp.solve(bvec, dT_sol.x.petsc_vec)
    dT_sol.x.scatter_forward()
    T.x.petsc_vec.axpy(-1.0, dT_sol.x.petsc_vec)
    T.x.scatter_forward()

    rnorm = bvec.norm()
    dnorm = dT_sol.x.petsc_vec.norm()
    print(f"  it {it:2d}: |residual| = {rnorm:.3e}   |correction| = {dnorm:.3e}")
    Amat.destroy()
    bvec.destroy()
    if rnorm < 1e-8:
        break


# --.. ..- .-.. .-.. ---  4. symbolic reference: dolfinx NewtonSolver  --.. ..- .-.. .-.. ---
T_sym = fem.Function(V, name="T_sym")
T_sym.interpolate(_linear_guess)
k_sym = A_COEF + B_COEF * T_sym
F_sym = k_sym * ufl.inner(ufl.grad(T_sym), ufl.grad(v)) * dx
problem = NonlinearProblem(
    F_sym, T_sym, bcs=bcs, petsc_options_prefix="zsym_",
    petsc_options={"snes_rtol": 1e-10, "snes_atol": 1e-12,
                   "ksp_type": "preonly", "pc_type": "lu"},
)
problem.solve()
snes = problem.solver
print(f"\nsymbolic SNES: converged_reason={snes.getConvergedReason()} "
      f"in {snes.getIterationNumber()} iterations")


# --.. ..- .-.. .-.. ---  5. compare against analytic + symbolic  --.. ..- .-.. .-.. ---
coords = V.tabulate_dof_coordinates()
xq = coords[:, 0]
T_nn_arr = T.x.array
T_sym_arr = T_sym.x.array
T_an_arr = T_analytic(xq)

err_nn_analytic = np.max(np.abs(T_nn_arr - T_an_arr))
err_nn_symbolic = np.max(np.abs(T_nn_arr - T_sym_arr))
err_sym_analytic = np.max(np.abs(T_sym_arr - T_an_arr))

print("\n--- comparison (max |dT| over all dofs, K) ---")
print(f"  NN vs analytic   : {err_nn_analytic:.4e}")
print(f"  NN vs symbolic   : {err_nn_symbolic:.4e}")
print(f"  symbolic vs anal : {err_sym_analytic:.4e}  (pure discretisation error)")


# --.. ..- .-.. .-.. ---  6. Taylor-remainder test on the residual (certify dk/dT)  --.. ..- .-.. .-.. ---
# At base state T_b, with perturbation direction p (zero on Dirichlet dofs):
#   R0 = |R(T_b + h p) - R(T_b)|                 -> O(h)   (slope ~1)
#   R1 = |R(T_b + h p) - R(T_b) - h J(T_b) p|    -> O(h^2) (slope ~2, certifies J=dR/dT)
rng = np.random.default_rng(0)
base = T.x.array.copy()
p = fem.Function(V)
p.x.array[:] = rng.standard_normal(base.shape)
for bc in bcs:
    p.x.array[bc._cpp_object.dof_indices()[0]] = 0.0   # homogeneous on Dirichlet dofs


def residual_at(state):
    T.x.array[:] = state
    T.x.scatter_forward()
    ev = evaluate_operands(F_ops)
    evaluate_external_operators(F_ops, ev)
    b = assemble_vector(F_form)
    b.ghostUpdate(addv=PETSc.InsertMode.ADD, mode=PETSc.ScatterMode.REVERSE)
    set_bc(b, bcs, T.x.petsc_vec, -1.0)
    out = b.array.copy()
    b.destroy()
    return out


# base residual + Jacobian action J(T_b) p
T.x.array[:] = base
T.x.scatter_forward()
update_operators()
Amat = assemble_matrix(J_form, bcs=bcs)
Amat.assemble()
R_base = residual_at(base)
Jp = Amat.createVecRight()
Amat.mult(p.x.petsc_vec, Jp)
Jp_arr = Jp.array.copy()

hs = np.array([1e-1, 1e-2, 1e-3, 1e-4, 1e-5])
R0, R1 = [], []
for h in hs:
    Rh = residual_at(base + h * p.x.array)
    R0.append(np.linalg.norm(Rh - R_base))
    R1.append(np.linalg.norm(Rh - R_base - h * Jp_arr))
R0, R1 = np.array(R0), np.array(R1)
slope0 = np.diff(np.log(R0)) / np.diff(np.log(hs))
slope1 = np.diff(np.log(R1)) / np.diff(np.log(hs))
T.x.array[:] = base  # restore
T.x.scatter_forward()

print("\n--- Taylor-remainder test (residual) ---")
print("   h         R0 (no grad)   R1 (with grad)")
for i, h in enumerate(hs):
    print(f"  {h:.0e}    {R0[i]:.3e}      {R1[i]:.3e}")
print(f"  mean slope R0 = {np.mean(slope0):.2f}  (expect ~1)")
print(f"  mean slope R1 = {np.mean(slope1):.2f}  (expect ~2  -> dk/dT certified)")


# --.. ..- .-.. .-.. ---  7. figure + verdict  --.. ..- .-.. .-.. ---
order = np.argsort(xq)
fig, ax = plt.subplots(figsize=(7, 4.5))
xx = np.linspace(0, W, 200)
ax.plot(xx, T_analytic(xx), "k-", lw=2, label="analytic")
ax.plot(xq[order], T_sym_arr[order], "C0--", lw=1.5, label="symbolic k(T) (NewtonSolver)")
ax.plot(xq[order], T_nn_arr[order], "C3o", ms=3, label="NN k(T) (external operator)")
ax.set_xlabel("x [m]")
ax.set_ylabel("T [K]")
ax.set_title("NN conductivity verification: steady slab, k(T)=a+bT")
ax.legend()
fig.tight_layout()
fig_path = os.path.join(HERE, "nn_conductivity_verify.png")
fig.savefig(fig_path, dpi=130)
print(f"\nfigure -> {fig_path}")

# verdict: NN must match symbolic to ~ the NN fit tolerance, and the tangent must be order 2
ok_phys = err_nn_symbolic < 1.0          # within 1 K of the symbolic discrete solution
ok_tang = np.mean(slope1) > 1.8
print("\n=== VERDICT ===")
print(f"  physics  (NN == symbolic, <1 K): {'PASS' if ok_phys else 'FAIL'}  ({err_nn_symbolic:.3e} K)")
print(f"  tangent  (Taylor slope > 1.8)  : {'PASS' if ok_tang else 'FAIL'}  ({np.mean(slope1):.2f})")
assert ok_phys and ok_tang, "verification FAILED"
print("  ALL CHECKS PASSED")
