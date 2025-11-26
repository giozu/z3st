# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.1.0 (2025)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

import matplotlib.pyplot as plt
import numpy as np
import dolfinx
import ufl
import scipy.sparse
import scipy.sparse.linalg as sla
from petsc4py import PETSc


# --. Assembly + analysis utilities --..
def assemble_and_analyze_matrix(a, L=None, bcs=None, title="Matrix", max_dense_size=5000, tol=1e-5):
    """
    Assemble PETSc matrix from variational form, apply boundary conditions,
    analyze sparsity and estimate condition number (sparse-safe).
    """
    if L is not None:
        F = a - L
        a, _ = ufl.system(F)

    a_compiled = dolfinx.fem.form(a)
    A = dolfinx.fem.petsc.assemble_matrix(a_compiled, bcs=bcs or [])
    A.assemble()

    # --- Basic info
    n = A.getSize()[0]
    nnz = A.getInfo()["nz_used"]
    density = nnz / (n**2)
    print(f"Matrix size: {n}x{n} ({n} DOFs)")
    print(f"Non-zero entries: {nnz:,}")
    print(f"Sparsity density: {density:.2e}")

    # --- Symmetry check
    is_sym = is_symmetric(A, tol)
    print(f"Symmetric (tolerance {tol:.0e}): {'Yes' if is_sym else 'No'}")

    # --- Condition number
    cond_est = None
    if n <= max_dense_size:
        cond_est = dense_condition(A)
        print(f"Condition number (dense): {cond_est:.3e}")
    else:
        print(f"[INFO] Skipping dense cond(A): {n} DOFs > {max_dense_size}")
        cond_est = estimate_condition_number_sparse(A)
        print(f"Estimated cond(A) (iterative): {cond_est:.3e}")
        if cond_est and cond_est > 1e8:
            print(f"[WARN] High condition number ({cond_est:.2e}) → possible under-constraining or soft modes.")

    print("-------------------------------------------------")
    return A


def is_symmetric(A, tol=1e-8):
    """Check approximate matrix symmetry."""
    A_t = A.copy()
    A_t.transpose()
    diff = A_t
    diff.axpy(-1.0, A)
    norm_diff = diff.norm(PETSc.NormType.FROBENIUS)
    norm_A = A.norm(PETSc.NormType.FROBENIUS)
    return norm_diff / (norm_A + 1e-16) < tol


def dense_condition(A):
    """Compute condition number using dense conversion (small systems only)."""
    A_dense = A.convert("dense")
    arr = A_dense.getDenseArray()
    s = np.linalg.svd(arr, compute_uv=False)
    return s[0] / s[-1]


def estimate_condition_number_sparse(A, max_iter=200):
    """
    Estimate condition number for large sparse matrices using PETSc CG eigenvalues.
    cond(A) ≈ λ_max / λ_min (estimated iteratively)
    """
    ksp = PETSc.KSP().create()
    ksp.setOperators(A)
    ksp.setType("cg")
    ksp.setTolerances(rtol=1e-6, max_it=max_iter)
    ksp.setComputeEigenvalues(True)

    # PETSc richiede vettori validi per almeno una solve()
    x = A.createVecRight()
    b = A.createVecLeft()
    b.setRandom()
    ksp.solve(b, x)

    try:
        eigs = np.array(ksp.computeEigenvalues())
        eigs = eigs[np.isfinite(eigs) & (eigs > 0)]
        if len(eigs) > 0:
            return np.max(eigs) / np.min(eigs)
        else:
            print("[WARN] No valid eigenvalues computed.")
            return np.nan
    except Exception as e:
        print(f"[WARN] Eigenvalue estimation failed: {e}")
        return np.nan


# --. Matrix-level analysis (standalone) --..
def analyze_matrix(A, title="Matrix", max_dense_size=5000, plot=False, tol=1e-5):
    """
    Analyze PETSc matrix properties: size, sparsity, symmetry, condition number.
    """
    size = A.getSize()[0]
    print(f"Matrix size: {size}x{size} ({size} DOFs)")

    if size <= max_dense_size:
        print("Matrix is small enough for dense analysis.")
        try:
            ai, aj, av = A.getValuesCSR()
            A_sp = scipy.sparse.csr_matrix((av, aj, ai), shape=A.getSize())

            if plot:
                plt.figure(figsize=(6, 6))
                plt.spy(A_sp, markersize=1, aspect='equal')
                plt.title(f"Sparsity Pattern: {title}")
                plt.xlabel("Column index")
                plt.ylabel("Row index")
                plt.grid(True)
                plt.savefig('output/sparsity_dense_matrix.png')
                plt.close()

            print("Converting to dense to calculate condition number...")
            A_np = A_sp.toarray()
            cond_number = np.linalg.cond(A_np)
            print(f"Condition number (dense): {cond_number:.2e}")

        except Exception as e:
            print(f"Could not perform dense analysis. Reason: {e}")
    else:
        print(f"Matrix is too large for dense analysis (>{max_dense_size} DOFs).")
        try:
            _, _, av = A.getValuesCSR()
            nnz = len(av)
            density = nnz / (size * size)
            print(f"Non-zero entries: {nnz:,}")
            print(f"Sparsity density: {density:.2e}")
        except Exception as e:
            print(f"Could not retrieve number of non-zero entries. Reason: {e}")

    is_sym = "Yes" if A.isSymmetric(tol=tol) else "No"
    print(f"Symmetric (tolerance 1e-8): {is_sym}")
    print("-" * (len(title) + 22) + "\n")


# --. Main mechanical constraint diagnostic --..
def check_mechanical_constraints(problem, max_dense_size=5000, eig_check=False, plot=False, tol=1e-5):
    """
    Check whether the mechanical problem is correctly constrained
    (under-constrained, well-constrained, or over-constrained).
    """
    print("\n[DIAGNOSTIC] Checking mechanical constraints...")

    u = ufl.TrialFunction(problem.V_m)
    v = ufl.TestFunction(problem.V_m)

    # Build stiffness form over all subdomains
    a_form = 0
    for label, mat in problem.materials.items():
        tag = problem.label_map[label]
        dx = ufl.Measure("dx", domain=problem.mesh, subdomain_data=problem.volume_tags, subdomain_id=tag)
        a_form += ufl.inner(problem.sigma_mech(u, mat), problem.epsilon(v)) * dx

    # Collect Dirichlet BCs
    bcs = []
    if hasattr(problem, "dirichlet_mechanical"):
        for bc_list in problem.dirichlet_mechanical.values():
            bcs.extend(bc_list)

    A = assemble_and_analyze_matrix(a_form, bcs=bcs, title="Mechanical stiffness matrix", tol=tol)

    # Convert PETSc to SciPy CSR for further diagnostics
    ai, aj, av = A.getValuesCSR()
    K = scipy.sparse.csr_matrix((av, aj, ai), shape=A.getSize())

    # Sparsity plot
    if plot:
        plt.figure(figsize=(6, 6))
        plt.spy(K, markersize=1, aspect='equal')
        plt.title("Mechanical stiffness matrix sparsity")
        plt.xlabel("Column index")
        plt.ylabel("Row index")
        plt.grid(True)
        plt.savefig('output/sparsity_matrix.png')
        plt.close()

    # Condition number (only if small)
    cond_number = None
    if K.shape[0] <= max_dense_size:
        K_dense = K.toarray()
        cond_number = np.linalg.cond(K_dense)
        print(f"Condition number: {cond_number:.2e}")
        diff_norm = np.linalg.norm(K_dense - K_dense.T)
        print("‖K - Kᵀ‖ =", diff_norm)
    else:
        print(f"[INFO] Skipping dense cond(A): {K.shape[0]} DOFs > {max_dense_size}")

    # Eigenvalues (optional)
    eigvals = None
    if eig_check:
        try:
            print("[INFO] Estimating smallest eigenvalues (rigid-body check)...")
            # k = min(6, K.shape[0]-1)
            k = 2

            eigvals = sla.eigsh(K, k=k, which='SM', return_eigenvectors=False)
            print("Eigenvalues (smallest):", np.round(eigvals, 6))
        except Exception as e:
            print(f"[WARN] Eigenvalue analysis failed: {e}")

    # Symmetry
    is_sym = A.isSymmetric(tol=tol)
    print(f"Symmetric: {'YES' if is_sym else 'NO'}")

    # Classification heuristic
    status = "UNKNOWN"
    if cond_number is not None:
        if cond_number > 1e12 or (eigvals is not None and np.min(eigvals) < 1e-8):
            status = "UNDER-CONSTRAINED"
        elif cond_number < 1e3 and not is_sym:
            status = "OVER-CONSTRAINED"
        else:
            status = "WELL-CONSTRAINED"

    print(f"\n[RESULT] Mechanical system appears: {status}")

    return {
        "size": K.shape[0],
        "cond": cond_number,
        "symmetry": is_sym,
        "eigvals": eigvals,
        "status": status,
    }


# --. Fixed DOFs visualization --..
def check_fixed_dofs(problem, fname):
    """
    Visualize all nodes constrained by mechanical Dirichlet BCs.
    Correct mapping from vector-valued DOFs to mesh node coordinates.
    """
    import pyvista as pv

    print("\n[DIAGNOSTIC] Visualizing fixed mechanical DOFs...")

    V = problem.V_m
    mesh = problem.mesh
    gdim = mesh.geometry.dim

    fixed_dofs = []

    # Collect all Dirichlet BC indices
    for bc_list in problem.dirichlet_mechanical.values():
        for bc in bc_list:
            try:
                dofs = bc.dof_indices()
                if isinstance(dofs, (list, tuple)):
                    for d in dofs:
                        fixed_dofs.extend(np.atleast_1d(d).tolist())
                else:
                    fixed_dofs.extend(np.atleast_1d(dofs).tolist())
            except Exception as e:
                print(f"  [WARN] Cannot extract dof indices from {bc}: {e}")
                continue

    if len(fixed_dofs) == 0:
        print("  [INFO] No Dirichlet mechanical BCs found.")
        return

    facets = problem.ft.find(problem.label_map[fname])
    print("Facets with Dirichlet BCs:", len(facets))

    fixed_dofs = np.unique(np.array(fixed_dofs, dtype=int).ravel())

    # Map vector DOFs to corresponding nodes
    n_nodes_total = mesh.geometry.x.shape[0]
    node_indices = np.floor_divide(fixed_dofs, gdim)
    node_indices = np.unique(node_indices)
    coords = mesh.geometry.x[node_indices]

    # Visualize
    plotter = pv.Plotter()
    plotter.add_mesh(pv.PolyData(mesh.geometry.x),
                     color="lightgrey", point_size=3, render_points_as_spheres=True)
    plotter.add_points(coords, color="red", point_size=10, render_points_as_spheres=True)
    plotter.add_legend([("Fixed nodes", "red")])
    plotter.show()

    print(f"  [INFO] Displayed {len(coords)} fixed nodes (out of {n_nodes_total} total).")


# --. Debug BC summary --..
def debug_dirichlet_mechanical(problem):
    """
    Print a structured summary of all mechanical Dirichlet BCs.
    Robust to scalar/vector spaces.
    """
    print("\n[DEBUG] Mechanical Dirichlet boundary conditions summary")
    print("──────────────────────────────────────────────────────────")

    # Summary per material
    for mat_name, bc_list in problem.dirichlet_mechanical.items():
        print(f"  - Material '{mat_name}' : {len(bc_list)} Dirichlet BC(s)")
    if not problem.dirichlet_mechanical:
        print("  [INFO] No mechanical Dirichlet BCs found.")
        return

    print("\n[DETAILS PER BC OBJECT]")
    print("────────────────────────")

    all_dofs = []
    for mat_name, bc_list in problem.dirichlet_mechanical.items():
        for i, bc in enumerate(bc_list):
            try:
                dofs_raw = bc.dof_indices()
                flat = []
                if isinstance(dofs_raw, (list, tuple)):
                    for p in dofs_raw:
                        flat.extend(np.atleast_1d(p).ravel().tolist())
                else:
                    flat.extend(np.atleast_1d(dofs_raw).ravel().tolist())
                all_dofs.extend(flat)
                print(f"  : {mat_name}[{i}]  |  {len(flat):5d} DOFs  "
                      f"(first 10): {flat[:10]}")
            except Exception as e:
                print(f"  [WARN] Could not read dof_indices: {e}")

    all_dofs = np.unique(np.array(all_dofs, dtype=int))

    try:
        V = problem.V_m
        gdim = problem.mesh.geometry.dim
        x_dofs = V.tabulate_dof_coordinates().reshape((-1, gdim))

        print(f"[DEBUG] V_m.dim = {V.dofmap.index_map.size_global}")
        print(f"[DEBUG] tabulated coordinates shape = {x_dofs.shape}")
        print(f"[DEBUG] max DOF requested = {all_dofs.max()}")

        valid_mask = all_dofs < x_dofs.shape[0]
        if not np.all(valid_mask):
            offenders = all_dofs[~valid_mask]
            print(f"[WARN] {offenders.size} DOFs out of bounds (ignored): {offenders[:10]}")
            all_dofs = all_dofs[valid_mask]

        coords = x_dofs[all_dofs]
    except Exception as e:
        print(f"  [FATAL] Failed to get DOF coordinates: {e}")
        return

    print("\n[COORDINATE RANGE OF CONSTRAINED DOFs]")
    print("────────────────────────────────────────")
    for i, name in enumerate(["x", "y", "z"][:gdim]):
        print(f"  {name}-range: {coords[:, i].min():.4f} → {coords[:, i].max():.4f}")

    try:
        print(f"  [CHECK] Unique Dirichlet DOFs : {all_dofs.size}")
        tol = 1e-12
        mask = np.abs(coords[:, 0] - 0.0) > tol
        if np.any(mask):
            offenders = coords[mask]
            print(f"  [WARN] {offenders.shape[0]} constrained DOFs off-plane x=0:")
            print(np.round(offenders, 8))
        else:
            print(f"  [CHECK] All constrained DOFs lie on x=0 (±{tol})")
    except Exception as e:
        print(f"  [WARN] DOF uniqueness/off-plane check failed: {e}")

    print("──────────────────────────────────────────────────────────\n")
