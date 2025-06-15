import matplotlib.pyplot as plt
import numpy as np
import dolfinx
import ufl
import scipy.sparse
# It's good practice to have all imports at the top

def assemble_and_analyze_matrix(a, L=None, bcs=None, title="Matrix", max_dense_size=5000):
    """
    Assemble PETSc matrix from variational form, apply boundary conditions, 
    plot sparsity, and estimate condition number if feasible.

    Parameters:
    - a: ufl.Form
        Bilinear form (LHS of the variational problem).
    - L: ufl.Form, optional
        Linear form (RHS of the variational problem). If provided, it is used
        with 'a' to define the full residual for non-linear problems.
    - bcs: list of dolfinx.fem.DirichletBC, optional
        List of Dirichlet boundary conditions to apply to the matrix.
    - title: str, optional
        A label for plots and printed output.
    - max_dense_size: int, optional
        The maximum number of degrees of freedom for which a dense matrix
        analysis (condition number, spy plot) will be attempted.
    """
    # If a linear form L is provided, extract the bilinear form (Jacobian)
    # from the full residual F = a - L. This is useful for non-linear problems.
    if L is not None:
        F = a - L
        a, _ = ufl.system(F)

    # Compile and assemble the matrix A from the bilinear form a
    a_compiled = dolfinx.fem.form(a)
    A = dolfinx.fem.petsc.assemble_matrix(a_compiled, bcs=bcs or [])
    A.assemble()

    # Pass the assembled matrix to the analysis function
    analyze_matrix(A, title, max_dense_size)
    return A

def analyze_matrix(A, title="Matrix", max_dense_size=5000, plot = False):
    """
    Analyzes a PETSc matrix for its size, sparsity, symmetry, and conditioning.

    This function prints key properties of the matrix. For small matrices
    (size <= max_dense_size), it generates a sparsity plot and computes the
    condition number. For larger matrices, it reports the number of non-zero
    entries and sparsity density to avoid excessive memory usage.

    Parameters:
    - A: dolfinx.fem.petsc.Mat
        The assembled PETSc matrix to be analyzed.
    - title: str, optional
        A label for plots and printed output.
    - max_dense_size: int, optional
        The threshold for performing dense matrix analysis.
    """
    size = A.getSize()[0]
    print(f"Matrix size: {size}x{size} ({size} DOFs)")

    # Perform dense analysis only for reasonably small matrices
    if size <= max_dense_size:
        print("Matrix is small enough for dense analysis.")
        try:
            # Convert PETSc matrix to a SciPy CSR matrix for plotting and analysis
            ai, aj, av = A.getValuesCSR()
            A_sp = scipy.sparse.csr_matrix((av, aj, ai), shape=A.getSize())

            # 1. Sparsity pattern plot
            if plot is True:
                plt.figure(figsize=(6, 6))
                plt.spy(A_sp, markersize=1, aspect='equal')
                plt.title(f"Sparsity Pattern: {title}")
                plt.xlabel("Column index")
                plt.ylabel("Row index")
                plt.grid(True)
                plt.show()

            # 2. Condition number calculation
            print("Converting to dense to calculate condition number...")
            A_np = A_sp.toarray()
            cond_number = np.linalg.cond(A_np)
            print(f"Condition number (dense): {cond_number:.2e}")

        except Exception as e:
            print(f"Could not perform dense analysis. Reason: {e}")
    else:
        # For large matrices, avoid creating a dense version
        print(f"Matrix is too large for dense analysis (>{max_dense_size} DOFs).")
        try:
            # FIX: A more robust way to get the number of non-zero entries is to
            # retrieve the CSR data and check the length of the values array.
            # This is more reliable than A.getInfo().
            _, _, av = A.getValuesCSR()
            nnz = len(av)
            density = nnz / (size * size)
            print(f"Non-zero entries: {nnz:,}") # Use comma for thousands separator
            print(f"Sparsity density: {density:.2e}")
        except Exception as e:
            # Fallback if CSR data isn't available for some reason
            print(f"Could not retrieve number of non-zero entries. Reason: {e}")

    # 3. Symmetry check (works for both dense and sparse)
    is_sym = "Yes" if A.isSymmetric(tol=1e-8) else "No"
    print(f"Symmetric (tolerance 1e-8): {is_sym}")
    print("-" * (len(title) + 22) + "\n")