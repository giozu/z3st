import os
import subprocess
import re
import numpy as np
from z3st.utils.utils_extract_vtu import extract_field

# Configuration
GEO_FILE = "mesh.geo"
VTU_FILE = os.path.join("output", "fields.vtu")
nx_values = [6, 11, 16, 21, 26, 31, 36, 41]
results = []

# Data
Ri, Ro, Lz = 2.0, 2.1, 20
k, E, nu, alpha = 48.1, 1.77e11, 0.3, 1.7e-5
Ti, To = 490.0, 500.0
q0, mu = 2.0e6, 24.0
Lx = Ro - Ri

def update_geo(nx):
    with open(GEO_FILE, "r") as f:
        content = f.read()
    content = re.sub(r"nx\s*=\s*\d+;", f"nx = {nx};", content)
    with open(GEO_FILE, "w") as f:
        f.write(content)

def analytic_T(x):
    """Analytical solution for Dirichlet-Dirichlet cylindrical shell (slab approx)."""
    x_rel = x - Ri
    term1 = Ti + (To - Ti) * (x_rel / Lx)
    term2 = (q0 / (mu**2 * k)) * ((x_rel / Lx) * (np.exp(-mu * Lx) - 1) - (np.exp(-mu * x_rel) - 1))
    return term1 + term2

# Sensitivity loop
for nx in nx_values:
    print(f"Testing mesh: {nx-1} radial elements...")
    
    # Mesh & Solve
    update_geo(nx)
    subprocess.run(["gmsh", GEO_FILE, "-2", "-o", "mesh.msh", "-format", "msh2"], check=True, stdout=subprocess.DEVNULL)
    subprocess.run(["python3", "-m", "z3st"], check=True, stdout=subprocess.DEVNULL)
    
    # Extraction
    x_all, y_all, z_all, T_all = extract_field(VTU_FILE, field_name="Temperature")
    
    # L2 Error Calculation
    T_ref_nodes = analytic_T(x_all)
    diff_sq = (T_all - T_ref_nodes)**2
    l2_err = np.sqrt(np.sum(diff_sq) / np.sum(T_ref_nodes**2))
    
    print(f"   -> Nodes: {len(x_all)} | Global L2 Error: {l2_err:.4e}")
    
    # Append [mesh_size_h, error]
    results.append([1.0/(nx-1), l2_err])

# Save
np.savetxt("convergence_data.txt", np.array(results))
print("\n[INFO] convergence_data.txt generated successfully.")