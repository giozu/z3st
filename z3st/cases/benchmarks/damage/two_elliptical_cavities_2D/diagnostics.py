import os
import dolfinx
import ufl
from mpi4py import MPI

CASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_FILE = os.path.join(CASE_DIR, "output", "stress_strain.txt")

def _init_diagnostics(problem):
    ymax_tag = problem.label_map.get("ymax")
    
    ds_ymax = ufl.Measure(
        "ds",
        domain=problem.mesh,
        subdomain_data=problem.facet_tags,
        subdomain_id=ymax_tag,
    )

    mat_name = next(iter(problem.materials))
    
    # On extrait la contrainte yy (indice [1, 1]) et le déplacement u_y (indice [1])
    sigma_yy = problem.stress[mat_name][1, 1]
    u_y = problem.u[1]

    per_step.force_form = dolfinx.fem.form(sigma_yy * ds_ymax)
    per_step.disp_form = dolfinx.fem.form(u_y * ds_ymax)
    
    # Forme pour calculer l'aire/longueur du bord
    unity = dolfinx.fem.Constant(problem.mesh, dolfinx.default_scalar_type(1.0))
    per_step.area_form = dolfinx.fem.form(unity * ds_ymax)
    
    per_step.initialized = True

    if problem.mesh.comm.rank == 0:
        os.makedirs(os.path.dirname(OUTPUT_FILE), exist_ok=True)
        with open(OUTPUT_FILE, "w") as f:
            f.write("Step\tStrain_yy\tStress_yy_MPa\n")

def per_step(problem, step, t):
    if not hasattr(per_step, "initialized"):
        _init_diagnostics(problem)

    # Assemblage et intégration exacte sur le bord
    F_yy = problem.mesh.comm.allreduce(dolfinx.fem.assemble_scalar(per_step.force_form), op=MPI.SUM)
    U_y = problem.mesh.comm.allreduce(dolfinx.fem.assemble_scalar(per_step.disp_form), op=MPI.SUM)
    A_ymax = problem.mesh.comm.allreduce(dolfinx.fem.assemble_scalar(per_step.area_form), op=MPI.SUM)

    if A_ymax > 0:
        avg_sigma_yy = F_yy / A_ymax
        avg_u_y = U_y / A_ymax
    else:
        avg_sigma_yy, avg_u_y = 0.0, 0.0

    Ly = getattr(problem.mgr, "Ly", 0.560499)
    
    strain_yy = avg_u_y / Ly if Ly > 0 else 0.0
    stress_yy_MPa = avg_sigma_yy * 1e-6

    if problem.mesh.comm.rank == 0:
        with open(OUTPUT_FILE, "a") as f:
            f.write(f"{step}\t{strain_yy:.6e}\t{stress_yy_MPa:.6e}\n")