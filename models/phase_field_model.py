# models/phase_field_model.py
import numpy as np
import ufl, dolfinx
from petsc4py import PETSc
import yaml

class PhaseFieldModel:
    """
    Phase-field fracture (AT-2) with traction/compression split (Amor) and
    irreversibility enforced via a history field H(x) = max_tau Psi_e^+(tau).

    Strong form implemented (with s ≡ d, damage):
        2*ell * Δd + (1 - d)/(2*ell) = (∂g/∂d)(d) * Psi_e^+ / Gc,
    with g_d(d) = (1 - d)^2  →  ∂g/∂d = -2(1 - d).

    Equivalent weak form used here (standard AT-2 scaling):
        ∫ [ (Gc/ell) d w + Gc*ell ∇d·∇w + 2 H d w ] dx  =  ∫ 2 H w dx

    PUBLIC API used by the solver:
      - setup_phasefield_space()
      - update_history(u)           # updates H from current u
      - build_damage_forms()        # builds a_d_form, L_d_form
      - assemble_solve_damage()     # solves for d with CG+AMG
      - sigma_damage_linear(u_trial, d, material) # degraded, linear in u_trial
    """

    # -------- init & parameters --------
    def __init__(self):
        print("__PhaseFieldModel initializer__")

        input_file = "input.yaml"
        with open(input_file, "r") as f:
            data = yaml.safe_load(f)

        pf = data.get("fracture", {})
        if not pf:
            raise ValueError("[PhaseFieldModel] 'fracture' missing in input.yaml.")

        self.model = pf.get("model", "AT2")
        self.Gc = pf.get("Gc", 1000.0)
        self.ell = pf.get("ell", 0.01)
        self.kappa = pf.get("kappa", 1e-8)
        self.split = pf.get("split", "none")
        self.irreversibility = pf.get("irreversibility", "none")
        self.sigma_c = pf.get("sigma_c", 10.0)

        print("\n[PhaseFieldModel] parameters loaded from input.yaml:")
        print(f"  model            : {self.model}")
        print(f"  Gc               : {self.Gc}")
        print(f"  ell              : {self.ell}")
        print(f"  kappa            : {self.kappa}")
        print(f"  split            : {self.split}")
        print(f"  irreversibility  : {self.irreversibility}\n")
        print(f"  sigma_c          : {self.sigma_c}")

        self.H_threshold = 0.5 * self.sigma_c**2 * self.ell/ self.Gc
        print(f"  Damage threshold energy H_threshold = {self.H_threshold:.2e}")


        # FE containers
        self.V_d = None
        self.d   = None           # damage field
        self.d_n = None           # previous damage (optional)
        self.H   = None           # history field (Psi_e^+ max-in-time)
        self.d_trial = None       # TrialFunction(d)
        self.w_test  = None       # TestFunction(w)

        self.dirichlet_damage = []   # optional Dirichlet BCs on d
        self.a_d_form = None
        self.L_d_form = None

    # -------- helpers --------
    def _I(self): return ufl.Identity(self.gdim)

    @staticmethod
    def _pos(x):  # <x>_+
        return (x + ufl.sqrt(x*x)) / 2.0

    def _dev(self, eps):
        return eps - (ufl.tr(eps)/3.0) * self._I()

    def g_d(self, d):
        """Degradation function with small residual stiffness."""
        return (1.0 - d)**2 + self.kappa

    # -------- function spaces & fields --------
    def setup_phasefield_space(self):
        """Create P1 scalar space for d and H and initialize to zero."""
        self.V_d = dolfinx.fem.functionspace(self.mesh, ("Lagrange", 1))

        self.d   = dolfinx.fem.Function(self.V_d, name="Damage")
        self.d.x.array[:] = 0.0

        self.d_n = dolfinx.fem.Function(self.V_d, name="Damage_n")
        self.d_n.x.array[:] = 0.0

        self.H   = dolfinx.fem.Function(self.V_d, name="History")
        self.H.x.array[:] = 0.0

        self.d_trial = ufl.TrialFunction(self.V_d)
        self.w_test  = ufl.TestFunction(self.V_d)

    # -------- elastic split (Amor) --------
    def psi_plus(self, u, material):
        """
        Psi_e^+(u) = 0.5*K*<tr ε>_+^2 + G * ||dev ε||^2  (Amor dev/vol split)
        """
        eps = self.epsilon(u)            # from MechanicalModel mixin
        tr  = ufl.tr(eps)
        dev = self._dev(eps)

        lam = material["lmbda"]
        G   = material["G"]
        K   = lam + 2.0*G/3.0

        tr_p = self._pos(tr)
        return 0.5*K*(tr_p*tr_p) + G*ufl.inner(dev, dev)

    # Indicator chi_vol to freeze tension/compression split for linearized mech step
    def update_split_indicator(self, u_ref):
        tr_prev = ufl.tr(self.epsilon(u_ref))
        expr = dolfinx.fem.Expression(tr_prev, self.V_d.element.interpolation_points())
        tmp  = dolfinx.fem.Function(self.V_d)
        tmp.interpolate(expr)

        if not hasattr(self, "chi_vol"):
            self.chi_vol = dolfinx.fem.Function(self.V_d, name="ChiVol")
        self.chi_vol.x.array[:] = (tmp.x.array >= 0.0).astype(self.chi_vol.x.array.dtype)

    def sigma_damage_linear(self, u_trial, d, material):
        eps = self.epsilon(u_trial)
        tr  = ufl.tr(eps)
        dev = self._dev(eps)

        lam = material["lmbda"]
        G   = material["G"]
        K   = lam + 2.0*G/3.0
        chi = self.chi_vol  # Function(V_d), frozen

        gdeg = (1.0 - d)**2 + self.kappa  # <-- qui
        return gdeg * (K * chi * tr * self._I() + 2.0 * G * dev) \
            + (1.0 - chi) * K * tr * self._I()

    def psi_split(self, u, material):
        """
        Compute the positive and negative parts of the elastic energy.

        Parameters:
            u : displacement field
            material : material object or dictionary

        Returns:
            psi_plus : positive part (drives damage)
            psi_minus : negative part (ignored in damage evolution)
        """
        eps   = ufl.sym(ufl.grad(u))
        lam   = material["lmbda"]
        G     = material["G"]

        tr_eps  = ufl.tr(eps)
        dev_eps = eps - tr_eps / 3.0 * self._I()

        tr_p = ufl.max_value(tr_eps, 0.0)
        tr_m = ufl.min_value(tr_eps, 0.0)

        psi_pos = 0.5 * (lam + 2.0*G/3.0) * (tr_p*tr_p) + G * ufl.inner(dev_eps, dev_eps)
        psi_neg = 0.5 * (lam + 2.0*G/3.0) * (tr_m*tr_m)

        return psi_pos, psi_neg


    def update_history(self, u):
        """
        H(x) := max( H(x), Psi_e^+(u)(x) ) su tutti i DOF P1,
        ma solo nelle celle del materiale corrente (mascherate).
        """
        V0 = dolfinx.fem.functionspace(self.mesh, ("DG", 0))      # cell-wise
        V1 = self.V_d                                            # P1

        H_new = dolfinx.fem.Function(V1)
        H_new.x.array[:] = self.H.x.array

        for name, mat in self.materials.items():
            tag = self.label_map[name]
            cells = self.volume_tags.find(tag)

            # 1) indicatore cell-wise
            chi0 = dolfinx.fem.Function(V0)
            chi0.x.array[:] = 0.0
            chi0.x.array[cells] = 1.0

            # 2) portiamo la maschera su P1 (interpolazione ai punti di interpolazione di V1)
            chi1 = dolfinx.fem.Function(V1)
            chi1.interpolate(dolfinx.fem.Expression(chi0, V1.element.interpolation_points()))

            # 3) ψ⁺ su P1
            psi_p = self.psi_plus(u, mat)  # <-- usa proprio psi_plus
            psi1  = dolfinx.fem.Function(V1)
            psi1.interpolate(dolfinx.fem.Expression(psi_p, V1.element.interpolation_points()))

            # 4) update con maschera
            mask = (chi1.x.array > 0.5)
            H_arr   = H_new.x.array
            psi_arr = psi1.x.array
            np.maximum(H_arr[mask], psi_arr[mask], out=H_arr[mask])


        self.H.x.array[:] = H_new.x.array

    # -------- damage forms --------
    def build_damage_forms(self):
        """
        Build bilinear and linear forms of the AT-2 damage subproblem:
          a(d,w) = ∫ [ (Gc/ell) d w + Gc*ell ∇d·∇w + 2 H d w ] dx
          L(w)   = ∫ 2 H w dx
        """
        dx = ufl.Measure("dx", domain=self.mesh)

        a_bulk = (self.Gc/self.ell) * self.d_trial * self.w_test \
                 + (self.Gc*self.ell) * ufl.dot(ufl.grad(self.d_trial), ufl.grad(self.w_test))

        a_hist = 2.0 * self.H * self.d_trial * self.w_test

        self.a_d_form = dolfinx.fem.form((a_bulk + a_hist) * dx)

        H_eff = ufl.max_value(self.H - self.H_threshold, 0.0)
        self.L_d_form = dolfinx.fem.form((2.0 * H_eff * self.w_test) * dx)

        print("[PF] damage forms built.")

    # -------- linear solve for d --------
    def assemble_solve_damage(self):
        """
        Solve the linear system for d with (optional) Dirichlet BCs in self.dirichlet_damage.
        """
        from dolfinx.fem import petsc as fem_petsc

        A = fem_petsc.assemble_matrix(self.a_d_form, bcs=self.dirichlet_damage)
        A.assemble()
        b = fem_petsc.assemble_vector(self.L_d_form)
        fem_petsc.apply_lifting(b, [self.a_d_form], [self.dirichlet_damage])
        b.ghostUpdate(addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE)
        fem_petsc.set_bc(b, self.dirichlet_damage)

        # KSP: symmetric SPD → CG with AMG preconditioner
        ksp = PETSc.KSP().create(self.mesh.comm)
        ksp.setOperators(A)
        ksp.setType("cg")
        pc = ksp.getPC()
        pc.setType("hypre")  # or "gamg"
        ksp.setFromOptions()

        d_sol = dolfinx.fem.Function(self.V_d, name="Damage_tmp")
        ksp.solve(b, d_sol.x.petsc_vec)
        d_sol.x.scatter_forward()
        np.clip(d_sol.x.array, 0.0, 1.0, out=d_sol.x.array)

        # irreversibility: d >= d_n
        if self.d_n is not None:
            np.maximum(d_sol.x.array, self.d_n.x.array, out=d_sol.x.array)

        self.d.x.array[:] = d_sol.x.array
        if self.d_n is not None:
            self.d_n.x.array[:] = self.d.x.array