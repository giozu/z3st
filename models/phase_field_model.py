# models/phase_field_model.py
import numpy as np
import ufl, dolfinx
from petsc4py import PETSc

class PhaseFieldModel:
    """
    Frattura a campo di fase (AT-2) con split tipo Amor e irreversibilità via history-field.
    Si integra con ZEST come mixin (vedi spine.py / solver.py).
    """
    def __init__(self):
        print("__PhaseFieldModel initializer__")
        # Parametri da config (impostati in Config.load_input)
        pf = self.fracture_params = getattr(self, "fracture_params", {
            "Gc": 1000,     # energia di frattura
            "ell": 0.01,   # length-scale
            "kappa": 1e-8, # residual stiffness
            "model": "AT2",
            "split": "amor_dev_vol",
            "irreversibility": "history"
        })
        self.Gc = float(pf["Gc"])
        self.ell = float(pf["ell"])
        self.kappa = float(pf["kappa"])
        self.cw = 0.5  # AT-2: w(d)=d^2 -> c_w = 1/2  (Gerasimov & De Lorenzis) 

        self.V_d = None
        self.d = None       # funzione di danno
        self.d_n = None     # danno al passo precedente (se ti serve)
        self.H = None       # history field (max nel tempo di psi_plus)
        self.eta = None     # test
        self.wd = None      # trial

        # Contenitori forme d
        self.a_damage = {}
        self.L_damage = {}
        self.dirichlet_damage = []

    # ---------- helper ----------
    def _I(self): return ufl.Identity(self.tdim)

    @staticmethod
    def _pos(x):  # <x>_+ 
        return (x + ufl.sqrt(x*x)) / 2.0

    @staticmethod
    def _neg(x):  # <x>_- 
        return (x - ufl.sqrt(x*x)) / 2.0

    def _dev(self, eps):
        return eps - (ufl.tr(eps)/3.0) * self._I()

    def g(self, d):
        # funzione di degrado (residuo kappa per evitare singolarità)
        return (1.0 - d)**2 + self.kappa

    # ---------- spazi & campi ----------
    def setup_phasefield_space(self):
        # spazio scalare P1 per d e H
        self.V_d = dolfinx.fem.functionspace(self.mesh, ("Lagrange", 1))
        self.d   = dolfinx.fem.Function(self.V_d, name="Damage")
        self.d.x.array[:] = 0.0
        self.d_n = dolfinx.fem.Function(self.V_d, name="Damage_n")
        self.d_n.x.array[:] = 0.0

        self.H   = dolfinx.fem.Function(self.V_d, name="History")
        self.H.x.array[:] = 0.0

        self.wd  = ufl.TrialFunction(self.V_d)
        self.eta = ufl.TestFunction(self.V_d)

    # ---------- split energetico (Amor) ----------
    def psi_split(self, u, material):
        """
        ψ⁺ = 0.5*K <tr ε>_+^2 + G ||dev ε||^2
        ψ⁻ = 0.5*K <tr ε>_-^2
        """
        eps = self.epsilon(u)  # da MechanicalModel
        tr  = ufl.tr(eps)
        dev = self._dev(eps)

        lmbda = material["lmbda"]
        G     = material["G"]
        K     = lmbda + 2.0*G/3.0

        tr_p = self._pos(tr)
        tr_m = self._neg(tr)

        psi_plus  = 0.5*K*(tr_p*tr_p) + G*ufl.inner(dev, dev)
        psi_minus = 0.5*K*(tr_m*tr_m)
        return psi_plus, psi_minus

    def g_d(self, d):
        """
        Degradation function g(d) for the phase-field fracture model.
        For AT2: g(d) = (1 - d)^2
        """
        return (1.0 - d)**2

    # 2) Indicatore di trazione "congelato" dal campo u_ref (qui: u_n)
    def update_split_indicator(self, u_ref):
        """
        Costruisce χ_vol(x) = 1 se tr(eps(u_ref)) >= 0, altrimenti 0.
        Salvata in self.chi_vol (Function su V_d).
        """
        import numpy as np
        # tr(eps(u_ref)) su V_d
        tr_prev = ufl.tr(self.epsilon(u_ref))
        expr = dolfinx.fem.Expression(tr_prev, self.V_d.element.interpolation_points())
        tmp  = dolfinx.fem.Function(self.V_d)
        tmp.interpolate(expr)

        # χ_vol = Heaviside(tr_prev) (0/1 come array numpy)
        if not hasattr(self, "chi_vol"):
            self.chi_vol = dolfinx.fem.Function(self.V_d, name="ChiVol")
        self.chi_vol.x.array[:] = (tmp.x.array >= 0.0).astype(self.chi_vol.x.array.dtype)

    # 3) Stress "linearizzato" per il sottoproblema meccanico (usa χ_vol noto)
    def sigma_damage_linear(self, u_trial, d, material):
        """
        σ = g(d) * [ K χ_vol tr(eps(u_trial)) I + 2G dev(eps(u_trial)) ]
            + (1-χ_vol) K tr(eps(u_trial)) I
        -> forma lineare in u_trial (χ_vol è un coefficiente noto dall'iterazione precedente).
        """
        eps = self.epsilon(u_trial)
        tr  = ufl.tr(eps)
        dev = self._dev(eps)

        lmbda = material["lmbda"]
        G     = material["G"]
        K     = lmbda + 2.0*G/3.0

        I = self._I()
        chi = self.chi_vol  # Function(V_d) noto

        return self.g_d(d) * (K * chi * tr * I + 2.0 * G * dev) + (1.0 - chi) * K * tr * I

    def sigma_damage(self, u, d, material):
        """
        σ = g(d) [K <tr ε>_+ I + 2G dev ε] + [K <tr ε>_- I]
        """
        eps = self.epsilon(u)
        tr  = ufl.tr(eps)
        dev = self._dev(eps)

        lmbda = material["lmbda"]
        G     = material["G"]
        K     = lmbda + 2.0*G/3.0

        tr_p = self._pos(tr)
        tr_m = self._neg(tr)

        return self.g_d(d) * (K*tr_p*self._I() + 2.0*G*dev) + (K*tr_m)*self._I()

    def build_damage_forms(self):
        """
        a_d(wd, eta) = (Gc/cw) * ( 1/ell * wd*eta + ell * grad wd · grad eta ) + 2 H wd eta
        L_d(eta)     = 2 H eta
        (AT-2 con history field H; tutte le forme sono integrate su dx!)
        """
        # Misura di integrazione sul dominio
        dx = ufl.Measure("dx", domain=self.mesh)

        # Parte bulk + gradiente (AT-2) — BILINEARE
        a_bulk = (self.Gc / self.cw) * (
            (1.0 / self.ell) * self.wd * self.eta
            + self.ell * ufl.dot(ufl.grad(self.wd), ufl.grad(self.eta))
        )

        # Termine con H (linearizzato in wd) — resta nel lato sinistro perché è bilineare in wd*eta
        a_hist = 2.0 * self.H * self.wd * self.eta

        # FORM bilineare e lineare (con *dx!)
        self.a_d_form = dolfinx.fem.form((a_bulk + a_hist) * dx)
        self.L_d_form = dolfinx.fem.form((2.0 * self.H * self.eta) * dx)

        print("[PF] forms ok → a_d (rank={} terms), L_d (rank={} terms)".format(
            self.a_d_form.rank, self.L_d_form.rank))

    def update_history(self, u):
        """
        Aggiorna H = max(H, ψ⁺(u)) su tutta la mesh (irreversibilità).
        """
        psi_plus_global = dolfinx.fem.Function(self.V_d)
        # somma ψ⁺ sui materiali (se più materiali)
        psi_plus_global.x.array[:] = 0.0
        for name, mat in self.materials.items():
            psi_plus, _ = self.psi_split(u, mat)
            expr = dolfinx.fem.Expression(psi_plus, self.V_d.element.interpolation_points())
            tmp  = dolfinx.fem.Function(self.V_d)
            tmp.interpolate(expr)
            psi_plus_global.x.array[:] += tmp.x.array

        # H = max(H, ψ⁺)
        self.H.x.array[:] = np.maximum(self.H.x.array, psi_plus_global.x.array)

    def assemble_solve_damage(self):
        """
        Risolve il sottoproblema lineare di d (AT-2) su V_d:
        (Gc/cw)( 1/ell d eta + ell grad d·grad eta ) + 2 H d eta = 2 H eta
        con eventuali Dirichlet su d in self.dirichlet_damage.
        """
        from dolfinx.fem import petsc as fem_petsc
        from petsc4py import PETSc
        import numpy as np

        # Assembla matrice e vettore
        A = fem_petsc.assemble_matrix(self.a_d_form, bcs=self.dirichlet_damage)
        A.assemble()
        b = fem_petsc.assemble_vector(self.L_d_form)
        fem_petsc.apply_lifting(b, [self.a_d_form], [self.dirichlet_damage])
        b.ghostUpdate(addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE)
        fem_petsc.set_bc(b, self.dirichlet_damage)

        # KSP
        ksp = PETSc.KSP().create(self.mesh.comm)
        ksp.setOperators(A)
        ksp.setType("cg")
        pc = ksp.getPC()
        pc.setType("hypre")  # o "gamg"
        ksp.setFromOptions()

        # Soluzione
        uh = dolfinx.fem.Function(self.V_d, name="Damage_tmp")
        ksp.solve(b, uh.x.petsc_vec)         # <-- QUI il fix: usa .x.petsc_vec
        uh.x.scatter_forward()

        # Clamp 0 ≤ d ≤ 1
        np.clip(uh.x.array, 0.0, 1.0, out=uh.x.array)
        self.d.x.array[:] = uh.x.array

