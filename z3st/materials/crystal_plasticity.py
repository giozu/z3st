# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.4.1 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
Small-strain, rate-dependent crystal plasticity on multiple slip systems.

Card keys (with ``constitutive: custom`` and ``plasticity.mode: custom``)::

    stress_function: materials.crystal_plasticity.crystal_stress
    crystal: fcc              # fcc: {111}<110> (12); bcc: {110}<111> + {112}<111> (24)
    euler_angles: [0, 0, 0]   # Bunge (phi1, Phi, phi2), degrees
    g0: 200.0e+6              # initial slip resistance (Pa), same on every system
    gamma0: 0.001             # reference slip rate (1/s)
    n_pow: 5.0                # rate exponent

Model::

    tau_a   = sigma_tr : P_a,         sigma_tr = C : (eps - ep_n)
    dgamma_a/dt = gamma0 * (tau_a / g0) * |tau_a / g0|^(n - 1)
    ep      = ep_n + dt * sum_a dgamma_a/dt * P_a
    sigma   = C : (eps - ep)
    p       = p_n + dt * sum_a |dgamma_a/dt|        (cumulative slip)

Two time integrations:

- explicit (default): the slip rates are evaluated on the trial stress, so
  the update is implicit in the total strain and explicit in the plastic
  strain. Pure UFL; the tangent comes from UFL differentiation of ``sigma``.
  First-order in time: under a constant strain rate the steady stress lies
  below the exact saturation stress by about E * (strain increment per step),
  and the update needs dt * n * E * (strain rate) / stress below about 1.
- implicit (backward Euler), selected by adding to the card::

      stress_operator: materials.crystal_plasticity.crystal_stress_operator

  The slip rates are evaluated on the end-of-step stress. The local problem
  is solved by Newton at every quadrature point in numpy and enters the weak
  form as an external operator (``dolfinx-external-operator``), with the
  consistent tangent d sigma / d eps = lambda I(x)I + 2 mu J^-1, J being the
  Jacobian of the local problem. The steady stress is the exact saturation
  stress for any time step.

  Slip-system hardening (implicit law only), selected by ``hardening``.

  ``hardening: saturation`` (phenomenological)::

      h0: 5.0e+9        # hardening modulus (Pa)
      g_sat: 300.0e+6   # saturation slip resistance (Pa)
      a_hard: 2.25      # exponent
      q_latent: 1.4     # latent-hardening ratio (1 on the diagonal)

      dg_a/dt = sum_b q_ab h0 (1 - g_b/g_sat)^a_hard |dgamma_b/dt|,   g_a(0) = g0

  ``hardening: dislocation_density`` (Taylor resistance, Kocks-Mecking
  evolution; ``g0`` is not used)::

      tau0: 20.0e+6        # lattice friction (Pa)
      alpha_taylor: 0.3    # Taylor coefficient
      burgers: 2.5e-10     # Burgers vector (m)
      rho0: 1.0e+11        # initial dislocation density per system (1/m2)
      k1: 1.3e+8           # storage coefficient (1/m)
      k2: 10.0             # dynamic-recovery coefficient
      a_latent: 1.0        # interaction between different systems (1 on the diagonal)

      F_a         = sum_b A_ab rho_b
      g_a         = tau0 + alpha_taylor G burgers sqrt(F_a)
      drho_a/dt   = (k1 sqrt(F_a) - k2 rho_a) |dgamma_a/dt|

  The hardening variables (g_a or rho_a) join the plastic strain as unknowns
  of the local solve (9 + number of systems per point).

Isotropic elasticity, no temperature dependence.

``uniaxial_reference`` integrates either scheme at a material point with
numpy and scipy, for verification. ``python crystal_plasticity.py`` runs the self-check.
"""

import itertools

import numpy as np

# --. slip systems --..
_FAMILIES = {
    "fcc": [((1, 1, 1), (1, 1, 0))],
    "bcc": [((1, 1, 0), (1, 1, 1)), ((1, 1, 2), (1, 1, 1))],
}


def _family(v):
    """All distinct directions of the cubic family of v, up to sign."""
    out = set()
    for perm in itertools.permutations(v):
        for signs in itertools.product((1, -1), repeat=3):
            w = tuple(a * s for a, s in zip(perm, signs))
            if tuple(-x for x in w) not in out:
                out.add(w)
    return sorted(out)


def schmid_tensors(crystal="fcc", euler_angles=(0.0, 0.0, 0.0)):
    """Symmetric Schmid tensors P_a = sym(m (x) n) in the sample frame, shape (N, 3, 3)."""
    phi1, Phi, phi2 = np.radians(euler_angles)

    def rz(a):
        return np.array([[np.cos(a), np.sin(a), 0], [-np.sin(a), np.cos(a), 0], [0, 0, 1]])

    rx = np.array([[1, 0, 0], [0, np.cos(Phi), np.sin(Phi)], [0, -np.sin(Phi), np.cos(Phi)]])
    g = rz(phi2) @ rx @ rz(phi1)        # Bunge: crystal = g @ sample

    P = []
    for plane, direction in _FAMILIES[crystal.lower()]:
        for n in _family(plane):
            for m in _family(direction):
                if np.dot(n, m) == 0:
                    n_s = g.T @ (np.array(n) / np.linalg.norm(n))
                    m_s = g.T @ (np.array(m) / np.linalg.norm(m))
                    P.append(0.5 * (np.outer(m_s, n_s) + np.outer(n_s, m_s)))
    return np.array(P)


def _card(material):
    P = schmid_tensors(material.get("crystal", "fcc"), material.get("euler_angles", (0.0, 0.0, 0.0)))
    # g0 is absent with dislocation-density hardening (the resistance follows from rho0)
    return P, float(material["gamma0"]), float(material.get("g0", np.nan)), float(material["n_pow"])


# --. UFL law --..
def _update(u, material, model):
    """UFL expressions (eps, ep_new, p_new) for the current displacement."""
    import ufl

    P, gamma0, g0, n_pow = _card(material)
    if material.get("hardening", "none") != "none" and "stress_operator" not in material:
        raise ValueError("crystal plasticity: hardening needs the implicit law (stress_operator)")
    lmbda, mu = material["lmbda"], material["G"]

    # No time step exists before the first solve, when the output expressions
    # are compiled once: there dt = 0, so the written stress is C : (eps - ep_n)
    # with ep_n already updated, i.e. the converged stress.
    # With the implicit operator this function only serves the output, for
    # the same reason with dt = 0.
    dt = 0.0 if "stress_operator" in material else getattr(model, "dt", 0.0)

    eps = model.epsilon(u)
    eps_el_tr = eps - model.ep_n
    sigma_tr = lmbda * ufl.tr(eps_el_tr) * ufl.Identity(3) + 2 * mu * eps_el_tr

    dep = 0
    dp = 0
    for P_a in P:
        s = ufl.inner(sigma_tr, ufl.as_matrix(P_a.tolist())) / g0
        # s * |s|^(n-1), written so that its derivative is finite at s = 0
        rate = gamma0 * s * (s * s + 1e-30) ** ((n_pow - 1) / 2)
        dep += dt * rate * ufl.as_matrix(P_a.tolist())
        dp += dt * abs(rate)
    return eps, model.ep_n + dep, model.p_n + dp


def crystal_stress(u, T, material, model=None):
    """Cauchy stress, for ``stress_function``."""
    import ufl

    eps, ep_new, _ = _update(u, material, model)
    eps_el = eps - ep_new
    return material["lmbda"] * ufl.tr(eps_el) * ufl.Identity(3) + 2 * material["G"] * eps_el


def get_cp_internal_variables(u, T, material, model=None):
    """Updated (plastic strain tensor, cumulative slip), for PlasticityModel."""
    if "stress_operator" in material:
        # Called once per accepted step: state at the converged u, then the
        # hardening variables advance (ep_n and p_n are advanced by the caller).
        model._evaluate_external_operators()
        ep_new, p_new, z_new = model._cp_state[material["name"]]
        model.slip_state_n[material["name"]].x.array[:] = z_new.x.array
        return ep_new, p_new
    _, ep_new, p_new = _update(u, material, model)
    return ep_new, p_new


# --. implicit law (external operator) --..
def _interaction(n_sys, latent):
    A = np.full((n_sys, n_sys), float(latent))
    np.fill_diagonal(A, 1.0)
    return A


def _hardening(material, n_sys, mu):
    """Hardening law of the card: (z0, law).

    z is the internal variable per slip system (z0 its initial value) and
    ``law(z, absr)`` returns, for z and |slip rate| of shape (N, A):

        g    (N, A)     slip resistance
        G    (N, A, A)  d g_a / d z_c
        M    (N, A, A)  dz_a/dt = sum_b M_ab |rate_b|
        D    (N, A, A)  d (M |rate|)_a / d z_c at fixed rates
    """
    kind = material.get("hardening", "none")
    eye = np.eye(n_sys)

    if kind in ("none", "saturation"):
        # z = g:  dg_a/dt = sum_b q_ab h0 (1 - g_b/g_sat)^a_hard |rate_b|
        h0 = float(material["h0"]) if kind == "saturation" else 0.0
        g_sat = float(material.get("g_sat", material["g0"]))
        a_hard = float(material.get("a_hard", 1.0))
        q = _interaction(n_sys, material.get("q_latent", 1.0))

        def law(z, absr):
            if h0 > 0:
                x = 1 - z / g_sat
                H = h0 * np.abs(x) ** a_hard * np.sign(x)
                dH = -h0 * a_hard * np.abs(x) ** (a_hard - 1) / g_sat
            else:
                H = dH = np.zeros_like(z)
            return z, np.broadcast_to(eye, z.shape + (n_sys,)), q[None] * H[:, None, :], q[None] * (dH * absr)[:, None, :]

        return float(material["g0"]), law

    if kind == "dislocation_density":
        # z = rho:  g_a = tau0 + alpha mu b sqrt(F_a),  F_a = sum_b A_ab rho_b   (Taylor)
        #           drho_a/dt = (k1 sqrt(F_a) - k2 rho_a) |rate_a|               (Kocks-Mecking)
        tau0, k1, k2 = float(material["tau0"]), float(material["k1"]), float(material["k2"])
        c = float(material["alpha_taylor"]) * mu * float(material["burgers"])
        A = _interaction(n_sys, material.get("a_latent", 1.0))

        def law(z, absr):
            root_F = np.sqrt(z @ A.T)
            droot = A[None] / (2 * root_F[:, :, None])           # d sqrt(F_a) / d rho_c
            K = k1 * root_F - k2 * z
            return (tau0 + c * root_F, c * droot, eye[None] * K[:, :, None],
                    absr[:, :, None] * (k1 * droot - k2 * eye[None]))

        return float(material["rho0"]), law

    raise ValueError(f"crystal plasticity: unknown hardening '{kind}'")


def _implicit_update(eps, ep_n, z_n, P, gamma0, n_pow, lmbda, mu, dt, z_ref, law):
    """Backward-Euler local solve, batched over points.

    Unknowns per point: the plastic strain (9) and the hardening variables
    (A, in units of z_ref). eps, ep_n: (N, 3, 3); z_n: (N, A). Returns
    (sigma, tangent, ep, z, slip increment) with shapes (N, 3, 3),
    (N, 3, 3, 3, 3), (N, 3, 3), (N, A), (N,).
    """
    N, A = eps.shape[0], len(P)
    Pf = P.reshape(A, 9)
    e, ep_n = eps.reshape(N, 9), ep_n.reshape(N, 9)

    def newton(ep, z, dt):
        """Newton from (ep, z); returns None if it does not converge."""
        ep, z = ep.copy(), z.copy()
        J = np.zeros((N, 9 + A, 9 + A))
        for _ in range(50):
            tau = 2 * mu * (e - ep) @ Pf.T                     # P_a is traceless
            s = tau / law(z, tau * 0)[0]
            pw = np.abs(s) ** (n_pow - 1)
            rate = gamma0 * s * pw
            g, G, M, D = law(z, np.abs(rate))
            dr_dtau = gamma0 * n_pow * pw / g
            dr_dg = -n_pow * rate / g
            sgn = np.sign(rate)

            R = np.concatenate([ep - ep_n - dt * rate @ Pf,
                                (z - z_n - dt * np.einsum("nab,nb->na", M, np.abs(rate))) / z_ref], axis=1)
            if not np.isfinite(R).all() or (z <= 0).any():
                return None
            J_ee = 2 * mu * dt * np.einsum("na,ai,aj->nij", dr_dtau, Pf, Pf)
            J_ze = 2 * mu * dt * np.einsum("nab,nb,bi->nai", M, sgn * dr_dtau, Pf) / z_ref
            J[:, :9, :9] = np.eye(9) + J_ee
            J[:, :9, 9:] = -dt * np.einsum("ai,na,nac->nic", Pf, dr_dg, G) * z_ref
            J[:, 9:, :9] = J_ze
            J[:, 9:, 9:] = np.eye(A) - dt * (D + np.einsum("nab,nb,nbc->nac", M, sgn * dr_dg, G))
            if np.abs(R[:, :9]).max() < 1e-14 and np.abs(R[:, 9:]).max() < 1e-12:
                return ep, z, rate, J, J_ee, J_ze
            dx = np.linalg.solve(J, R[..., None])[..., 0]
            ep -= dx[:, :9]
            z -= dx[:, 9:] * z_ref
        return None

    def solve(dt, depth=0):
        # From the last converged state first. When the trial stress is far
        # above the slip resistance the power law defeats that guess: then
        # the solution for half the time step is the guess (continuation in
        # dt; for dt -> 0 the solution is the starting state itself).
        # NOTE: the whole batch is redone if one point fails.
        out = newton(ep_n, z_n, dt)
        if out is None and depth < 30:
            half = solve(dt / 2, depth + 1)
            out = newton(half[0], half[1], dt)
        if out is None:
            raise RuntimeError("crystal plasticity: local Newton not converged")
        return out

    ep, z, rate, J, J_ee, J_ze = solve(dt)

    # d(ep, z)/d eps = J^-1 [J_ee; J_ze], from the implicit function theorem
    dep = np.linalg.solve(J, np.concatenate([J_ee, J_ze], axis=1))[:, :9, :]
    eye = np.eye(3)
    sigma = lmbda * np.trace(eps, axis1=1, axis2=2)[:, None, None] * eye + 2 * mu * (e - ep).reshape(N, 3, 3)
    tangent = lmbda * np.einsum("ij,kl->ijkl", eye, eye) + 2 * mu * (np.eye(9) - dep).reshape(N, 3, 3, 3, 3)
    return sigma, tangent, ep.reshape(N, 3, 3), z, dt * np.abs(rate).sum(axis=1)


def crystal_stress_operator(u, material, model):
    """Cauchy stress as an external operator of the strain, for ``stress_operator``."""
    import basix.ufl
    import dolfinx
    from dolfinx_external_operator import FEMExternalOperator   # already imported by mechanical_model

    P, gamma0, _, n_pow = _card(material)
    E, nu = float(material["E"]), float(material["nu"])
    lmbda, mu = E * nu / ((1 + nu) * (1 - 2 * nu)), E / (2 * (1 + nu))
    z0, law = _hardening(material, len(P), mu)
    name = material["name"]

    # Hardening variables at the last converged step (slip resistance or
    # dislocation density), one per system and quadrature point. In
    # model.slip_state_n so that a failed step rolls them back.
    if not hasattr(model, "slip_state_n"):
        model.slip_state_n = {}
    if name not in model.slip_state_n:
        element = basix.ufl.quadrature_element(
            model.mesh.topology.cell_name(), value_shape=(len(P),), degree=model.q_degree)
        model.slip_state_n[name] = dolfinx.fem.Function(dolfinx.fem.functionspace(model.mesh, element))
        model.slip_state_n[name].x.array[:] = z0
    z_n = model.slip_state_n[name]

    # End-of-step state of the last evaluation, read by PlasticityModel.
    ep_new = dolfinx.fem.Function(model.V_pl)
    p_new = dolfinx.fem.Function(model.Q_pl)
    z_new = dolfinx.fem.Function(z_n.function_space)
    if not hasattr(model, "_cp_state"):
        model._cp_state = {}
    model._cp_state[name] = (ep_new, p_new, z_new)

    def solve(eps):
        # NOTE: solved once for the residual and once for the tangent;
        # share the result if the local solve ever shows up in a profile.
        sigma, tangent, ep, z, dp = _implicit_update(
            eps.reshape(-1, 3, 3), model.ep_n.x.array.reshape(-1, 3, 3), z_n.x.array.reshape(-1, len(P)),
            P, gamma0, n_pow, lmbda, mu, model.dt, z0, law)
        ep_new.x.array[:] = ep.reshape(-1)
        p_new.x.array[:] = model.p_n.x.array + dp
        z_new.x.array[:] = z.reshape(-1)
        return sigma, tangent

    def external(derivatives):
        if derivatives == (0,):
            return lambda eps: solve(eps)[0].reshape(-1)
        if derivatives == (1,):
            return lambda eps: solve(eps)[1].reshape(-1)
        raise NotImplementedError(derivatives)

    return FEMExternalOperator(model.epsilon(u), function_space=model.V_pl, external_function=external)


# --. material-point reference (numpy) --..
def uniaxial_reference(material, eps_zz, dt, implicit=False):
    """Either scheme at a material point under uniaxial stress along z.

    ``eps_zz`` is the axial strain at each step; the other five strain
    components are found so that the other five stress components vanish.
    With ``implicit`` the plastic strain is a further unknown of the same
    scipy root solve (it does not use ``_implicit_update``).
    Returns (sigma_zz, cumulative slip) per step.
    """
    from scipy.optimize import root

    P, gamma0, g0, n_pow = _card(material)
    E, nu = float(material["E"]), float(material["nu"])
    lmbda, mu = E * nu / ((1 + nu) * (1 - 2 * nu)), E / (2 * (1 + nu))
    free = [(0, 0), (1, 1), (0, 1), (0, 2), (1, 2)]

    def hooke(e):
        return lmbda * np.trace(e) * np.eye(3) + 2 * mu * e

    def step(x, ezz, ep_n):
        eps = np.zeros((3, 3))
        eps[2, 2] = ezz
        for (i, j), v in zip(free, x):
            eps[i, j] = eps[j, i] = v
        s = np.einsum("aij,ij->a", P, hooke(eps - ep_n)) / g0
        rate = gamma0 * s * np.abs(s) ** (n_pow - 1)
        ep = ep_n + dt * np.einsum("a,aij->ij", rate, P)
        return hooke(eps - ep), ep, dt * np.abs(rate).sum()

    # Hardening written out again here, on purpose independent of _hardening.
    kind = material.get("hardening", "none")
    assert implicit or kind == "none"
    A = len(P)
    if kind == "dislocation_density":
        z0 = float(material["rho0"])
        inter = np.full((A, A), float(material.get("a_latent", 1.0)))
        np.fill_diagonal(inter, 1.0)
        c = float(material["alpha_taylor"]) * mu * float(material["burgers"])

        def resistance(z):
            return float(material["tau0"]) + c * np.sqrt(np.maximum(inter @ z, 0))   # root may try z < 0

        def z_dot(z, absr):
            return (float(material["k1"]) * np.sqrt(np.maximum(inter @ z, 0)) - float(material["k2"]) * z) * absr
    else:
        z0 = g0
        inter = np.full((A, A), float(material.get("q_latent", 1.0)))
        np.fill_diagonal(inter, 1.0)
        h0 = float(material["h0"]) if kind == "saturation" else 0.0

        def resistance(z):
            return z

        def z_dot(z, absr):
            x = 1 - z / float(material.get("g_sat", g0))
            return inter @ (h0 * np.abs(x) ** float(material.get("a_hard", 1.0)) * np.sign(x) * absr)

    def step_implicit(y, ezz, ep_n, z_n, dt=dt):
        eps = np.zeros((3, 3))
        eps[2, 2] = ezz
        for (i, j), v in zip(free, y[:5]):
            eps[i, j] = eps[j, i] = v
        ep, z = y[5:14].reshape(3, 3), y[14:] * z0
        sigma = hooke(eps - ep)
        s = np.einsum("aij,ij->a", P, sigma) / resistance(z)
        rate = gamma0 * s * np.abs(s) ** (n_pow - 1)
        res = np.concatenate([[sigma[i, j] / E for i, j in free],
                              (ep - ep_n - dt * np.einsum("a,aij->ij", rate, P)).ravel(),
                              (z - z_n - dt * z_dot(z, np.abs(rate))) / z0])
        return res, sigma, ep, z, dt * np.abs(rate).sum()

    ep_n, z_n, p, x = np.zeros((3, 3)), np.full(A, z0), 0.0, np.zeros(5)
    sigma_zz, slip = [], []
    for ezz in eps_zz:
        if implicit:
            y0 = np.concatenate([x, ep_n.ravel(), z_n / z0])
            # root can stall when the step starts far above the slip resistance:
            # then ramp the time step up to dt, each solution the next guess.
            for ramp in ([dt], dt * 0.5 ** np.arange(20, -1, -1)):
                y = y0
                for dt_k in ramp:
                    y = root(lambda y: step_implicit(y, ezz, ep_n, z_n, dt_k)[0], y, tol=1e-14).x
                if np.abs(step_implicit(y, ezz, ep_n, z_n)[0]).max() < 1e-10:
                    break
            else:
                raise RuntimeError(f"uniaxial_reference: no solution at eps_zz = {ezz}")
            _, sigma, ep_n, z_n, dp = step_implicit(y, ezz, ep_n, z_n)
            x = y[:5]
        else:
            x = root(lambda y: [step(y, ezz, ep_n)[0][i, j] / E for i, j in free], x, tol=1e-12).x
            sigma, ep_n, dp = step(x, ezz, ep_n)
        p += dp
        sigma_zz.append(sigma[2, 2])
        slip.append(p)
    return np.array(sigma_zz), np.array(slip)


def saturation_stress(material, strain_rate):
    """Steady uniaxial stress along z, with the hardening variables at their steady value.

    strain_rate = gamma0 * sum_a |m_a|^(n+1) * (sigma / g_a)^n, with g_a = g0
    (no hardening), g_sat (saturation law) or the Taylor resistance at the
    steady densities k1 sqrt(F_a) = k2 rho_a of the systems that slip.
    """
    P, gamma0, g0, n_pow = _card(material)
    m = np.abs(P[:, 2, 2])
    kind = material.get("hardening", "none")
    g = np.full(len(P), g0)
    if kind == "saturation":
        g[:] = float(material["g_sat"])
    elif kind == "dislocation_density":
        mu = float(material["E"]) / (2 * (1 + float(material["nu"])))
        A = _interaction(len(P), material.get("a_latent", 1.0))
        rho = np.full(len(P), float(material["rho0"]))
        for _ in range(500):                                  # fixed point, a contraction
            rho[m > 1e-12] = (float(material["k1"]) * np.sqrt(A @ rho) / float(material["k2"]))[m > 1e-12]
        g = float(material["tau0"]) + float(material["alpha_taylor"]) * mu * float(material["burgers"]) * np.sqrt(A @ rho)
    return (strain_rate / (gamma0 * (m ** (n_pow + 1) / g ** n_pow).sum())) ** (1 / n_pow)


def _self_check():
    def factors(crystal, euler=(0, 0, 0)):
        m = np.abs(schmid_tensors(crystal, euler)[:, 2, 2])
        return sorted(np.round(m[m > 1e-12], 4).tolist())

    r6 = round(1 / np.sqrt(6), 4)
    # fcc along [001]: 8 systems at 1/sqrt(6)
    assert factors("fcc") == [r6] * 8
    # fcc along [111]: 6 systems at 2/(3 sqrt(6))
    to_111 = (0.0, np.degrees(np.arccos(1 / np.sqrt(3))), 45.0)
    assert factors("fcc", to_111) == [round(2 / (3 * np.sqrt(6)), 4)] * 6
    # bcc along [001]: {110}: 8 at 1/sqrt(6); {112}: 8 at 1/(3 sqrt(2)), 4 at sqrt(2)/3
    assert len(schmid_tensors("bcc")) == 24
    assert factors("bcc") == [round(1 / (3 * np.sqrt(2)), 4)] * 8 + [r6] * 8 + [round(np.sqrt(2) / 3, 4)] * 4

    # the material-point curve saturates below the exact value by about
    # E * (strain increment), and the gap halves with the time step
    mat = dict(crystal="fcc", E=200e9, nu=0.3, g0=200e6, gamma0=1e-3, n_pow=5.0)
    rate, gaps = 5e-3, []
    for dt in (0.01, 0.005):
        sigma, slip = uniaxial_reference(mat, rate * dt * np.arange(1, int(5 / dt) + 1), dt)
        gaps.append(saturation_stress(mat, rate) - sigma[-1])
        assert 0.5 < gaps[-1] / (mat["E"] * rate * dt) < 1.5, gaps
        assert slip[-1] > 0
    assert abs(gaps[0] / gaps[1] - 2) < 0.05, gaps

    # implicit scheme: the steady stress is the exact saturation stress
    dt = 0.05
    sigma, _ = uniaxial_reference(mat, rate * dt * np.arange(1, 201), dt, implicit=True)
    assert abs(sigma[-1] / saturation_stress(mat, rate) - 1) < 1e-8, sigma[-1]

    # with hardening the steady stress is the saturation stress at g_sat
    # (a_hard = 1: exponential approach to g_sat, complete within the run)
    hard = dict(mat, hardening="saturation", h0=2e10, g_sat=300e6, a_hard=1.0, q_latent=1.4)
    sigma, _ = uniaxial_reference(hard, rate * dt * np.arange(1, 801), dt, implicit=True)
    assert abs(sigma[-1] / saturation_stress(hard, rate) - 1) < 1e-6, sigma[-1]
    assert abs(saturation_stress(hard, rate) / saturation_stress(mat, rate) - 1.5) < 1e-12

    # dislocation densities: steady state k1 sqrt(F) = k2 rho. A numerical
    # test of the scheme: the strain needed is far beyond small strain.
    dd = dict(crystal="fcc", E=200e9, nu=0.3, gamma0=1e-3, n_pow=5.0, hardening="dislocation_density",
              tau0=20e6, alpha_taylor=0.3, burgers=2.5e-10, rho0=1e11, k1=1.3e8, k2=50.0, a_latent=1.0)
    sigma, _ = uniaxial_reference(dd, rate * 1.0 * np.arange(1, 801), 1.0, implicit=True)
    assert abs(sigma[-1] / saturation_stress(dd, rate) - 1) < 1e-6, (sigma[-1], saturation_stress(dd, rate))
    # closed form with a_latent = 1 along [001]: rho = (k1/k2)^2 (4 + sqrt(16 + 4 (k2/k1)^2 rho0))
    rho = (1.3e8 / 50) ** 2 * (4 + np.sqrt(16 + 4 * (50 / 1.3e8) ** 2 * 1e11))
    g = 20e6 + 0.3 * (200e9 / 2.6) * 2.5e-10 * np.sqrt(8 * rho + 4 * 1e11)
    assert abs(saturation_stress(dd, rate) / saturation_stress(dict(mat, g0=g), rate) - 1) < 1e-10

    # implicit local solve: tangent against central finite differences, for
    # the three hardening options, from a hardened state
    rng = np.random.default_rng(0)
    eps = rng.normal(scale=3e-3, size=(1, 3, 3))
    eps = 0.5 * (eps + eps.transpose(0, 2, 1))
    P = schmid_tensors("bcc", (10, 20, 30))
    for card, lo, hi in [(dict(mat), 200e6, 280e6),
                         (dict(mat, hardening="saturation", h0=5e9, g_sat=300e6, a_hard=2.25, q_latent=1.4), 200e6, 280e6),
                         (dict(dd, a_latent=0.5), 1e12, 5e13)]:
        z_ref, law = _hardening(card, len(P), 77e9)
        z_n = rng.uniform(lo, hi, size=(1, len(P)))
        args = (np.zeros((1, 3, 3)), z_n, P, 1e-3, 5.0, 115e9, 77e9, 0.05, z_ref, law)
        tangent = _implicit_update(eps, *args)[1][0]
        h = 1e-8
        for k in range(3):
            for l in range(3):
                d = np.zeros((1, 3, 3))
                d[0, k, l] = h
                fd = (_implicit_update(eps + d, *args)[0][0] - _implicit_update(eps - d, *args)[0][0]) / (2 * h)
                assert np.abs(fd - tangent[:, :, k, l]).max() < 1e-5 * np.abs(tangent).max(), (card.get("hardening"), k, l)
    print("crystal_plasticity self-check: OK")


if __name__ == "__main__":
    _self_check()
