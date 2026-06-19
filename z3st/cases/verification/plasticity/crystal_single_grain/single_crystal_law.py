import ufl
import numpy as np
import dolfinx

# Global flag to print Schmid factor only once
_schmid_printed = False

def single_crystal_stress(u, T, material, model=None):
    """
    Single Crystal Plasticity with implicit time integration.
    Uses Automatic Differentiation (AD) on slip-system mapping.

    Rate-dependent formulation (viscoplastic):
    ------------------------------------------
    γ̇ = γ₀ |τ/g₀|^n sign(τ)  (power law)

    For each time step Δt, the accumulated plastic slip is:
    Δγ = Δt * γ̇

    This gives plastic strain increment:
    Δε_p = Δγ * P = Δt * γ̇ * P

    The equilibrium is solved implicitly, so γ̇ is computed at the current state.
    This is equivalent to backward Euler integration.
    """

    # -- 1. Slip Systems (Favorably oriented for Z-axis tension)
    # Slip system: (111)[0-11] - favorably oriented for z-axis loading
    # For FCC: slip plane {111}, slip direction <110>
    n_alpha = ufl.as_vector([1, 1, 1]) / ufl.sqrt(3)   # (111) plane normal
    m_alpha = ufl.as_vector([0, 1, -1]) / ufl.sqrt(2)  # [0-11] slip direction

    # Calculate and print Schmid factor for z-axis loading (only once)
    global _schmid_printed
    if not _schmid_printed:
        # Schmid factor = (m·e_z)(n·e_z) for uniaxial loading along z
        n_np = np.array([1, 1, 1]) / np.sqrt(3)
        m_np = np.array([0, 1, -1]) / np.sqrt(2)
        e_z = np.array([0, 0, 1])

        # Method 1: Direct calculation for uniaxial z-loading
        schmid_factor_z = abs(np.dot(m_np, e_z) * np.dot(n_np, e_z))

        # Method 2: Using Schmid tensor P_zz component
        # P = sym(m ⊗ n), so P_zz gives the projection
        P_np = 0.5 * (np.outer(m_np, n_np) + np.outer(n_np, m_np))
        schmid_tensor_zz = abs(P_np[2, 2])  # z-component

        # Method 3: Full tensor contraction with uniaxial stress
        sigma_unit_z = np.zeros((3, 3))
        sigma_unit_z[2, 2] = 1.0  # Unit stress in z-direction
        schmid_full = abs(np.sum(P_np * sigma_unit_z))  # τ = P : σ

        print("\n" + "="*70)
        print("CRYSTAL PLASTICITY - SCHMID FACTOR CALCULATION")
        print("="*70)
        print(f"Slip system: (111)[0-11]")
        print(f"  Plane normal n = {n_np}")
        print(f"  Slip direction m = {m_np}")
        print(f"\nLoading direction: e_z = [0, 0, 1]")
        print(f"\nSchmid factor calculations:")
        print(f"  Method 1 (direct):        μ = |m·e_z| × |n·e_z| = {schmid_factor_z:.6f}")
        print(f"  Method 2 (tensor P_zz):   μ = |P_zz|           = {schmid_tensor_zz:.6f}")
        print(f"  Method 3 (full tensor):   μ = |P:σ_zz|         = {schmid_full:.6f}")
        print(f"\nSchmid tensor P:")
        for i in range(3):
            print(f"  {P_np[i, :]}")
        print(f"\nFor uniaxial stress σ_zz:")
        print(f"  Resolved shear stress: τ = {schmid_full:.6f} × σ_zz")
        print("="*70 + "\n")

        _schmid_printed = True

    # Schmid tensor P = sym(m ⊗ n)
    P = 0.5 * (ufl.outer(m_alpha, n_alpha) + ufl.outer(n_alpha, m_alpha))

    # -- 2. Elasticity
    E = material.get("E", 200e9)
    nu = material.get("nu", 0.3)
    lmbda = E * nu / ((1 + nu) * (1 - 2*nu))
    mu = E / (2 * (1 + nu))

    def epsilon(u):
        return ufl.sym(ufl.grad(u))

    eps = epsilon(u)
    dim = 3

    # -- 3. Get plastic strain history from previous time step
    if model is not None and hasattr(model, 'ep_n'):
        # Use plastic strain from previous converged solution
        eps_p_old = model.ep_n
    else:
        # Fallback: no history available, start from zero
        eps_p_old = ufl.as_tensor([[0, 0, 0], [0, 0, 0], [0, 0, 0]])

    # -- 4. Get time step size from solver
    if model is not None and hasattr(model, 'dt'):
        dt = model.dt
    else:
        dt = 0.05  # Fallback value

    # -- 5. Plasticity Parameters
    gamma0 = float(material.get("gamma0", 0.001))  # Reference slip rate (1/s)
    g0 = float(material.get("g0", 200e6))          # Slip resistance (Pa)
    n_pow = float(material.get("n_pow", 5))        # Power law exponent

    # -- 6. Trial Elastic Stress (using old plastic strain)
    # σ_trial = C:(ε_total - ε_p^n)
    eps_elastic_trial = eps - eps_p_old
    sigma_trial = lmbda * ufl.tr(eps_elastic_trial) * ufl.Identity(dim) + 2 * mu * eps_elastic_trial

    # Resolved shear stress on slip system
    tau = ufl.inner(sigma_trial, P)

    # -- 7. Slip Rate (Power Law with Automatic Differentiation)
    # γ̇ = γ₀ |τ/g₀|ⁿ sign(τ)
    tau_var = ufl.variable(tau)
    gamma_dot = gamma0 * (abs(tau_var/g0))**n_pow * ufl.sign(tau_var)

    # -- 8. Plastic Strain Rate
    # ε̇_p = γ̇ · P
    eps_p_dot = gamma_dot * P

    # -- 9. Time Integration (Backward Euler)
    # ε_p^{n+1} = ε_p^n + Δt · ε̇_p^{n+1}
    eps_p_new = eps_p_old + dt * eps_p_dot

    # -- 10. Final Stress (using NEW plastic strain)
    # σ = C:(ε_total - ε_p^{n+1})
    sigma = lmbda * ufl.tr(eps - eps_p_new) * ufl.Identity(dim) + 2 * mu * (eps - eps_p_new)

    return sigma


def get_cp_internal_variables(u, T, material, model=None):
    """
    Updated plastic strain tensor for crystal plasticity.

    Mirrors single_crystal_stress(); called by PlasticityModel after
    convergence to update ep_n -> ep.

    Parameters:
        u: Displacement field
        T: Temperature field (not used in this model)
        material: Material dictionary with parameters
        model: PlasticityModel instance (provides ep_n, dt)

    Returns:
        eps_p_new: Updated plastic strain tensor (3x3)

    Note:
        The cumulative plastic strain p is calculated in PlasticityModel as:
            p = sqrt(1.5 * ep_new : ep_new)

        For single-slip crystal plasticity, this gives p ≈ 0.866 * γ instead of p = γ.
        This is acceptable since p is not used in the constitutive law (no hardening).

        If hardening is added in the future, p should be computed as:
            p = ∫|γ̇| dt  (cumulative slip on the active slip system)
    """
    import ufl

    # Slip system (same as in single_crystal_stress)
    n_alpha = ufl.as_vector([1, 1, 1]) / ufl.sqrt(3)
    m_alpha = ufl.as_vector([0, 1, -1]) / ufl.sqrt(2)
    P = 0.5 * (ufl.outer(m_alpha, n_alpha) + ufl.outer(n_alpha, m_alpha))

    # Elasticity
    E = material.get("E", 200e9)
    nu = material.get("nu", 0.3)
    lmbda = E * nu / ((1 + nu) * (1 - 2*nu))
    mu = E / (2 * (1 + nu))

    # Strain
    eps = ufl.sym(ufl.grad(u))
    dim = 3

    # Get old plastic strain
    if model is not None and hasattr(model, 'ep_n'):
        eps_p_old = model.ep_n
    else:
        eps_p_old = ufl.as_tensor([[0, 0, 0], [0, 0, 0], [0, 0, 0]])

    # Get time step
    if model is not None and hasattr(model, 'dt'):
        dt = model.dt
    else:
        dt = 0.05

    # Plasticity parameters
    gamma0 = float(material.get("gamma0", 0.001))
    g0 = float(material.get("g0", 200e6))
    n_pow = float(material.get("n_pow", 5))

    # Trial stress
    eps_elastic_trial = eps - eps_p_old
    sigma_trial = lmbda * ufl.tr(eps_elastic_trial) * ufl.Identity(dim) + 2 * mu * eps_elastic_trial
    tau = ufl.inner(sigma_trial, P)

    # Slip rate
    tau_var = ufl.variable(tau)
    gamma_dot = gamma0 * (abs(tau_var/g0))**n_pow * ufl.sign(tau_var)

    # Updated plastic strain
    eps_p_dot = gamma_dot * P
    eps_p_new = eps_p_old + dt * eps_p_dot

    return eps_p_new
