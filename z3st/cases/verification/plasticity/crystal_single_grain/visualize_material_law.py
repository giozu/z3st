#!/usr/bin/env python3
"""
Visualize the crystal plasticity material law and its derivative, motivating
the use of automatic differentiation.
"""

import numpy as np
import matplotlib.pyplot as plt

def plot_slip_law():
    """
    Plot the power-law slip rate and its derivative (strongly nonlinear).
    """

    # Parameters from single_crystal_law.py
    gamma0 = 0.001  # Reference slip rate
    g0 = 200e6      # Slip resistance (Pa)
    n_pow = 5       # Power law exponent

    # Resolved shear stress range
    tau = np.linspace(-150e6, 150e6, 1000)

    # Slip rate
    gamma_dot = gamma0 * np.abs(tau / g0)**n_pow * np.sign(tau)

    # Derivative dgamma_dot/dtau
    dgamma_dtau = gamma0 * n_pow * np.abs(tau / g0)**(n_pow - 1) / g0

    # Create figure
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # ========== Plot 1: Slip Rate Law ==========
    ax1 = axes[0, 0]
    ax1.plot(tau * 1e-6, gamma_dot, linewidth=2.5, color='#2E86AB')
    ax1.axvline(x=g0 * 1e-6, color='red', linestyle='--', linewidth=2,
                label=f'Slip resistance g₀ = {g0*1e-6:.0f} MPa')
    ax1.axvline(x=-g0 * 1e-6, color='red', linestyle='--', linewidth=2)
    ax1.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
    ax1.set_xlabel('Resolved Shear Stress τ (MPa)', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Slip Rate γ̇', fontsize=12, fontweight='bold')
    ax1.set_title(f'Material Law: γ̇ = γ₀ |τ/g₀|ⁿ sign(τ) with n={n_pow}',
                  fontsize=13, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=11)

    # ========== Plot 2: Derivative (Tangent Stiffness) ==========
    ax2 = axes[0, 1]
    ax2.semilogy(tau * 1e-6, np.abs(dgamma_dtau) + 1e-20, linewidth=2.5, color='#E63946')
    ax2.axvline(x=g0 * 1e-6, color='red', linestyle='--', linewidth=2,
                label='Critical stress')
    ax2.axvline(x=-g0 * 1e-6, color='red', linestyle='--', linewidth=2)
    ax2.set_xlabel('Resolved Shear Stress τ (MPa)', fontsize=12, fontweight='bold')
    ax2.set_ylabel('|dγ̇/dτ| (log scale)', fontsize=12, fontweight='bold')
    ax2.set_title('Derivative: What Z3ST Computes Automatically',
                  fontsize=13, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=11)

    # ========== Plot 3: Comparison with Linear Law ==========
    ax3 = axes[1, 0]
    gamma_dot_linear = gamma0 * tau / g0  # Linear approximation
    ax3.plot(tau * 1e-6, gamma_dot, linewidth=2.5, color='#2E86AB', label=f'Power law (n={n_pow})')
    ax3.plot(tau * 1e-6, gamma_dot_linear, '--', linewidth=2, color='gray', alpha=0.7, label='Linear (n=1)')
    ax3.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
    ax3.set_xlabel('Resolved Shear Stress τ (MPa)', fontsize=12, fontweight='bold')
    ax3.set_ylabel('Slip Rate γ̇', fontsize=12, fontweight='bold')
    ax3.set_title('Power Law vs Linear Comparison', fontsize=13, fontweight='bold')
    ax3.grid(True, alpha=0.3)
    ax3.legend(fontsize=11)
    ax3.set_xlim([-150, 150])

    # ========== Plot 4: Derivative Variation ==========
    ax4 = axes[1, 1]
    # Show derivative variation over stress range
    tau_positive = tau[tau > 0]
    dgamma_positive = gamma0 * n_pow * (tau_positive / g0)**(n_pow - 1) / g0
    ax4.loglog(tau_positive * 1e-6, dgamma_positive, linewidth=2.5, color='#E63946')
    ax4.axvline(x=g0 * 1e-6, color='red', linestyle='--', linewidth=2, label=f'g₀ = {g0*1e-6:.0f} MPa')
    ax4.set_xlabel('Resolved Shear Stress τ (MPa)', fontsize=12, fontweight='bold')
    ax4.set_ylabel('dγ̇/dτ (log-log)', fontsize=12, fontweight='bold')
    ax4.set_title('Derivative Variation (Why Manual Jacobian Fails)', fontsize=13, fontweight='bold')
    ax4.grid(True, alpha=0.3, which='both')
    ax4.legend(fontsize=11)

    plt.tight_layout()
    plt.savefig('output/material_law_visualization.png', dpi=300, bbox_inches='tight')
    print("\n[INFO] Material law visualization saved to: output/material_law_visualization.png")

    # Print derivative analysis
    print("\n" + "="*70)
    print("MATERIAL LAW COMPLEXITY ANALYSIS")
    print("="*70)
    print(f"\nPower law exponent: n = {n_pow}")
    print(f"Slip resistance: g₀ = {g0*1e-6:.1f} MPa")

    # Calculate derivative at critical stress
    tau_crit = g0
    dgamma_crit = gamma0 * n_pow * (tau_crit / g0)**(n_pow - 1) / g0
    gamma_crit = gamma0 * (tau_crit / g0)**n_pow

    # Calculate derivative at 2x critical stress
    tau_high = 2 * g0
    dgamma_high = gamma0 * n_pow * (tau_high / g0)**(n_pow - 1) / g0
    gamma_high = gamma0 * (tau_high / g0)**n_pow

    print(f"\nAt τ = {tau_crit*1e-6:.1f} MPa (critical stress):")
    print(f"  Slip rate: γ̇ = {gamma_crit:.2e}")
    print(f"  Derivative: dγ̇/dτ = {dgamma_crit:.2e} Pa⁻¹")

    print(f"\nAt τ = {tau_high*1e-6:.1f} MPa (2× critical):")
    print(f"  Slip rate: γ̇ = {gamma_high:.2e}")
    print(f"  Derivative: dγ̇/dτ = {dgamma_high:.2e} Pa⁻¹")

    derivative_ratio = dgamma_high / dgamma_crit
    print(f"\nDerivative variation: {derivative_ratio:.1e}×")
    print("\n→ This extreme variation makes manual Jacobians error-prone!")
    print("→ UFL computes exact derivatives symbolically - no approximations!")
    print("="*70 + "\n")

if __name__ == "__main__":
    plot_slip_law()
