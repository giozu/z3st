#!/usr/bin/env python3
"""
Script to create graphs displaying contact pressure as a function of time
and other variables from history.csv data.
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "output")
HISTORY_CSV = os.path.join(OUT, "history.csv")

def load_history_data():
    """Load data from history.csv file."""
    if not os.path.exists(HISTORY_CSV):
        raise FileNotFoundError(f"{HISTORY_CSV} not found. Run the case first.")
    
    data = np.genfromtxt(HISTORY_CSV, delimiter=",", names=True)
    return data

def f1(t, w0, xi, K, offset=360):
    """
    Réponse temporelle d'un système du second ordre sur-amorti (ξ > 1) à un échelon unitaire,
    avec un offset de 360 jours.

    Formule :
    s(t) = K * (1 - (1 / (2 * sqrt(ξ² - 1))) *
                     (exp(w0 * (ξ - sqrt(ξ² - 1)) * (t - offset)) / (ξ - sqrt(ξ² - 1)) -
                      exp(w0 * (ξ + sqrt(ξ² - 1)) * (t - offset)) / (ξ + sqrt(ξ² - 1)))) * u(t - offset)

    Args:
        t : Temps (scalaire ou tableau NumPy).
        K : Gain statique.
        xi : Coefficient d'amortissement (ξ > 1 pour un système sur-amorti).
        w0 : Pulsation naturelle (ω₀).
        offset : Délai avant le début de la réponse (en jours, par défaut 360).

    Returns:
        s(t) : Réponse temporelle du système.
    """
    t = np.asarray(t)  # Convertir en tableau NumPy
    s = np.zeros_like(t)  # Initialiser le résultat à 0

    # Masque pour les temps >= offset
    mask = t >= offset
    t_relative = t[mask] - offset  # Temps relatif après l'offset

    if np.any(mask):  # Si au moins un élément satisfait la condition
        sqrt_term = np.sqrt(xi**2 - 1)
        denominator = 2 * sqrt_term

        # Calcul des exponentielles
        exp1 = np.exp(-w0 * (xi - sqrt_term) * t_relative)
        exp2 = np.exp(-w0 * (xi + sqrt_term) * t_relative)

        # Calcul des termes fractionnaires
        term1 = exp1 / (xi - sqrt_term)
        term2 = exp2 / (xi + sqrt_term)

        # Combinaison des termes
        s[mask] = K * (1 - (1 / denominator) * (term1 - term2))

    return s

def f2(t, w_n, zeta, K, offset=360):
    """
    Réponse temporelle d'un système du second ordre sous-amorti à un échelon unitaire,
    avec un offset de 360 jours. Cette version est vectorisée pour accepter des tableaux NumPy.

    Args:
        t : Temps (en jours), peut être un scalaire ou un tableau NumPy.
        w_n : Pulsation naturelle (en rad/jour).
        zeta : Coefficient d'amortissement (0 < zeta < 1).
        K : Gain statique (valeur finale de la réponse, en MPa).
        offset : Délai avant le début de la réponse (en jours).

    Returns:
        y(t) : Réponse temporelle (pression en MPa).
    """
    t = np.asarray(t)  # Convertir en tableau NumPy si ce n'est pas déjà le cas
    result = np.zeros_like(t)

    # Masque pour les temps >= offset
    mask = t >= offset
    t_relative = t[mask] - offset

    if np.any(mask):  # Si au moins un élément satisfait la condition
        omega_d = w_n * np.sqrt(1 - zeta**2)
        phi = np.acos(zeta)
        result[mask] = K * (1 - (np.exp(-zeta * w_n * t_relative) / np.sqrt(1 - zeta**2)) * np.sin(omega_d * t_relative + phi))
    return result


def optimize_second_order_response(data):
    """Optimize the parameters of the second-order response function to fit the pressure data."""
    time_days = data["time_days"]
    pressure = data["contact_pressure_MPa"]
    
    # Filter data for times greater than or equal to 360 days
    mask = time_days >= 360
    time_filtered = time_days[mask]
    pressure_filtered = pressure[mask]

    # Initial guess for parameters: [w_n, zeta, K]
    initial_guess = [0.01, 2, 25.0]

    # Least squares optimization
    popt, pcov = curve_fit(
        lambda t, w_n, zeta, K: f1(t, w_n, zeta, K, 602.8571),
        time_filtered,
        pressure_filtered,
        p0=initial_guess,
        maxfev=10000  # Increase the number of iterations if necessary
    )

    # Retrieve optimized parameters
    w_n_opt, zeta_opt, K_opt = popt
    print(f"Optimized parameters:")
    print(f"  ω_n (natural frequency) = {w_n_opt:.6f} rad/day")
    print(f"  ζ (damping ratio) = {zeta_opt:.6f}")
    print(f"  K (static gain) = {K_opt:.6f} MPa") 
    return popt





def plot_pressure_vs_time(data):
    """Create graph with contact pressure on y-axis and time on x-axis."""
    time_days = data["time_days"]
    pressure = data["contact_pressure_MPa"]
    w_n_opt, zeta_opt, K_opt = optimize_second_order_response(data)
    analytical_pressure = f1(time_days, w_n_opt, zeta_opt, K_opt, 602.8571) 
    error=0
    for i in range (len(analytical_pressure)):
        if pressure[i]!=0:
            error+=np.abs(analytical_pressure[i]-pressure[i])/pressure[i]
    print(f"Average relative error between analytical and numerical pressure: {error/len(analytical_pressure):.6f}")
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(time_days, pressure, "o-", color="#C44E52", lw=2, markersize=6)
    ax.plot(time_days, analytical_pressure, "+-", color="gray", lw=2, markersize=6)
    ax.set_xlabel("Time (days)", fontsize=12)
    ax.set_ylabel("Contact Pressure (MPa)", fontsize=12)
    ax.set_title("Contact Pressure vs Time", fontsize=14, fontweight="bold")
    ax.grid(alpha=0.3)
    
    # Highlight the onset of contact (first non-zero pressure)
    contact_onset = next((i-1 for i, p in enumerate(pressure) if p > 0), None)
    if contact_onset is not None:
        ax.axvline(time_days[contact_onset], color="gray", linestyle="--", alpha=0.7)
        ax.text(time_days[contact_onset], pressure.max() * 0.9, 
                f"Contact onset: {time_days[contact_onset]:.1f} days",
                rotation=90, verticalalignment="top", fontsize=10)
    
    fig.tight_layout()
    output_path = os.path.join(OUT, "pressure_vs_time.png")
    fig.savefig(output_path, dpi=150)
    plt.close(fig)
    print(f"Created {output_path}")
    return output_path

def plot_pressure_vs_all_variables(data):
    """Create multiple subplots showing pressure as function of other variables."""
    # Extract variables
    pressure = data["contact_pressure_MPa"]
    burnup_avg = data["burnup_avg_MWdkgU"]
    burnup_max = data["burnup_max_MWdkgU"]
    gap = data["gap_um"]
    t_max = data["T_max_K"]
    t_min = data["T_min_K"]
    
    # Create figure with subplots
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    fig.suptitle("Contact Pressure as Function of Various Variables", 
                 fontsize=16, fontweight="bold")
    
    # Plot 1: Pressure vs Average Burnup
    axes[0, 0].plot(burnup_avg, pressure, "o-", color="#4C72B0", lw=2, markersize=6)
    axes[0, 0].set_xlabel("Average Burnup (MWd/kgU)", fontsize=11)
    axes[0, 0].set_ylabel("Contact Pressure (MPa)", fontsize=11)
    axes[0, 0].set_title("Pressure vs Average Burnup", fontsize=12)
    axes[0, 0].grid(alpha=0.3)
    
    # Plot 2: Pressure vs Maximum Burnup
    axes[0, 1].plot(burnup_max, pressure, "s-", color="#DD8452", lw=2, markersize=6)
    axes[0, 1].set_xlabel("Maximum Burnup (MWd/kgU)", fontsize=11)
    axes[0, 1].set_ylabel("Contact Pressure (MPa)", fontsize=11)
    axes[0, 1].set_title("Pressure vs Maximum Burnup", fontsize=12)
    axes[0, 1].grid(alpha=0.3)
    
    # Plot 3: Pressure vs Gap
    axes[0, 2].plot(gap, pressure, "o-", color="#55A868", lw=2, markersize=6)
    axes[0, 2].axvline(0, color="gray", linestyle="--", alpha=0.7)
    axes[0, 2].set_xlabel("Gap (μm)", fontsize=11)
    axes[0, 2].set_ylabel("Contact Pressure (MPa)", fontsize=11)
    axes[0, 2].set_title("Pressure vs Gap", fontsize=12)
    axes[0, 2].grid(alpha=0.3)
    
    # Plot 4: Pressure vs Maximum Temperature
    axes[1, 0].plot(t_max, pressure, "^-", color="#C44E52", lw=2, markersize=6)
    axes[1, 0].set_xlabel("Maximum Temperature (K)", fontsize=11)
    axes[1, 0].set_ylabel("Contact Pressure (MPa)", fontsize=11)
    axes[1, 0].set_title("Pressure vs Maximum Temperature", fontsize=12)
    axes[1, 0].grid(alpha=0.3)
    
    # Plot 5: Pressure vs Minimum Temperature
    axes[1, 1].plot(t_min, pressure, "v-", color="#911B1F", lw=2, markersize=6)
    axes[1, 1].set_xlabel("Minimum Temperature (K)", fontsize=11)
    axes[1, 1].set_ylabel("Contact Pressure (MPa)", fontsize=11)
    axes[1, 1].set_title("Pressure vs Minimum Temperature", fontsize=12)
    axes[1, 1].grid(alpha=0.3)
    
    # Plot 6: Pressure vs Temperature Difference
    temp_diff = t_max - t_min
    axes[1, 2].plot(temp_diff, pressure, "d-", color="#8172B2", lw=2, markersize=6)
    axes[1, 2].set_xlabel("Temperature Difference (K)", fontsize=11)
    axes[1, 2].set_ylabel("Contact Pressure (MPa)", fontsize=11)
    axes[1, 2].set_title("Pressure vs Temperature Difference", fontsize=12)
    axes[1, 2].grid(alpha=0.3)
    
    fig.tight_layout()
    output_path = os.path.join(OUT, "pressure_vs_all_variables.png")
    fig.savefig(output_path, dpi=150)
    plt.close(fig)
    print(f"Created {output_path}")
    return output_path

def main():
    """Main function to generate all pressure analysis plots."""
    print("[pressure_analysis] Loading data from history.csv...")
    data = load_history_data()
    
    print(f"[pressure_analysis] Found {len(data)} data points")
    print(f"[pressure_analysis] Time range: {data['time_days'][0]:.1f} to {data['time_days'][-1]:.1f} days")
    print(f"[pressure_analysis] Pressure range: {data['contact_pressure_MPa'].min():.2f} to {data['contact_pressure_MPa'].max():.2f} MPa")
    
    print("[pressure_analysis] Creating pressure vs time plot...")
    plot_pressure_vs_time(data)
    
    print("[pressure_analysis] Creating pressure vs all variables plot...")
    plot_pressure_vs_all_variables(data)
    
    print("[pressure_analysis] All plots generated successfully!")

if __name__ == "__main__":
    main()