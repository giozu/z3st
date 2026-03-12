#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: pressure_vessel_2D

non-regression script
---------------------
Extracts stress fields for a 2D axisymmetric pressure vessel.
Geometry:
- Cylindrical shell (Ri=1.52, t=0.052, L=2.0)
- Hemi-spherical head (Ri=1.52, t=0.026)

Extractions:
1. Mid-height of the cylindrical shell (y = -1.0 m)
2. Middle of the hemi-spherical head (45 degrees, y = x)

"""

import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
import pandas as pd

# --.. ..- .-.. .-.. --- configuration --.. ..- .-.. .-.. ---
CASE_DIR = os.path.dirname(os.path.abspath(__file__))
VTU_FILE = os.path.join(CASE_DIR, "output", "fields.vtu")
OUT_DIR = os.path.join(CASE_DIR, "output")
os.makedirs(OUT_DIR, exist_ok=True)

# Geometry
Ri = 1.52  # m
t_cyl = 0.052  # m
t_head = 0.026  # m
L_cyl = 2.0  # m
Ro_cyl = Ri + t_cyl
Ro_head = Ri + t_head

# Parameters
TOLERANCE = 1.0e-3  # m

def extract_mid_cylinder(grid):
    """
    Extract stress at mid-height of cylinder (y = -L_cyl / 2).
    """
    print(f"\n[INFO] Extracting stress at cylinder mid-height (y = {-L_cyl/2})")
    
    # Target Y
    y_target = -L_cyl / 2.0
    
    # Get Cell Centers and Data
    # Prefer cell data for stress usually, but let's check what we have
    if "Stress_steel (cells)" in grid.cell_data:
        data_source = "cell"
        coords = grid.cell_centers().points
        stress_flat = grid.cell_data["Stress_steel (cells)"]
    elif "Stress_steel (points)" in grid.point_data:
        data_source = "point"
        coords = grid.points
        stress_flat = grid.point_data["Stress_steel (points)"]
    else:
        print("[ERROR] 'Stress_steel' not found in fields.")
        return

    # Filter by Y position
    mask = np.abs(coords[:, 1] - y_target) < TOLERANCE
    
    if np.sum(mask) == 0:
        print(f"[WARNING] No data found at y ≈ {y_target}")
        return

    x_vals = coords[mask, 0] # Radial coordinate
    stress_vals = stress_flat[mask]
    
    # Reshape stress to Nx3x3 if needed, or if it's already flat 9-component
    if stress_vals.ndim == 1:
        # Assuming 9 components
        # Verify shape
        pass 
    elif stress_vals.shape[1] == 9:
        pass
    else:
        print(f"[ERROR] Unexpected stress shape: {stress_vals.shape}")
        return

    # Sort by radius (x)
    sort_idx = np.argsort(x_vals)
    r = x_vals[sort_idx]
    s = stress_vals[sort_idx]
    
    # Components (Assuming tensor stored as xx, xy, xz, yx, yy, yz, zx, zy, zz)
    # Check PyVista/VTK ordering. Usually: XX, YY, ZZ, XY, YZ, XZ or similar.
    # Check shape first. If it is 9 comps, usually it is full tensor.
    # In 2D Axisymmetric (X=r, Y=z/axial):
    # s_xx -> Radial
    # s_yy -> Axial
    # s_zz -> Hoop (check if generalized plane strain or axi formulation usage)
    # s_xy -> Shear
    
    # Let's assume standard VTK symmetric tensor order or full 9
    # For now, let's just plot components 0 (XX), 4 (YY), 8 (ZZ) and maybe 1 (XY)
    
    sigma_rr = s[:, 0]
    sigma_yy = s[:, 4]
    sigma_zz = s[:, 8] # Hoop
    
    # Save CSV
    df = pd.DataFrame({
        "r": r,
        "sigma_rr": sigma_rr,
        "sigma_axial": sigma_yy,
        "sigma_hoop": sigma_zz
    })
    csv_path = os.path.join(OUT_DIR, "stress_cylinder_mid.csv")
    df.to_csv(csv_path, index=False)
    print(f"[INFO] Saved CSV: {csv_path}")
    
    # Plot
    plt.figure(figsize=(8, 6))
    plt.plot(r, sigma_rr * 1e-6, '-o', label=r'$\sigma_{rr}$ (Radial)')
    plt.plot(r, sigma_yy * 1e-6, '-s', label=r'$\sigma_{yy}$ (Axial)')
    plt.plot(r, sigma_zz * 1e-6, '-^', label=r'$\sigma_{\theta\theta}$ (Hoop)')
    
    plt.xlabel("Radius [m]")
    plt.ylabel("Stress [MPa]")
    plt.title(f"Stress at Cylinder Mid-Height (y={y_target} m)")
    plt.grid(True, linestyle="--", alpha=0.7)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "stress_cylinder_mid.png"), dpi=300)
    plt.close()


def extract_head_middle(grid):
    """
    Extract stress at middle of head (45 degrees: y = x, x > 0).
    """
    print(f"\n[INFO] Extracting stress at Head Middle (45 deg)")
    
    # Get Data
    if "Stress_steel (cells)" in grid.cell_data:
        coords = grid.cell_centers().points
        stress_flat = grid.cell_data["Stress_steel (cells)"]
    elif "Stress_steel (points)" in grid.point_data:
        coords = grid.points
        stress_flat = grid.point_data["Stress_steel (points)"]
    else:
        return

    x = coords[:, 0]
    y = coords[:, 1]
    
    # Filter for x > 0 and y > 0 (Head quadrant)
    # And close to y = x line
    
    # Distance from line x - y = 0 is |x - y| / sqrt(2)
    dist = np.abs(x - y) / np.sqrt(2)
    
    mask = (dist < TOLERANCE) & (x > 0.1) & (y > 0.1) # Avoid origin singularity if any
    
    if np.sum(mask) == 0:
        print(f"[WARNING] No data found at 45 deg line")
        return

    x_s = x[mask]
    y_s = y[mask]
    s_s = stress_flat[mask]
    
    # Spherical Radius coordinate from origin
    r_sph = np.sqrt(x_s**2 + y_s**2)
    
    # Sort
    sort_idx = np.argsort(r_sph)
    r_sph = r_sph[sort_idx]
    s_s = s_s[sort_idx]
    x_s = x_s[sort_idx]
    y_s = y_s[sort_idx]
    
    # Global Components
    sig_xx = s_s[:, 0] # Radial global
    sig_yy = s_s[:, 4] # Axial global
    sig_zz = s_s[:, 8] # Hoop global (also Hoop local)
    sig_xy = s_s[:, 1]
    
    # Transformation to Local Coordinates (Normal, Meridional)
    # Angle alpha = 45 deg = pi/4
    # Normal vector n = (cos a, sin a) = (1/sqrt(2), 1/sqrt(2))
    # Tangent vector t = (sin a, -cos a) = (1/sqrt(2), -1/sqrt(2)) (Wait, check direction)
    # Tangent to circle is perpendicular to radius. 
    # If radius is at 45 deg, tangent is at -45 deg (315 deg) or 135 deg.
    
    # Let's use rotation matrix for 45 deg.
    # x' (Normal direction, parallel to radius)
    # y' (Meridional direction, tangent to radius)
    
    theta = np.pi / 4
    c = np.cos(theta)
    s = np.sin(theta)
    
    # Sigma_normal (should be approximate radial stress through thickness)
    # n = (c, s)
    sig_normal = sig_xx * c**2 + sig_yy * s**2 + 2 * sig_xy * c * s
    
    # Sigma_meridional (tangential in plane)
    # t = (-s, c)
    sig_meridional = sig_xx * s**2 + sig_yy * c**2 - 2 * sig_xy * c * s
    
    # Save CSV
    df = pd.DataFrame({
        "r_spherical": r_sph,
        "x": x_s,
        "y": y_s,
        "sigma_normal": sig_normal,
        "sigma_meridional": sig_meridional,
        "sigma_hoop": sig_zz
    })
    csv_path = os.path.join(OUT_DIR, "stress_head_mid.csv")
    df.to_csv(csv_path, index=False)
    print(f"[INFO] Saved CSV: {csv_path}")
    
    # Plot
    plt.figure(figsize=(8, 6))
    plt.plot(r_sph, sig_normal * 1e-6, '-o', label=r'$\sigma_{normal}$ (Through-thickness)')
    plt.plot(r_sph, sig_meridional * 1e-6, '-s', label=r'$\sigma_{meridional}$ (Tangential)')
    plt.plot(r_sph, sig_zz * 1e-6, '-^', label=r'$\sigma_{hoop}$')
    
    plt.xlabel("Spherical Radius [m]")
    plt.ylabel("Stress [MPa]")
    plt.title("Stress at Hemi-Sphere Middle (45 deg)")
    plt.grid(True, linestyle="--", alpha=0.7)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "stress_head_mid.png"), dpi=300)
    plt.close()


def main():
    if not os.path.exists(VTU_FILE):
        print(f"[ERROR] VTU file not found: {VTU_FILE}")
        sys.exit(1)
        
    print(f"[INFO] Reading {VTU_FILE}...")
    grid = pv.read(VTU_FILE)
    
    extract_mid_cylinder(grid)
    extract_head_middle(grid)
    
    print("\n[INFO] Done.")

if __name__ == "__main__":
    main()
