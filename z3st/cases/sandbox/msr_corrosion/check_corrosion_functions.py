#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Checks of the corrosion model's pure functions against Bhave et al. 2023 values."""

from math import isclose

from z3st.models.corrosion_model import (
    atomic_fraction_cr,
    eq50_mass_loss,
    hart_diffusivity,
    mass_fraction_cr,
    profile_mass_loss,
    redox_surface_fraction,
)

cr = dict(D0_bulk=4.28e-4, Q_bulk=2.852e5, D0_gb=4.85e-2, Q_gb=1.962e5, gb_width=5e-10)
assert isclose(hart_diffusivity({**cr, "grain_size": 16.22e-6}, 973.0), 4.401e-17, rel_tol=3e-3)  # reference used R = 8.314
assert isclose(hart_diffusivity({**cr, "grain_size": 4.64e-6}, 973.0), 1.533e-16, rel_tol=3e-3)
assert isclose(mass_fraction_cr(0.0561), 0.0500, abs_tol=1e-4)
assert isclose(mass_fraction_cr(0.2201), 0.2000, abs_tol=1e-4)
assert isclose(atomic_fraction_cr(1.8e-3), 2.0315e-3, rel_tol=1e-3)
# Authors' 1D inputs (ref. [64]): salt potentials w_Ni = -0.45065, w_Cr = -0.9699 eV give
# x_Cr = 1.7726e-3 at the metal surface; the residual is rounding of kB, F and H_xs
xs, xn = redox_surface_fraction(973.0, E_F=3.3607, E0_F=2.871)
assert isclose(xs, 1.7726e-3, rel_tol=2e-3) and isclose(xn, 0.9982, abs_tol=1e-4)
# Bhave's input constants (Q/R = 34311.274 K and 23576.457 K, d = 16.1122 um) give Table 2's D
bh = dict(D0_bulk=4.278214e-4, Q_bulk=2.852798e5, D0_gb=4.849144e-2, Q_gb=1.960256e5, gb_width=5e-10)
assert isclose(hart_diffusivity({**bh, "grain_size": 16.1122e-6}, 973.0), 4.531e-17, rel_tol=5e-3)
assert isclose(hart_diffusivity({**bh, "grain_size": 4.64e-6}, 973.0), 1.569e-16, rel_tol=5e-3)
# Reference values for V1 and V2, kg/m2 at 1000 h
assert isclose(profile_mass_loss(3.6e6, 4.401e-17, 0.0561, 2.0315e-3, 8700.0), 5.9735e-3, rel_tol=2e-3)
assert isclose(profile_mass_loss(3.6e6, 1e-16, 0.2, 0.0, 8570.0), 3.3587e-2, rel_tol=2e-3)
assert isclose(eq50_mass_loss(0.0500, 1.8e-3, 8700.0, 4.401e-17, 3.6e6), 5.958e-3, rel_tol=2e-3)

# A Cr-saturated surface must be rejected
try:
    redox_surface_fraction(973.0, E_F=3.36, E0_F=2.871, x_Cr_salt=0.1)
except ValueError:
    pass
else:
    raise AssertionError("redox_surface_fraction accepted a Cr-saturated surface")

print("check_corrosion_functions: all passed")
