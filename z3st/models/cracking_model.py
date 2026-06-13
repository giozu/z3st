# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
Isotropic softening model for cracking due to thermal stresses in oxide materials (Barani et al., NED 342, 2019).

In oxide fuel, the temperature gradient cracks the pellet as
soon as the thermal tensile stress exceeds the fracture stress.
Resolving the crack pattern is neither affordable in an
engineering rod; instead, the cracked material is represented as an
ISOTROPICALLY SOFTENED solid: the elastic constants are rescaled as a function
of the number of macroscopic cracks n, conserving principal strains and
minimising the squared deviation of the principal stresses between the cracked
(anisotropic) and the equivalent isotropic description.

Scaling of the elastic constants, applied from the VIRGIN
constants (E, nu):

    f(nu)      = (2/3) * (2 - nu)/(2 + nu) * 1/(1 - nu)
    E_iso(n)   = f(nu)^n * E
    nu_iso(n)  = nu / (2^n + (2^n - 1) * nu)

Number of cracks, an empirical correlation on the rod-average:

    n = 0                                                    LHR <  LHR0
    n = n0 + (n_inf - n0) * (1 - exp(-(LHR - LHR0)/tau))     LHR >= LHR0

with LHR0 = 5 kW/m (first crack, n0 = 1), n_inf = 12, tau = 21 kW/m. Crack
healing is not modelled, so n is driven by the MAXIMUM rod-average LHR seen so
far in the power history (irreversible).

Card (material):

    cracking: isotropic       # opt-in
    cracking_lhr0:  5.0e3     # (W/m)  optional, default 5 kW/m
    cracking_n0:    1.0       # (-)    optional
    cracking_n_inf: 12.0      # (-)    optional
    cracking_tau:   21.0e3    # (W/m)  optional, default 21 kW/m

The rescale is applied once per time step (spine.parameters, i.e. before the
solve); the mechanical weak form is rebuilt every staggered iteration, so the
updated lmbda/G are picked up with no further plumbing. The derived constants
(lmbda, G, bulk_modulus) are recomputed from (E_iso, nu_iso); the virgin
values are kept in the card as E_virgin / nu_virgin.

"""

import math


class CrackingModel:

    def cracking_active(self, material):
        return str(material.get("cracking", "")).lower() == "isotropic"

    @staticmethod
    def _n_cracks(lhr, lhr0, n0, n_inf, tau):
        """Number of cracks at rod-average LHR (W/m)."""
        if lhr < lhr0:
            return 0.0
        return n0 + (n_inf - n0) * (1.0 - math.exp(-(lhr - lhr0) / tau))

    def update_cracking(self):
        """Per-step rescale of the elastic constants of every cracking
        material from the maximum rod-average LHR seen so far. Called from
        spine.parameters once per time step."""
        lhr = float(getattr(self, "lhr", 0.0))
        for name, mat in self.materials.items():
            if not self.cracking_active(mat):
                continue

            # Virgin constants, captured on first call
            if "E_virgin" not in mat:
                mat["E_virgin"] = float(mat["E"])
                mat["nu_virgin"] = float(mat["nu"])
                mat["_lhr_max"] = 0.0

            mat["_lhr_max"] = max(float(mat["_lhr_max"]), lhr)

            lhr0 = float(mat.get("cracking_lhr0", 5.0e3))
            n0 = float(mat.get("cracking_n0", 1.0))
            n_inf = float(mat.get("cracking_n_inf", 12.0))
            tau = float(mat.get("cracking_tau", 21.0e3))

            n = self._n_cracks(mat["_lhr_max"], lhr0, n0, n_inf, tau)

            E0, nu0 = mat["E_virgin"], mat["nu_virgin"]
            if n <= 0.0:
                E_iso, nu_iso = E0, nu0
            else:
                f = (2.0 / 3.0) * (2.0 - nu0) / (2.0 + nu0) / (1.0 - nu0)
                E_iso = (f ** n) * E0
                two_n = 2.0 ** n
                nu_iso = nu0 / (two_n + (two_n - 1.0) * nu0)

            mat["E"] = E_iso
            mat["nu"] = nu_iso
            mat["bulk_modulus"] = E_iso / (3 * (1 - 2 * nu_iso))

            lmbda_new = E_iso * nu_iso / ((1 + nu_iso) * (1 - 2 * nu_iso))
            G_new = E_iso / (2 * (1 + nu_iso))
            if hasattr(mat["lmbda"], "value"):
                mat["lmbda"].value = lmbda_new
                mat["G"].value = G_new
            else:
                mat["lmbda"] = lmbda_new
                mat["G"] = G_new

            print(
                f"  [cracking] {name}: LHR_max = {mat['_lhr_max']/1e3:.1f} kW/m → "
                f"n = {n:.2f} cracks, E_iso/E = {E_iso/E0:.4f}, nu_iso = {nu_iso:.4f}"
            )
