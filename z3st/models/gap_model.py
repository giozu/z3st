# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.2.0 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---


import dolfinx
import numpy as np
from mpi4py import MPI
from petsc4py import PETSc


class GapModel:
    def __init__(self):
        print("[GapModel] initializer")

    def set_gap_conductance(self, T_i):

        # Gap condutance
        if self.gap_model == "Fixed":
            h_gap_value = self.h_gap_value
            print(f"\n[INFO] Applying fixed gap conductance (h_gap = {h_gap_value:.2f} W/m²K)")

        elif self.gap_model == "Gas":

            self.set_gap_temperature(T_i)

            gas_thermal_conductivity = (
                self.h_gap_value * 1e-4 * self.gap_temperature**0.79
            )  # (W/m-K)

            # Deformed gap
            live_gap = getattr(self, "_last_gap", None)
            if live_gap is not None:
                floor = float(getattr(self, "gap_contact_thickness", 1.0e-6))
                gap_size = max(float(live_gap), floor)
                gap_label = "live, deformed"
            else:
                name_a = getattr(self, "gap_surface_a", "lateral_1")
                name_b = getattr(self, "gap_surface_b", "inner_2")
                missing = [n for n in (name_a, name_b) if n not in self.label_map]
                if missing:
                    raise KeyError(
                        f"Gas gap conductance: surface label(s) {missing} not in "
                        "the geometry label map — set models.gap_conductance."
                        "surface_a / surface_b to the paired facet labels."
                    )
                gap_size = self.average_gap_distance(
                    self.mesh,
                    self.facet_tags,
                    label_a=self.label_map[name_a],
                    label_b=self.label_map[name_b],
                )
                gap_label = "cold, undeformed"

            h_gap_value = gas_thermal_conductivity / gap_size  # (W/m2-K)

            print(f"  → k_gas       = {gas_thermal_conductivity:.2f} W/m·K")
            print(f"  → gap_size    = {gap_size*1e3:.3f} mm ({gap_label})")
            print(f"  → h_gap       = {h_gap_value:.2f} W/m²K")

        else:
            h_gap_value = 0.0

        # Contact-coupled conductance (Todreas & Kazimi, Nuclear Systems I,
        # 3rd ed., Eq. 8.141/8.142): on gap closure the pellet-clad contact
        # pressure adds a contact term to the open-gap value.
        if self.on.get("contact", False) and getattr(self, "gap_contact_coupling", False):
            h_contact = self.contact_conductance()
            if h_contact > 0.0:
                print(
                    f"  → h_open      = {h_gap_value:.2f} W/m²K, "
                    f"h_contact = {h_contact:.2f} W/m²K (gap closed)"
                )
                h_gap_value += h_contact

        # Under-relax h between staggered iterations (models.gap_conductance.
        # relax, default 1.0 = off): damps the contact-pressure ↔ conductance
        # ↔ temperature feedback at its source. _h_gap_prev is reset at every
        # time step (solve_staggered), so the damping never lags across steps.
        omega_h = float(getattr(self, "gap_relax", 1.0))
        h_prev = getattr(self, "_h_gap_prev", None)
        if omega_h < 1.0 and h_prev is not None:
            h_gap_value = omega_h * h_gap_value + (1.0 - omega_h) * h_prev
            print(f"  → h_gap (relaxed, ω={omega_h:.2f}) = {h_gap_value:.2f} W/m²K")
        self._h_gap_prev = h_gap_value

        # Persistent Constant so the cached thermal form sees the update in
        # place (no form rebuild needed between staggered iterations).
        if not hasattr(self, "_h_gap_const"):
            self._h_gap_const = dolfinx.fem.Constant(self.mesh, PETSc.ScalarType(h_gap_value))
        else:
            self._h_gap_const.value = PETSc.ScalarType(h_gap_value)

        return self._h_gap_const

    def contact_conductance(self):
        """
        Solid-contact gap conductance on gap closure, Todreas & Kazimi,
        *Nuclear Systems I*, 3rd ed., Eq. 8.141 (Ross-Stoute form):

            h_contact = C * (2 k_f k_c)/(k_f + k_c) * P_i / (H * sqrt(delta_g)),

        with the pellet-clad contact pressure P_i taken from the penalty
        :class:`ContactModel`. The empirical constant C = 10 ft^-1/2 is here
        expressed in SI as C_SI = 18.11 m^-1/2 so that, with k in W/(m.K),
        delta_g in m and the dimensionless ratio P_i/H, the result is W/(m^2.K).
        Returns 0 when the surfaces are not in contact (P_i = 0).
        """
        P_i = float(getattr(self, "_last_pressure", 0.0))  # Pa, from ContactModel
        if P_i <= 0.0:
            return 0.0

        # harmonic mean of the two solid conductivities; a symbolic k(T) card
        # is evaluated at the current gap temperature (UFL folds constants, so
        # k_func(float) returns a plain number). Prefer the two materials of
        # the gap pair itself (recorded by set_gap_temperature) — a rod model
        # can carry other conducting cards (spring, coolant, structure) whose
        # dict position must not decide the contact conductance.
        pair = getattr(self, "_gap_pair_labels", None)
        mats = ([self.materials[n] for n in pair if n in self.materials]
                if pair is not None else [])
        if len(mats) < 2:
            mats = list(self.materials.values())
        ks = [k for k in (self._k_at_gap(m) for m in mats)
              if k is not None]
        if len(ks) < 2:
            print(
                "  [WARNING] contact_coupling enabled but fewer than two materials "
                "carry a usable 'k' (numeric card or symbolic k(T)); "
                "h_contact = 0 — gap conductance stays at the open-gap value."
            )
            return 0.0
        k_f, k_c = ks[0], ks[1]
        k_harm = 2.0 * k_f * k_c / (k_f + k_c)

        H = self.gap_meyer_hardness            # Pa (Meyer hardness, softer solid)
        delta_g = self.gap_contact_thickness   # m (roughness-based gas space)
        C_SI = 18.11                           # m^-1/2  (= 10 ft^-1/2, Eq. 8.141)

        return C_SI * k_harm * (P_i / H) / (delta_g ** 0.5)

    def _k_at_gap(self, material):
        """Numeric conductivity of one material for the Ross-Stoute harmonic
        mean: a plain card value, or a symbolic k(T) card evaluated at the
        current mean gap temperature (falls back to the card's T_initial when
        the gap temperature is not yet available, e.g. Fixed gap model)."""
        k_card = material.get("k")
        if isinstance(k_card, (int, float)):
            return float(k_card)
        k_func = material.get("_k_func")
        if k_func is not None:
            T_gap = float(
                getattr(self, "gap_temperature",
                        material.get("T_initial", material.get("T_ref", 293.15)))
            )
            try:
                return float(k_func(T_gap))
            except (TypeError, ValueError):
                return None
        return None

    def set_gap_temperature(self, T_i):
        """Mean gap temperature from the first Robin BC that carries a
        ``pair:`` entry (convective Robin entries — h_conv/T_ext — are skipped,
        and every material's BC list is scanned until a pair is found).

        MPI-safe: each surface mean is a global sum over *owned* dofs divided
        by the global count, so all ranks agree on ``gap_temperature`` even
        when the paired surfaces are split across partitions.
        """
        comm = self.mesh.comm
        n_owned = self.V_t.dofmap.index_map.size_local

        def global_mean(dofs):
            dofs = np.asarray(dofs)
            owned = dofs[dofs < n_owned] if dofs.size else dofs
            s = comm.allreduce(float(T_i.x.array[owned].sum()) if owned.size else 0.0,
                               op=MPI.SUM)
            n = comm.allreduce(int(owned.size), op=MPI.SUM)
            return s / n if n > 0 else None

        for label in self.materials:
            for bc_info in self.robin_thermal.get(label, []):
                if "pair" not in bc_info:
                    continue  # convective Robin BC, not a gap coupling
                region_id = bc_info["id"]
                pair_region = bc_info["pair"]

                mean_here = global_mean(self.mgr.locate_facets_dofs(region_id, self.V_t))
                mean_other = global_mean(
                    self.mgr.locate_facets_dofs(self.label_map[pair_region], self.V_t)
                )
                if mean_here is None or mean_other is None:
                    print("  [WARNING] gap pair surface has no dofs anywhere; "
                          "gap temperature not updated.")
                    return

                self.gap_temperature = 0.5 * (mean_here + mean_other)
                # Remember the paired regions so the Ross-Stoute contact term
                # uses the conductivities of the actual contacting materials.
                self._gap_pair_labels = (label, pair_region)

                print(
                    f"  → Average gap temperature between {label} and {pair_region}: {self.gap_temperature:.2f} K"
                )
                return

        print("  [WARNING] Gas gap model active but no Robin BC with a 'pair:' "
              "entry was found; gap temperature not updated.")

    def average_gap_distance(self, mesh, ft, label_a, label_b):

        facets_a = ft.find(label_a)  # e.g., facets of cyl_1_interface
        facets_b = ft.find(label_b)  # e.g., facets of cyl_2_interface

        x = mesh.geometry.x
        topology = mesh.topology
        fdim = topology.dim - 1

        topology.create_connectivity(fdim, 0)  # facet → vertex
        # dolfinx computes entity midpoints in C++; the Python loop this replaces
        # did the same thing (mean of the facet's vertex coordinates) one facet at
        # a time.
        centroids_a = dolfinx.mesh.compute_midpoints(mesh, fdim, facets_a)
        centroids_b = dolfinx.mesh.compute_midpoints(mesh, fdim, facets_b)

        # MPI: gather both surfaces globally — with the mesh partitioned, a rank
        # can own facets of only one side (empty-tree crash) or miss the true
        # nearest neighbour across a rank boundary. Gap surfaces are small, so
        # allgathering their centroids is cheap and makes the result identical
        # on every rank.
        comm = mesh.comm
        if comm.size > 1:
            gdim = x.shape[1]
            centroids_a = np.concatenate(
                [c.reshape(-1, gdim) for c in comm.allgather(centroids_a)]
            )
            centroids_b = np.concatenate(
                [c.reshape(-1, gdim) for c in comm.allgather(centroids_b)]
            )

        if centroids_a.size == 0 or centroids_b.size == 0:
            raise RuntimeError(
                f"average_gap_distance: no facets found for gap labels "
                f"{label_a}/{label_b} — check the geometry label map."
            )

        from scipy.spatial import cKDTree

        tree_b = cKDTree(centroids_b)
        distances, _ = tree_b.query(centroids_a, k=1)

        return distances.mean()
