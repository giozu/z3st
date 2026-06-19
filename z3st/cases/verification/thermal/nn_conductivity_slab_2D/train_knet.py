#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST: train k(T) = NN(T) --.. ..- .-.. .-.. ---
"""Offline training of the thermal-conductivity network.

Fits a small smooth MLP to k(T) (from nn_conductor.yaml::k_law) over the
temperature range spanned by the Dirichlet BCs, then saves a checkpoint
(knet.pt) in the format expected by z3st.models.nn_conductivity.NNConductivity.

Run once before the solve:  python3 train_knet.py
"""

import os

import numpy as np
import yaml
import torch

from z3st.models.nn_conductivity import _build_mlp

CASE = os.path.dirname(os.path.abspath(__file__))
HIDDEN = [32, 32]
ACTIVATION = "tanh"          # smooth -> continuous dk/dT

# --. read the reference law and the temperature range --..
with open(os.path.join(CASE, "nn_conductor.yaml")) as f:
    mat = yaml.safe_load(f)
a, b = float(mat["k_law"]["a"]), float(mat["k_law"]["b"])

with open(os.path.join(CASE, "boundary_conditions.yaml")) as f:
    bcs = yaml.safe_load(f)
temps = [float(e["temperature"]) for e in bcs["thermal"]["slab"] if e["type"] == "Dirichlet"]
T_min, T_max = min(temps), max(temps)

PAD = 0.20  # fraction of the BC span added on each side
span = T_max - T_min
T_lo, T_hi = T_min - PAD * span, T_max + PAD * span
T0, Tscale = 0.5 * (T_lo + T_hi), 0.5 * (T_hi - T_lo)
print(f"fitting k(T) = 1/({a} + {b}*T) over [{T_lo:.0f}, {T_hi:.0f}] K "
      f"(BC span [{T_min:.0f}, {T_max:.0f}] padded by {PAD:.0%})")

# --. train --..
torch.manual_seed(0)
net = _build_mlp(HIDDEN, ACTIVATION)
T_train = torch.linspace(T_lo, T_hi, 400).unsqueeze(-1)
# k_train = a + b * T_train   # k(T)
k_train = 1/(a + b * T_train) # k(T)
Tn = (T_train - T0) / Tscale
opt = torch.optim.Adam(net.parameters(), lr=1e-2)
for step in range(4000):
    opt.zero_grad()
    loss = torch.mean((net(Tn) - k_train) ** 2)
    loss.backward()
    opt.step()
print(f"final training MSE = {loss.item():.3e}")

# fit quality
net.double().eval()
with torch.no_grad():
    Tg = torch.linspace(T_lo, T_hi, 200, dtype=torch.float64).unsqueeze(-1)
    k_true = 1.0 / (a + b * Tg)
    rel = ((net((Tg - T0) / Tscale) - k_true).abs().max().item()) / k_true.max().item()
print(f"max relative fit error = {rel:.2e}")

# --. save checkpoint --..
out = os.path.join(CASE, "knet.pt")
torch.save(
    {"state_dict": net.state_dict(), "hidden": HIDDEN,
     "activation": ACTIVATION, "norm": [T0, Tscale]},
    out,
)
print(f"saved -> {out}")
