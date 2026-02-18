#!/usr/bin/env python3
"""Generate reference SOAP power spectrum values for the HCN test system.

This script computes SOAP descriptors using featomic (Python) and writes
reference values that the C++ SoapDescriptorEngine test can compare against.

Usage:
    pixi run -e cuda python gen_soap_reference.py

The hyperparameters here must exactly match those used in the C++ test
(FeatomicTest.cpp).
"""

import json

import numpy as np
import torch
from featomic.torch import SoapPowerSpectrum

# --- HCN system (must match pos.con) ---
# Positions in Angstroms (from the .con file)
positions = torch.tensor(
    [
        [12.49920530576333633, 12.49960080579726807, 12.54209698318736699],  # C
        [12.50069048683807615, 12.50039919420273193, 11.38365673946223389],  # N
        [12.50079469423666367, 12.50006476947390155, 13.61634326053776611],  # H
    ],
    dtype=torch.float64,
    requires_grad=True,
)

# Atomic numbers: C=6, N=7, H=1
types = torch.tensor([6, 7, 1], dtype=torch.int32)

# Cell: 25 Ang cubic
cell = torch.tensor(
    [[25.0, 0.0, 0.0], [0.0, 25.0, 0.0], [0.0, 0.0, 25.0]], dtype=torch.float64
)

# PBC: match C++ matterToTensors() which derives PBC from cell norms.
# A 25 Ang vacuum box has non-zero cell vectors -> PBC=[True, True, True].
pbc = torch.tensor([True, True, True])

# --- SOAP hyperparameters (must match C++ test) ---
CUTOFF = 6.0
SMOOTHING_WIDTH = 0.5
DENSITY_WIDTH = 0.3
MAX_ANGULAR = 4
MAX_RADIAL = 6

hypers = {
    "cutoff": {
        "radius": CUTOFF,
        "smoothing": {"type": "ShiftedCosine", "width": SMOOTHING_WIDTH},
    },
    "density": {"type": "Gaussian", "width": DENSITY_WIDTH},
    "basis": {
        "type": "TensorProduct",
        "max_angular": MAX_ANGULAR,
        "radial": {"type": "Gto", "max_radial": MAX_RADIAL},
    },
}

calculator = SoapPowerSpectrum(**hypers)

# --- Build a metatomic System from tensors ---
from metatomic.torch import System

system = System(types=types, positions=positions, cell=cell, pbc=pbc)

# --- Compute descriptor with position gradients ---
descriptor = calculator.compute(
    [system],
    gradients=["positions"],
)

# Register autograd so we can backprop through positions
from featomic.torch import register_autograd

descriptor = register_autograd([system], descriptor, forward_gradients=["positions"])

# --- Average over atoms (same logic as C++ SoapDescriptorEngine) ---
keys = descriptor.keys
n_blocks = len(keys)

block_means = []
for b in range(n_blocks):
    block = descriptor.block(b)
    values = block.values  # [n_atoms_in_block, n_components, n_props]
    n_samples = values.shape[0]
    flat = values.reshape(n_samples, -1)  # [n_atoms, D_block]
    block_sum = flat.sum(dim=0)  # [D_block]
    block_means.append(block_sum)

avg_descriptor = torch.cat(block_means, dim=0)  # [D_total]
n_atoms = positions.shape[0]
avg_descriptor = avg_descriptor / float(n_atoms)

D = avg_descriptor.shape[0]
print(f"Descriptor dimension D = {D}")
print(f"Descriptor norm = {avg_descriptor.norm().item():.15e}")

# --- Compute Jacobian [D, 3N] via autograd ---
n_cart = 3 * n_atoms
jacobian = torch.zeros(D, n_cart, dtype=torch.float64)

for d in range(D):
    if positions.grad is not None:
        positions.grad.zero_()
    avg_descriptor[d].backward(retain_graph=True)
    if positions.grad is not None:
        jacobian[d] = positions.grad.reshape(-1).clone()

print(f"Jacobian shape = {jacobian.shape}")
print(f"Jacobian Frobenius norm = {jacobian.norm().item():.15e}")

# --- Save reference data ---
ref = {
    "hypers": {
        "cutoff": CUTOFF,
        "smoothing_width": SMOOTHING_WIDTH,
        "density_width": DENSITY_WIDTH,
        "max_angular": MAX_ANGULAR,
        "max_radial": MAX_RADIAL,
    },
    "system": {
        "n_atoms": int(n_atoms),
        "types": types.tolist(),
        "positions": positions.detach().tolist(),
    },
    "descriptor_dim": D,
    "descriptor_norm": avg_descriptor.norm().item(),
    "descriptor": avg_descriptor.detach().tolist(),
    "jacobian_frobenius_norm": jacobian.norm().item(),
    # Store first few rows of Jacobian for spot-checking
    "jacobian_first_5_rows": jacobian[:5].tolist(),
}

# Find first 5 non-zero Jacobian rows for meaningful spot-checks
row_norms = jacobian.norm(dim=1)
nz_mask = row_norms > 1e-15
nz_indices = torch.where(nz_mask)[0].tolist()[:5]
ref["jacobian_nonzero_rows"] = {
    "indices": nz_indices,
    "values": [jacobian[i].tolist() for i in nz_indices],
    "norms": [row_norms[i].item() for i in nz_indices],
}

import pathlib

out_dir = pathlib.Path(__file__).resolve().parent

with open(out_dir / "soap_reference.json", "w") as f:
    json.dump(ref, f, indent=2)

print("Wrote soap_reference.json")

# --- Also write a compact binary for exact comparison ---
np.savez(
    out_dir / "soap_reference.npz",
    descriptor=avg_descriptor.detach().numpy(),
    jacobian=jacobian.detach().numpy(),
)
print("Wrote soap_reference.npz")
