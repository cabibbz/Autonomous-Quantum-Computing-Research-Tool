"""
Experiment 2a: CHSH Violation Landscape
Use statevector math (no sampling) to compute exact CHSH S across angle grid.
2 qubits, pure linear algebra — should be very fast.
"""

from qiskit.circuit import QuantumCircuit
from qiskit.quantum_info import Statevector
import numpy as np
import json, os, time

RESULTS_DIR = "results"
os.makedirs(RESULTS_DIR, exist_ok=True)

t0 = time.time()

# For Phi+ state with RY measurement rotations, correlation is:
# E(a, b) = cos(a - b)
# This is exact — no sampling needed for the landscape.

def chsh_exact(a0, a1, b0, b1):
    """Exact CHSH S for Phi+ state."""
    E00 = np.cos(a0 - b0)
    E01 = np.cos(a0 - b1)
    E10 = np.cos(a1 - b0)
    E11 = np.cos(a1 - b1)
    return E00 - E01 + E10 + E11

# Sweep: fix a0=0, vary a1 and b0 (b1 = b0 + pi/2)
n_points = 50  # Can do a fine grid since this is just math
a1_values = np.linspace(0, 2*np.pi, n_points)
b0_values = np.linspace(0, 2*np.pi, n_points)
S_grid = np.zeros((n_points, n_points))

for i, a1 in enumerate(a1_values):
    for j, b0 in enumerate(b0_values):
        S_grid[i, j] = chsh_exact(0, a1, b0, b0 + np.pi/2)

# Also do a 1D sweep: fix optimal structure, vary one parameter
# a0=0, a1=pi/2, sweep b0 (b1 = b0+pi/2)
b_sweep = np.linspace(0, 2*np.pi, 200)
S_sweep = [chsh_exact(0, np.pi/2, b, b + np.pi/2) for b in b_sweep]

# Verify: also compute via statevector for a few points to confirm analytic formula
verification = []
for a1, b0 in [(np.pi/2, np.pi/4), (np.pi/3, np.pi/6), (np.pi, np.pi/2)]:
    # Numeric check via statevector
    E_numeric = {}
    for ai, a in enumerate([0, a1]):
        for bi, b in enumerate([b0, b0 + np.pi/2]):
            qc = QuantumCircuit(2)
            qc.h(0)
            qc.cx(0, 1)
            qc.ry(-a, 0)
            qc.ry(-b, 1)
            sv = Statevector.from_instruction(qc)
            probs = sv.probabilities_dict()
            E = sum((1-2*int(k[1]))*(1-2*int(k[0]))*p for k, p in probs.items())
            E_numeric[f"a{ai}_b{bi}"] = E
    S_numeric = E_numeric["a0_b0"] - E_numeric["a0_b1"] + E_numeric["a1_b0"] + E_numeric["a1_b1"]
    S_analytic = chsh_exact(0, a1, b0, b0 + np.pi/2)
    verification.append({
        "a1": float(a1), "b0": float(b0),
        "S_analytic": float(S_analytic),
        "S_numeric": float(S_numeric),
        "match": bool(abs(S_analytic - S_numeric) < 1e-10)
    })

dt = time.time() - t0

# Analysis
max_idx = np.unravel_index(np.argmax(S_grid), S_grid.shape)
min_idx = np.unravel_index(np.argmin(S_grid), S_grid.shape)
violation_mask = np.abs(S_grid) > 2.0

results = {
    "experiment": "chsh_angle_landscape",
    "method": "exact_analytic (E(a,b)=cos(a-b) for Phi+)",
    "grid_size": n_points,
    "a0_fixed": 0.0,
    "a1_range": [0, float(2*np.pi)],
    "b0_range": [0, float(2*np.pi)],
    "b1_rule": "b0 + pi/2",
    "S_grid": S_grid.tolist(),
    "max_S": float(S_grid.max()),
    "max_at": {"a1": float(a1_values[max_idx[0]]), "b0": float(b0_values[max_idx[1]])},
    "min_S": float(S_grid.min()),
    "min_at": {"a1": float(a1_values[min_idx[0]]), "b0": float(b0_values[min_idx[1]])},
    "violation_fraction": float(np.mean(violation_mask)),
    "S_statistics": {
        "mean": float(np.mean(S_grid)),
        "std": float(np.std(S_grid)),
    },
    "b_sweep_optimal": {
        "b_values": b_sweep.tolist(),
        "S_values": S_sweep,
        "description": "S vs b0 with a0=0, a1=pi/2, b1=b0+pi/2"
    },
    "verification": verification,
    "runtime_seconds": dt,
}

with open(os.path.join(RESULTS_DIR, "sprint_002a_chsh_landscape.json"), "w") as f:
    json.dump(results, f, indent=2)

print(f"Runtime: {dt:.2f}s")
print(f"Grid: {n_points}x{n_points}")
print(f"S range: [{S_grid.min():.4f}, {S_grid.max():.4f}]")
print(f"Max S = {S_grid.max():.4f} at a1={a1_values[max_idx[0]]:.3f}, b0={b0_values[max_idx[1]]:.3f}")
print(f"Min S = {S_grid.min():.4f} at a1={a1_values[min_idx[0]]:.3f}, b0={b0_values[min_idx[1]]:.3f}")
print(f"Fraction of grid violating |S|>2: {np.mean(violation_mask):.1%}")
print(f"\nVerification (analytic vs statevector):")
for v in verification:
    print(f"  a1={v['a1']:.3f}, b0={v['b0']:.3f}: analytic={v['S_analytic']:.6f}, numeric={v['S_numeric']:.6f}, match={v['match']}")
print(f"\nSaved to results/sprint_002a_chsh_landscape.json")
