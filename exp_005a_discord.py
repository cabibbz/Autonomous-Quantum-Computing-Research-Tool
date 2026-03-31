"""
Experiment 5a: Quantum Discord — Separating Quantum from Classical Correlations

Compute pairwise quantum discord for GHZ, W, and Cluster states (6 qubits).
Discord = MI - Classical Correlations.
Discord > 0 means irreducibly quantum correlations.

For 2-qubit reduced state rho_AB:
  D(A|B) = S(rho_A) - S(rho_AB) + min_{meas on B} sum_k p_k S(rho_A|k)

The minimization is over projective measurements on B, parameterized by
(theta, phi) on the Bloch sphere.
"""

from qiskit.circuit import QuantumCircuit
from qiskit.quantum_info import Statevector, partial_trace, entropy
import numpy as np
from scipy.optimize import minimize
import json, os, time

RESULTS_DIR = "results"
os.makedirs(RESULTS_DIR, exist_ok=True)

t0 = time.time()
N = 6

def make_ghz(n):
    qc = QuantumCircuit(n)
    qc.h(0)
    for i in range(1, n):
        qc.cx(0, i)
    return Statevector.from_instruction(qc)

def make_w(n):
    qc = QuantumCircuit(n)
    qc.x(0)
    for i in range(n - 1):
        theta = 2 * np.arccos(np.sqrt(1 / (n - i)))
        qc.cry(theta, i, i + 1)
        qc.cx(i + 1, i)
    return Statevector.from_instruction(qc)

def make_cluster(n):
    qc = QuantumCircuit(n)
    for i in range(n):
        qc.h(i)
    for i in range(n - 1):
        qc.cz(i, i + 1)
    return Statevector.from_instruction(qc)

def von_neumann_entropy(rho_matrix):
    """Von Neumann entropy from a density matrix (numpy array)."""
    eigenvalues = np.linalg.eigvalsh(rho_matrix)
    eigenvalues = eigenvalues[eigenvalues > 1e-15]
    return -np.sum(eigenvalues * np.log2(eigenvalues))

def get_2qubit_rdm(sv, i, j, n_total):
    """Get the 2-qubit reduced density matrix for qubits i,j."""
    trace_out = [k for k in range(n_total) if k not in (i, j)]
    if not trace_out:
        return np.array(sv.to_operator())
    rho = partial_trace(sv, trace_out)
    return np.array(rho)

def measurement_projectors(theta, phi):
    """Projective measurement on a qubit along (theta, phi) on Bloch sphere.
    Returns two 2x2 projectors |n><n| and |n_perp><n_perp|."""
    cos_t2 = np.cos(theta / 2)
    sin_t2 = np.sin(theta / 2)
    # |n> = cos(t/2)|0> + e^{i*phi}*sin(t/2)|1>
    n_ket = np.array([cos_t2, np.exp(1j * phi) * sin_t2])
    # |n_perp> = sin(t/2)|0> - e^{i*phi}*cos(t/2)|1>
    n_perp = np.array([sin_t2, -np.exp(1j * phi) * cos_t2])
    P0 = np.outer(n_ket, n_ket.conj())
    P1 = np.outer(n_perp, n_perp.conj())
    return [P0, P1]

def conditional_entropy_after_measurement(rho_AB, theta, phi):
    """S(A|B) after projective measurement on B parameterized by (theta, phi).
    rho_AB is a 4x4 density matrix with A=first qubit, B=second qubit."""
    projectors = measurement_projectors(theta, phi)
    S_cond = 0.0
    for P in projectors:
        # (I_A tensor P_B) rho (I_A tensor P_B)
        I_A = np.eye(2)
        M = np.kron(I_A, P)
        rho_post = M @ rho_AB @ M
        p_k = np.real(np.trace(rho_post))
        if p_k < 1e-15:
            continue
        rho_post /= p_k
        # Trace out B to get rho_A|k
        rho_A_given_k = rho_post[:2, :2] + rho_post[2:, 2:]  # partial trace over B for 2-qubit
        # Actually need proper partial trace for 4x4
        # rho_A|k = Tr_B(rho_post) for a 2x2 system
        # For 4x4 matrix with basis |00>,|01>,|10>,|11>:
        # rho_A = [[rho[0,0]+rho[1,1], rho[0,2]+rho[1,3]], [rho[2,0]+rho[3,1], rho[2,2]+rho[3,3]]]
        rho_A_k = np.array([
            [rho_post[0, 0] + rho_post[1, 1], rho_post[0, 2] + rho_post[1, 3]],
            [rho_post[2, 0] + rho_post[3, 1], rho_post[2, 2] + rho_post[3, 3]]
        ])
        S_cond += p_k * von_neumann_entropy(np.real(rho_A_k))
    return S_cond

def quantum_discord(rho_AB):
    """Compute quantum discord D(A|B) = S(A) - S(AB) + min_{meas B} S(A|meas B).
    rho_AB is 4x4 with qubit ordering: A is first (rows 0,1 of tensor), B is second."""
    # S(A) — trace out B
    rho_A = np.array([
        [rho_AB[0, 0] + rho_AB[1, 1], rho_AB[0, 2] + rho_AB[1, 3]],
        [rho_AB[2, 0] + rho_AB[3, 1], rho_AB[2, 2] + rho_AB[3, 3]]
    ])
    S_A = von_neumann_entropy(np.real(rho_A))
    S_AB = von_neumann_entropy(np.real(rho_AB))

    # Minimize conditional entropy over measurement angles
    def objective(params):
        theta, phi = params
        return conditional_entropy_after_measurement(rho_AB, theta, phi)

    # Try multiple starting points
    best_val = float('inf')
    for t0_init in [0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi]:
        for p0_init in [0, np.pi/2, np.pi, 3*np.pi/2]:
            res = minimize(objective, [t0_init, p0_init],
                          method='Nelder-Mead',
                          options={'xatol': 1e-8, 'fatol': 1e-10, 'maxiter': 500})
            if res.fun < best_val:
                best_val = res.fun

    discord = S_A - S_AB + best_val
    classical_corr = S_A - best_val  # J(A|B) = S(A) - min S(A|meas B)

    # MI for comparison
    rho_B = np.array([
        [rho_AB[0, 0] + rho_AB[2, 2], rho_AB[0, 1] + rho_AB[2, 3]],
        [rho_AB[1, 0] + rho_AB[3, 2], rho_AB[1, 1] + rho_AB[3, 3]]
    ])
    S_B = von_neumann_entropy(np.real(rho_B))
    MI = S_A + S_B - S_AB

    return discord, classical_corr, MI

# Build states
states = {
    "GHZ": make_ghz(N),
    "W": make_w(N),
    "Cluster": make_cluster(N),
}

results = {}

for state_name, sv in states.items():
    print(f"\n=== {state_name} (n={N}) — Quantum Discord ===")

    discord_matrix = np.zeros((N, N))
    classical_matrix = np.zeros((N, N))
    mi_matrix = np.zeros((N, N))
    quantum_fraction = np.zeros((N, N))

    for i in range(N):
        for j in range(i + 1, N):
            rho_ij = get_2qubit_rdm(sv, i, j, N)
            d, cc, mi = quantum_discord(rho_ij)
            discord_matrix[i, j] = discord_matrix[j, i] = d
            classical_matrix[i, j] = classical_matrix[j, i] = cc
            mi_matrix[i, j] = mi_matrix[j, i] = mi
            if mi > 1e-10:
                quantum_fraction[i, j] = quantum_fraction[j, i] = d / mi
            print(f"  ({i},{j}): MI={mi:.4f}  Discord={d:.4f}  Classical={cc:.4f}  Q_frac={d/mi:.3f}" if mi > 1e-10
                  else f"  ({i},{j}): MI={mi:.4f}  Discord={d:.4f}  (no correlation)")

    # Summary stats
    offdiag = discord_matrix[np.triu_indices(N, k=1)]
    mi_offdiag = mi_matrix[np.triu_indices(N, k=1)]
    qf_offdiag = quantum_fraction[np.triu_indices(N, k=1)]
    correlated = mi_offdiag > 1e-10

    print(f"\n  Mean discord: {np.mean(offdiag):.4f}")
    print(f"  Mean MI: {np.mean(mi_offdiag):.4f}")
    if np.any(correlated):
        print(f"  Mean quantum fraction: {np.mean(qf_offdiag[correlated]):.4f}")

    results[state_name] = {
        "discord_matrix": discord_matrix.tolist(),
        "classical_matrix": classical_matrix.tolist(),
        "mi_matrix": mi_matrix.tolist(),
        "quantum_fraction_matrix": quantum_fraction.tolist(),
        "mean_discord": float(np.mean(offdiag)),
        "mean_mi": float(np.mean(mi_offdiag)),
        "mean_quantum_fraction": float(np.mean(qf_offdiag[correlated])) if np.any(correlated) else 0.0,
    }

dt = time.time() - t0

print(f"\n=== Summary ===")
for name in ["GHZ", "W", "Cluster"]:
    r = results[name]
    print(f"{name:8s}: mean_discord={r['mean_discord']:.4f}, mean_MI={r['mean_mi']:.4f}, "
          f"quantum_fraction={r['mean_quantum_fraction']:.3f}")

results["n_qubits"] = N
results["runtime_seconds"] = dt

with open(os.path.join(RESULTS_DIR, "sprint_005a_discord.json"), "w") as f:
    json.dump(results, f, indent=2)

print(f"\nRuntime: {dt:.2f}s")
print(f"Saved to results/sprint_005a_discord.json")
