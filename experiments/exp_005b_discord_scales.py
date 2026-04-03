"""
Experiment 5b: Discord at Different Scales

If cluster state quantum correlations are invisible pairwise but appear in I3,
do they show up as discord at the 1-vs-2 qubit level?

Compute D(A|BC) — discord between one qubit and a pair — for cluster, GHZ, W.
The measurement optimization is now over 2-qubit projective measurements on BC,
which is parameterized by a unitary on C^4 (much harder optimization).

Simpler approach: use the "geometric discord" approximation, or compute
D(A|BC) where we measure only B and C independently (local measurements).

We'll use local measurements: measure B with (theta_B, phi_B), measure C with
(theta_C, phi_C). This gives a lower bound on classical correlations and thus
an upper bound on discord. But for our comparison purposes, it's sufficient.
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
    eigenvalues = np.linalg.eigvalsh(rho_matrix)
    eigenvalues = eigenvalues[eigenvalues > 1e-15]
    return -np.sum(eigenvalues * np.log2(eigenvalues))

def get_3qubit_rdm(sv, i, j, k, n_total):
    """Get 3-qubit reduced density matrix for qubits i,j,k."""
    trace_out = [q for q in range(n_total) if q not in (i, j, k)]
    if not trace_out:
        rho = sv.to_operator()
        return np.array(rho)
    rho = partial_trace(sv, trace_out)
    return np.array(rho)

def measurement_projectors(theta, phi):
    cos_t2 = np.cos(theta / 2)
    sin_t2 = np.sin(theta / 2)
    n_ket = np.array([cos_t2, np.exp(1j * phi) * sin_t2])
    n_perp = np.array([sin_t2, -np.exp(1j * phi) * cos_t2])
    P0 = np.outer(n_ket, n_ket.conj())
    P1 = np.outer(n_perp, n_perp.conj())
    return [P0, P1]

def partial_trace_single(rho_8x8, trace_qubit):
    """Partial trace of a 3-qubit state (8x8) over one qubit.
    trace_qubit: 0, 1, or 2 (in the 3-qubit system A, B, C)."""
    result = np.zeros((4, 4), dtype=complex)
    if trace_qubit == 0:  # Trace out A, keep BC
        for a in range(2):
            for b1 in range(2):
                for c1 in range(2):
                    for b2 in range(2):
                        for c2 in range(2):
                            i = a * 4 + b1 * 2 + c1
                            j = a * 4 + b2 * 2 + c2
                            result[b1*2+c1, b2*2+c2] += rho_8x8[i, j]
    elif trace_qubit == 1:  # Trace out B, keep AC
        for b in range(2):
            for a1 in range(2):
                for c1 in range(2):
                    for a2 in range(2):
                        for c2 in range(2):
                            i = a1 * 4 + b * 2 + c1
                            j = a2 * 4 + b * 2 + c2
                            result[a1*2+c1, a2*2+c2] += rho_8x8[i, j]
    elif trace_qubit == 2:  # Trace out C, keep AB
        for c in range(2):
            for a1 in range(2):
                for b1 in range(2):
                    for a2 in range(2):
                        for b2 in range(2):
                            i = a1 * 4 + b1 * 2 + c
                            j = a2 * 4 + b2 * 2 + c
                            result[a1*2+b1, a2*2+b2] += rho_8x8[i, j]
    return result

def partial_trace_two(rho_8x8, keep_qubit):
    """Trace out 2 qubits, keep 1. Returns 2x2."""
    result = np.zeros((2, 2), dtype=complex)
    if keep_qubit == 0:  # Keep A
        for b in range(2):
            for c in range(2):
                for a1 in range(2):
                    for a2 in range(2):
                        result[a1, a2] += rho_8x8[a1*4+b*2+c, a2*4+b*2+c]
    elif keep_qubit == 1:  # Keep B
        for a in range(2):
            for c in range(2):
                for b1 in range(2):
                    for b2 in range(2):
                        result[b1, b2] += rho_8x8[a*4+b1*2+c, a*4+b2*2+c]
    elif keep_qubit == 2:  # Keep C
        for a in range(2):
            for b in range(2):
                for c1 in range(2):
                    for c2 in range(2):
                        result[c1, c2] += rho_8x8[a*4+b*2+c1, a*4+b*2+c2]
    return result

def discord_1vs2(rho_ABC, measure_qubits=(1, 2)):
    """Discord D(A|BC) with local measurements on B and C.
    rho_ABC is 8x8 with qubit ordering A=0, B=1, C=2.
    We measure B and C independently."""

    # S(A)
    rho_A = partial_trace_two(rho_ABC, keep_qubit=0)
    S_A = von_neumann_entropy(np.real(rho_A))

    # S(ABC)
    S_ABC = von_neumann_entropy(np.real(rho_ABC))

    def objective(params):
        theta_B, phi_B, theta_C, phi_C = params
        projs_B = measurement_projectors(theta_B, phi_B)
        projs_C = measurement_projectors(theta_C, phi_C)

        S_cond = 0.0
        for PB in projs_B:
            for PC in projs_C:
                # M = I_A tensor PB tensor PC
                M = np.kron(np.kron(np.eye(2), PB), PC)
                rho_post = M @ rho_ABC @ M
                p = np.real(np.trace(rho_post))
                if p < 1e-15:
                    continue
                rho_post /= p
                # Trace out B and C to get rho_A|outcome
                rho_A_cond = partial_trace_two(rho_post, keep_qubit=0)
                S_cond += p * von_neumann_entropy(np.real(rho_A_cond))
        return S_cond

    best_val = float('inf')
    for tB in [0, np.pi/2, np.pi]:
        for pB in [0, np.pi/2]:
            for tC in [0, np.pi/2, np.pi]:
                for pC in [0, np.pi/2]:
                    res = minimize(objective, [tB, pB, tC, pC],
                                  method='Nelder-Mead',
                                  options={'xatol': 1e-6, 'fatol': 1e-8, 'maxiter': 300})
                    if res.fun < best_val:
                        best_val = res.fun

    # Correct discord formula: D(A|BC) = S(BC) - S(ABC) + min S(A|meas on BC)
    rho_BC = partial_trace_single(rho_ABC, trace_qubit=0)
    S_BC = von_neumann_entropy(np.real(rho_BC))

    discord = S_BC - S_ABC + best_val

    # MI(A:BC) for comparison
    MI_A_BC = S_A + S_BC - S_ABC

    # Classical correlations: J(A:BC) = S(A) - min S(A|meas on BC)
    classical_corr = S_A - best_val

    return discord, MI_A_BC, classical_corr

states = {
    "GHZ": make_ghz(N),
    "W": make_w(N),
    "Cluster": make_cluster(N),
}

results = {}

# Test specific triples that are interesting:
# For cluster: consecutive triples (where I3 was negative)
# For comparison: same triples in GHZ and W
triples = [(0,1,2), (1,2,3), (2,3,4)]

for state_name, sv in states.items():
    print(f"\n=== {state_name} (n={N}) — Discord D(A|BC) ===")
    state_results = {}

    for (a, b, c) in triples:
        rho = get_3qubit_rdm(sv, a, b, c, N)
        d, mi, cc = discord_1vs2(rho)
        qfrac = d / mi if mi > 1e-10 else 0.0
        print(f"  D(q{a}|q{b}q{c}): discord={d:.4f}, MI={mi:.4f}, classical={cc:.4f}, Q_frac={qfrac:.3f}")
        state_results[f"({a},{b},{c})"] = {
            "discord": float(d),
            "MI_A_BC": float(mi),
            "classical_corr": float(cc),
            "quantum_fraction": float(qfrac),
        }

    results[state_name] = state_results

dt = time.time() - t0

print(f"\n=== Scale Comparison ===")
print("                 Pairwise discord    1-vs-2 discord")
print(f"GHZ:             0.000               {results['GHZ']['(1,2,3)']['discord']:.4f}")
print(f"W:               0.282               {results['W']['(1,2,3)']['discord']:.4f}")
print(f"Cluster:         0.000               {results['Cluster']['(1,2,3)']['discord']:.4f}")

results["triples_tested"] = [list(t) for t in triples]
results["n_qubits"] = N
results["runtime_seconds"] = dt

with open(os.path.join(RESULTS_DIR, "sprint_005b_discord_scales.json"), "w") as f:
    json.dump(results, f, indent=2)

print(f"\nRuntime: {dt:.2f}s")
print(f"Saved to results/sprint_005b_discord_scales.json")
