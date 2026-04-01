"""
Sprint 017c: Information structure at the error correction threshold.
The threshold theorem predicts a phase transition at p_th. Does the
information structure of encoded+noisy states show this transition?

Track: half-cut entropy, pairwise MI, I3, and logical MI (how much
information about the logical state survives in the physical qubits)
for the 3-qubit repetition code under bit-flip noise.

Then compare with the 9-qubit (2-level concatenated) code to see if
concatenation sharpens the transition.
"""
import numpy as np
import json
from qiskit.circuit import QuantumCircuit
from qiskit_aer import AerSimulator
from qiskit.quantum_info import DensityMatrix, partial_trace, entropy

def von_neumann_entropy(rho):
    """Von Neumann entropy in bits."""
    return float(entropy(rho, base=2))

def mutual_information(rho_full, qubits_a, qubits_b, n_total):
    """MI between two subsystems."""
    all_other_a = [i for i in range(n_total) if i not in qubits_a]
    all_other_b = [i for i in range(n_total) if i not in qubits_b]
    all_other_ab = [i for i in range(n_total) if i not in qubits_a and i not in qubits_b]

    rho_a = partial_trace(rho_full, all_other_a)
    rho_b = partial_trace(rho_full, all_other_b)
    rho_ab = partial_trace(rho_full, all_other_ab)

    return von_neumann_entropy(rho_a) + von_neumann_entropy(rho_b) - von_neumann_entropy(rho_ab)

def tripartite_info(rho_full, qa, qb, qc, n_total):
    """I3(A:B:C) = I(A:B) + I(A:C) - I(A:BC)."""
    iab = mutual_information(rho_full, qa, qb, n_total)
    iac = mutual_information(rho_full, qa, qc, n_total)
    iabc = mutual_information(rho_full, qa, qb + qc, n_total)
    return iab + iac - iabc

def get_noisy_density_matrix(n_qubits, init_state, p_bitflip):
    """Get density matrix of encoded state after bit-flip noise.
    Uses density_matrix simulation method for exact results."""
    sim = AerSimulator(method='density_matrix')
    qc = QuantumCircuit(n_qubits)

    if init_state == 1:
        qc.x(0)

    # Encode: repetition
    for i in range(1, n_qubits):
        qc.cx(0, i)

    # Apply bit-flip noise via Kraus operators:
    # We'll use the save_density_matrix approach
    # But noise model needs measurement... Let's compute analytically instead.

    # For repetition code under bit-flip noise:
    # |0_L> = |000...0>, |1_L> = |111...1>
    # After independent bit-flip with probability p on each qubit:
    # rho = sum over all bit-flip patterns of their probabilities * resulting state

    # For |0_L> = |00...0>:
    # Each qubit flips independently with prob p
    # rho = sum_{x in {0,1}^n} p^|x| * (1-p)^(n-|x|) * |x><x|
    # This is a diagonal density matrix!

    # For |1_L> = |11...1>:
    # rho = sum_{x} p^|x_bar| * (1-p)^(n-|x_bar|) * |x><x| where x_bar = bitwise complement

    dim = 2**n_qubits
    rho = np.zeros((dim, dim), dtype=complex)

    for x in range(dim):
        bits = [(x >> i) & 1 for i in range(n_qubits)]
        if init_state == 0:
            n_flips = sum(bits)  # number of 1s = number of flips from |000>
        else:
            n_flips = n_qubits - sum(bits)  # number of 0s = flips from |111>

        prob = (p_bitflip ** n_flips) * ((1 - p_bitflip) ** (n_qubits - n_flips))
        rho[x, x] = prob

    return DensityMatrix(rho)

def logical_fidelity_from_rho(rho, init_state, n_qubits):
    """Compute logical fidelity using majority vote on the density matrix."""
    dim = 2**n_qubits
    fidelity = 0.0
    for x in range(dim):
        bits = [(x >> i) & 1 for i in range(n_qubits)]
        majority = 1 if sum(bits) > n_qubits / 2 else 0
        if majority == init_state:
            fidelity += rho.data[x, x].real
    return fidelity

def logical_mutual_information(p, n_qubits):
    """MI between logical qubit L and physical block P.

    Prepare L uniformly: P(L=0) = P(L=1) = 1/2.
    rho_LP = 1/2 |0><0|_L x rho_P(0) + 1/2 |1><1|_L x rho_P(1)

    MI(L:P) = H(L) + H(P) - H(LP)
    H(L) = 1 bit always.
    H(P) = entropy of (rho_P(0) + rho_P(1))/2
    H(LP) = 1 + (H(rho_P(0)) + H(rho_P(1)))/2  ... no, need to be more careful.

    Actually: rho_LP is block diagonal in L basis:
    rho_LP = 1/2 * diag(rho_P(0), rho_P(1))  in the {|0_L>, |1_L>} x P basis
    H(LP) = 1 + (H(rho_P(0)) + H(rho_P(1)))/2  since the blocks are orthogonal
    Wait, H of a block diagonal matrix = sum of (p_i * H(block_i/p_i)) + H(p_i)
    = 1/2 * H(2*rho_P(0)) + 1/2 * H(2*rho_P(1)) + H(1/2,1/2)
    = 1/2 * (H(rho_P(0)) - log2(1/2)) + ... no.

    For block diagonal: rho = p0 * |0><0| x sigma0 + p1 * |1><1| x sigma1
    eigenvalues are {p0 * eig(sigma0)} union {p1 * eig(sigma1)}
    H(rho) = -sum p0*lambda_i*log(p0*lambda_i) - sum p1*mu_j*log(p1*mu_j)
    = p0 * H(sigma0) - p0*log(p0) + p1 * H(sigma1) - p1*log(p1)
    = p0*H(sigma0) + p1*H(sigma1) + H(p0,p1)

    So H(LP) = 1/2 * H(rho_P(0)) + 1/2 * H(rho_P(1)) + 1  (since H(1/2,1/2) = 1)

    MI(L:P) = H(L) + H(P) - H(LP)
    = 1 + H((rho_P(0)+rho_P(1))/2) - (1/2*H(rho_P(0)) + 1/2*H(rho_P(1)) + 1)
    = H((rho_P(0)+rho_P(1))/2) - 1/2*H(rho_P(0)) - 1/2*H(rho_P(1))
    = Holevo information!
    """
    rho0 = get_noisy_density_matrix(n_qubits, 0, p)
    rho1 = get_noisy_density_matrix(n_qubits, 1, p)

    # Average state
    rho_avg = DensityMatrix((rho0.data + rho1.data) / 2)

    h_avg = von_neumann_entropy(rho_avg)
    h0 = von_neumann_entropy(rho0)
    h1 = von_neumann_entropy(rho1)

    # Holevo bound = MI(L:P)
    return h_avg - 0.5 * h0 - 0.5 * h1

# ============================================================
# Main experiment: information measures vs noise for 3-qubit code
# ============================================================

p_values = [0.0, 0.02, 0.05, 0.08, 0.10, 0.15, 0.20, 0.25, 0.30,
            0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.70, 0.80, 0.90, 1.0]

print("=== 3-QUBIT REPETITION CODE: INFORMATION AT THRESHOLD ===\n")

results_3q = {
    "half_cut_entropy_0": [], "half_cut_entropy_1": [],
    "mi_01": [], "mi_02": [], "mi_12": [],
    "i3_012": [],
    "logical_mi": [],
    "logical_fidelity_0": [], "logical_fidelity_1": [],
    "total_entropy_0": [], "total_entropy_1": [],
}

for p in p_values:
    print(f"  p={p:.2f}...")

    rho0 = get_noisy_density_matrix(3, 0, p)
    rho1 = get_noisy_density_matrix(3, 1, p)

    # Total entropy
    results_3q["total_entropy_0"].append(round(von_neumann_entropy(rho0), 6))
    results_3q["total_entropy_1"].append(round(von_neumann_entropy(rho1), 6))

    # Half-cut entropy (trace out qubit 2, entropy of [0,1])
    rho0_01 = partial_trace(rho0, [2])
    rho1_01 = partial_trace(rho1, [2])
    results_3q["half_cut_entropy_0"].append(round(von_neumann_entropy(rho0_01), 6))
    results_3q["half_cut_entropy_1"].append(round(von_neumann_entropy(rho1_01), 6))

    # Pairwise MI (use |0> state as representative)
    mi01 = mutual_information(rho0, [0], [1], 3)
    mi02 = mutual_information(rho0, [0], [2], 3)
    mi12 = mutual_information(rho0, [1], [2], 3)
    results_3q["mi_01"].append(round(mi01, 6))
    results_3q["mi_02"].append(round(mi02, 6))
    results_3q["mi_12"].append(round(mi12, 6))

    # I3
    i3 = tripartite_info(rho0, [0], [1], [2], 3)
    results_3q["i3_012"].append(round(i3, 6))

    # Logical MI (Holevo information)
    lmi = logical_mutual_information(p, 3)
    results_3q["logical_mi"].append(round(lmi, 6))

    # Logical fidelity
    fid0 = logical_fidelity_from_rho(rho0, 0, 3)
    fid1 = logical_fidelity_from_rho(rho1, 1, 3)
    results_3q["logical_fidelity_0"].append(round(fid0, 6))
    results_3q["logical_fidelity_1"].append(round(fid1, 6))

# Print 3-qubit results
print(f"\n{'p':>5} | {'S_tot':>7} | {'S_half':>7} | {'MI(0,1)':>8} | {'I3':>7} | {'L_MI':>7} | {'L_Fid':>7}")
print("-" * 65)
for i, p in enumerate(p_values):
    print(f"{p:5.2f} | {results_3q['total_entropy_0'][i]:7.4f} | {results_3q['half_cut_entropy_0'][i]:7.4f} | "
          f"{results_3q['mi_01'][i]:8.4f} | {results_3q['i3_012'][i]:7.4f} | "
          f"{results_3q['logical_mi'][i]:7.4f} | {results_3q['logical_fidelity_0'][i]:7.4f}")

# ============================================================
# 9-qubit (2-level concatenated) code
# ============================================================
print("\n=== 9-QUBIT CONCATENATED CODE: INFORMATION AT THRESHOLD ===\n")

# Coarser sweep for 9 qubits (512-dim density matrix)
p_values_9q = [0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 1.0]

results_9q = {
    "logical_mi": [],
    "logical_fidelity_0": [], "logical_fidelity_1": [],
    "total_entropy_0": [],
    "half_cut_entropy_0": [],
}

def logical_fidelity_9q(rho, init_state):
    """2-level concatenated majority vote."""
    dim = 2**9
    fidelity = 0.0
    for x in range(dim):
        bits = [(x >> i) & 1 for i in range(9)]
        # Inner decode: blocks of 3
        b0 = 1 if sum(bits[0:3]) > 1 else 0
        b1 = 1 if sum(bits[3:6]) > 1 else 0
        b2 = 1 if sum(bits[6:9]) > 1 else 0
        # Outer decode
        decoded = 1 if (b0 + b1 + b2) > 1 else 0
        if decoded == init_state:
            fidelity += rho.data[x, x].real
    return fidelity

for p in p_values_9q:
    print(f"  p={p:.2f}...")

    rho0 = get_noisy_density_matrix(9, 0, p)
    rho1 = get_noisy_density_matrix(9, 1, p)

    results_9q["total_entropy_0"].append(round(von_neumann_entropy(rho0), 6))

    # Half-cut: trace out qubits 5-8, keep 0-4
    rho0_half = partial_trace(rho0, [5, 6, 7, 8])
    results_9q["half_cut_entropy_0"].append(round(von_neumann_entropy(rho0_half), 6))

    # Logical MI
    rho_avg = DensityMatrix((rho0.data + rho1.data) / 2)
    lmi = von_neumann_entropy(rho_avg) - 0.5 * von_neumann_entropy(rho0) - 0.5 * von_neumann_entropy(rho1)
    results_9q["logical_mi"].append(round(lmi, 6))

    # Logical fidelity
    fid0 = logical_fidelity_9q(rho0, 0)
    fid1 = logical_fidelity_9q(rho1, 1)
    results_9q["logical_fidelity_0"].append(round(fid0, 6))
    results_9q["logical_fidelity_1"].append(round(fid1, 6))

print(f"\n{'p':>5} | {'S_tot':>7} | {'S_half':>7} | {'L_MI':>7} | {'L_Fid':>7}")
print("-" * 45)
for i, p in enumerate(p_values_9q):
    print(f"{p:5.2f} | {results_9q['total_entropy_0'][i]:7.4f} | {results_9q['half_cut_entropy_0'][i]:7.4f} | "
          f"{results_9q['logical_mi'][i]:7.4f} | {results_9q['logical_fidelity_0'][i]:7.4f}")

# ============================================================
# Compare: uncoded single qubit logical MI
# ============================================================
print("\n=== UNCODED LOGICAL MI ===")
uncoded_lmi = []
for p in p_values:
    rho0 = get_noisy_density_matrix(1, 0, p)
    rho1 = get_noisy_density_matrix(1, 1, p)
    rho_avg = DensityMatrix((rho0.data + rho1.data) / 2)
    lmi = von_neumann_entropy(rho_avg) - 0.5 * von_neumann_entropy(rho0) - 0.5 * von_neumann_entropy(rho1)
    uncoded_lmi.append(round(lmi, 6))
    print(f"  p={p:.2f}: LMI={lmi:.6f}")

# Save
all_results = {
    "p_values_3q": p_values,
    "p_values_9q": p_values_9q,
    "three_qubit": results_3q,
    "nine_qubit": results_9q,
    "uncoded_logical_mi": uncoded_lmi,
    "threshold": 0.5,
    "notes": "Bit-flip noise on repetition code. Logical MI = Holevo info = distinguishability of encoded |0> vs |1>."
}

with open("results/sprint_017c_threshold_information.json", "w") as f:
    json.dump(all_results, f, indent=2)

print("\nResults saved to results/sprint_017c_threshold_information.json")
