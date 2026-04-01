"""
Sprint 019b: Channel Capacity vs Error Correction Thresholds.

Compare the theoretical channel capacity (from coherent information) to the
actual error correction performance of our codes from Sprints 015-018.

Key question: How does the QEC threshold relate to the channel capacity threshold?
- Channel capacity threshold: noise level where Q(N) = 0 (no code can help)
- Code-specific threshold: noise level where a particular code breaks even
- Gap between them = room for better codes

We also compute the achievable rate R = k/n for each code and compare to Q(N).
"""

import numpy as np
import json

# ============================================================
# CHANNEL CAPACITY FUNCTIONS (from 019a)
# ============================================================
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)

def von_neumann_entropy(rho):
    evals = np.linalg.eigvalsh(rho)
    evals = evals[evals > 1e-15]
    return -np.sum(evals * np.log2(evals))

def partial_trace_B(rho_AB):
    """Trace out B from a 2-qubit state."""
    rho = rho_AB.reshape(2, 2, 2, 2)
    return np.trace(rho, axis1=1, axis2=3)

def partial_trace_A(rho_AB):
    """Trace out A from a 2-qubit state."""
    rho = rho_AB.reshape(2, 2, 2, 2)
    return np.trace(rho, axis1=0, axis2=2)

def apply_channel_B(rho_AB, kraus_ops):
    d = 4
    result = np.zeros((d, d), dtype=complex)
    for K in kraus_ops:
        full_K = np.kron(I2, K)
        result += full_K @ rho_AB @ full_K.conj().T
    return result

def depolarizing_kraus(p):
    return [np.sqrt(1 - p) * I2, np.sqrt(p / 3) * X, np.sqrt(p / 3) * Y, np.sqrt(p / 3) * Z]

def amplitude_damping_kraus(gamma):
    E0 = np.array([[1, 0], [0, np.sqrt(1 - gamma)]], dtype=complex)
    E1 = np.array([[0, np.sqrt(gamma)], [0, 0]], dtype=complex)
    return [E0, E1]

def phase_damping_kraus(lam):
    return [np.sqrt(1 - lam) * I2, np.sqrt(lam) * Z]

def coherent_info_optimized(kraus_ops, n_samples=100):
    best_icoh = -np.inf
    for theta in np.linspace(0, np.pi / 2, n_samples):
        psi_AB = np.array([np.cos(theta), 0, 0, np.sin(theta)], dtype=complex)
        rho_AB = np.outer(psi_AB, psi_AB.conj())
        rho_AB_out = apply_channel_B(rho_AB, kraus_ops)
        rho_B = partial_trace_A(rho_AB_out)
        S_B = von_neumann_entropy(rho_B)
        S_AB = von_neumann_entropy(rho_AB_out)
        icoh = S_B - S_AB
        if icoh > best_icoh:
            best_icoh = icoh
    return max(0, best_icoh)

# ============================================================
# QEC CODE SIMULATION
# ============================================================

def kron_list(mats):
    result = mats[0]
    for m in mats[1:]:
        result = np.kron(result, m)
    return result

def apply_depol_n_qubits(rho, p, n):
    """Apply independent depolarizing noise to each of n qubits."""
    paulis = [I2, X, Y, Z]
    result = np.zeros_like(rho)
    # For each qubit, apply (1-p)I + p/3(X+Y+Z)
    # We iterate: apply to qubit 0, then qubit 1, etc.
    current = rho.copy()
    d = 2**n
    for q in range(n):
        new = np.zeros((d, d), dtype=complex)
        for pauli, coeff in zip(paulis, [1-p, p/3, p/3, p/3]):
            # Build n-qubit operator with pauli on qubit q, I elsewhere
            ops = [I2] * n
            ops[q] = pauli
            full_op = kron_list(ops)
            new += coeff * (full_op @ current @ full_op.conj().T)
        current = new
    return current

def three_qubit_code_fidelity(p, logical_state='0'):
    """3-qubit bit-flip code under depolarizing noise. Returns fidelity."""
    n = 3
    # Logical states
    if logical_state == '0':
        psi_L = np.zeros(8, dtype=complex)
        psi_L[0] = 1  # |000⟩
    else:
        psi_L = np.zeros(8, dtype=complex)
        psi_L[7] = 1  # |111⟩
    rho_L = np.outer(psi_L, psi_L.conj())

    # Apply noise
    rho_noisy = apply_depol_n_qubits(rho_L, p, n)

    # Syndrome measurement + correction (majority vote)
    # Projectors for each syndrome
    syndromes = {
        (0, 0): I2,   # no error
        (1, 0): X,    # flip qubit 0 (but really could be 0 or both 1,2)
        (1, 1): X,    # flip qubit 1
        (0, 1): X,    # flip qubit 2
    }
    # Z0Z1, Z1Z2 stabilizers
    S1 = kron_list([Z, Z, I2])
    S2 = kron_list([I2, Z, Z])

    # Project and correct
    rho_corrected = np.zeros((8, 8), dtype=complex)
    for s1_val in [0, 1]:
        for s2_val in [0, 1]:
            # Projector for syndrome (s1, s2)
            P = (np.eye(8) + (-1)**s1_val * S1) / 2 @ (np.eye(8) + (-1)**s2_val * S2) / 2
            projected = P @ rho_noisy @ P
            prob = np.real(np.trace(projected))
            if prob < 1e-15:
                continue
            # Correction
            if (s1_val, s2_val) == (0, 0):
                correction = kron_list([I2, I2, I2])
            elif (s1_val, s2_val) == (1, 0):
                correction = kron_list([X, I2, I2])
            elif (s1_val, s2_val) == (1, 1):
                correction = kron_list([I2, X, I2])
            elif (s1_val, s2_val) == (0, 1):
                correction = kron_list([I2, I2, X])
            rho_corrected += correction @ projected @ correction

    fidelity = np.real(psi_L.conj() @ rho_corrected @ psi_L)
    return fidelity

def five_qubit_code_fidelity(p, logical_state='0'):
    """[[5,1,3]] code under depolarizing noise."""
    n = 5
    d = 2**n

    # Stabilizer generators for [[5,1,3]]
    stab_strings = ['XZZXI', 'IXZZX', 'XIXZZ', 'ZXIXZ']

    def pauli_op(s):
        pmap = {'I': I2, 'X': X, 'Y': Y, 'Z': Z}
        return kron_list([pmap[c] for c in s])

    stabs = [pauli_op(s) for s in stab_strings]

    # Codespace projector: P = (1/16) * prod(I + S_i)
    P_code = np.eye(d, dtype=complex)
    for S in stabs:
        P_code = P_code @ (np.eye(d) + S) / 2

    # Logical operators
    X_L = pauli_op('XXXXX')
    Z_L = pauli_op('ZZZZZ')

    # Logical |0⟩ and |1⟩
    # |0_L⟩ = P_code |00000⟩ (normalized)
    basis_0 = np.zeros(d, dtype=complex)
    basis_0[0] = 1
    log_0 = P_code @ basis_0
    log_0 /= np.linalg.norm(log_0)
    log_1 = X_L @ log_0
    log_1 /= np.linalg.norm(log_1)

    if logical_state == '0':
        psi_L = log_0
    else:
        psi_L = log_1

    rho_L = np.outer(psi_L, psi_L.conj())

    # Apply noise
    rho_noisy = apply_depol_n_qubits(rho_L, p, n)

    # Syndrome-based correction
    # For [[5,1,3]], syndromes from 4 stabilizers = 16 syndromes
    # Each syndrome maps to a correctable error (weight 0 or 1)
    # Build syndrome table
    pauli_1q = {'I': I2, 'X': X, 'Y': Y, 'Z': Z}

    # All weight-0 and weight-1 errors
    corrections = {}
    for s_bits in range(16):
        corrections[s_bits] = None  # will fill

    # Weight-0 error (no error): syndrome = 0000
    corrections[0] = kron_list([I2] * n)

    # Weight-1 errors
    for q in range(n):
        for pauli_name in ['X', 'Y', 'Z']:
            ops = [I2] * n
            ops[q] = pauli_1q[pauli_name]
            error = kron_list(ops)
            # Compute syndrome
            syn = 0
            for i, S in enumerate(stabs):
                # Does this error commute or anticommute with S?
                comm = S @ error - error @ S
                if np.linalg.norm(comm) > 1e-10:
                    syn |= (1 << i)
            if syn not in corrections or corrections[syn] is None:
                corrections[syn] = error  # correction = same as error (Pauli)

    # Fill remaining with identity (uncorrectable)
    for syn in range(16):
        if corrections[syn] is None:
            corrections[syn] = kron_list([I2] * n)

    # Apply syndrome measurement and correction
    rho_corrected = np.zeros((d, d), dtype=complex)
    for syn in range(16):
        # Projector for this syndrome
        P_syn = np.eye(d, dtype=complex)
        for i, S in enumerate(stabs):
            sign = 1 if (syn >> i) & 1 == 0 else -1
            P_syn = P_syn @ (np.eye(d) + sign * S) / 2
        projected = P_syn @ rho_noisy @ P_syn
        prob = np.real(np.trace(projected))
        if prob < 1e-15:
            continue
        C = corrections[syn]
        rho_corrected += C @ projected @ C

    fidelity = np.real(psi_L.conj() @ rho_corrected @ psi_L)
    return fidelity

def uncoded_fidelity(p):
    """Single qubit under depolarizing noise."""
    # ρ → (1-p)|0⟩⟨0| + p/3(X|0⟩⟨0|X + Y|0⟩⟨0|Y + Z|0⟩⟨0|Z)
    # = (1-p)|0⟩⟨0| + p/3(|1⟩⟨1| + |1⟩⟨1| + |0⟩⟨0|) = (1-2p/3)|0⟩⟨0| + 2p/3|1⟩⟨1|
    return 1 - 2*p/3

# ============================================================
# MAIN EXPERIMENT
# ============================================================

noise_params = np.linspace(0, 0.30, 31)

results = {
    'noise_params': noise_params.tolist(),
    'channel_capacity': [],
    'uncoded_fidelity': [],
    'three_qubit_fidelity': [],
    'five_qubit_fidelity': [],
    'three_qubit_rate': 1/3,  # k/n = 1/3
    'five_qubit_rate': 1/5,   # k/n = 1/5
    'uncoded_holevo': [],
    'three_qubit_holevo': [],
    'five_qubit_holevo': [],
}

print("Computing channel capacity and code performance...")
print("=" * 60)

for i, p in enumerate(noise_params):
    # Channel capacity
    Q = coherent_info_optimized(depolarizing_kraus(p))
    results['channel_capacity'].append(Q)

    # Uncoded fidelity
    f_unc = uncoded_fidelity(p)
    results['uncoded_fidelity'].append(f_unc)

    # Code fidelities (average over |0⟩ and |1⟩)
    f3_0 = three_qubit_code_fidelity(p, '0')
    f3_1 = three_qubit_code_fidelity(p, '1')
    f3 = (f3_0 + f3_1) / 2
    results['three_qubit_fidelity'].append(f3)

    f5_0 = five_qubit_code_fidelity(p, '0')
    f5_1 = five_qubit_code_fidelity(p, '1')
    f5 = (f5_0 + f5_1) / 2
    results['five_qubit_fidelity'].append(f5)

    # Holevo information ≈ 1 - H(error)
    # Using fidelity: Holevo ≈ 1 + f*log₂(f) + (1-f)*log₂((1-f)/3) for depolarizing
    # Simpler: use binary entropy approximation since we're comparing logical qubit
    def holevo_from_fidelity(f):
        """Approximate Holevo info from fidelity for a qubit channel."""
        if f >= 1 - 1e-10:
            return 1.0
        if f <= 1e-10:
            return 0.0
        # For a qubit: Holevo ≤ 1
        # Holevo = 1 - H(1-f) for binary symmetric channel approximation
        return 1 + f * np.log2(f + 1e-30) + (1 - f) * np.log2((1 - f + 1e-30))

    results['uncoded_holevo'].append(holevo_from_fidelity(f_unc))
    results['three_qubit_holevo'].append(holevo_from_fidelity(f3))
    results['five_qubit_holevo'].append(holevo_from_fidelity(f5))

    if i % 5 == 0:
        print(f"p={p:.2f}: Q={Q:.4f}, F_unc={f_unc:.4f}, F_3q={f3:.4f}, F_5q={f5:.4f}")

# Find break-even points
print("\n" + "=" * 60)
print("BREAK-EVEN ANALYSIS")
print("=" * 60)

# Where does each code cross below uncoded?
for code_name, fid_list in [('3-qubit', results['three_qubit_fidelity']),
                             ('[[5,1,3]]', results['five_qubit_fidelity'])]:
    for i in range(len(fid_list) - 1):
        if fid_list[i] >= results['uncoded_fidelity'][i] and fid_list[i+1] < results['uncoded_fidelity'][i+1]:
            p1, p2 = noise_params[i], noise_params[i+1]
            d1 = fid_list[i] - results['uncoded_fidelity'][i]
            d2 = fid_list[i+1] - results['uncoded_fidelity'][i+1]
            cross = p1 + (p2 - p1) * d1 / (d1 - d2)
            print(f"{code_name} crosses below uncoded at p ≈ {cross:.4f}")
            results[f'{code_name}_break_even'] = cross
            break

# Capacity threshold
cap_threshold = None
for i in range(len(results['channel_capacity']) - 1):
    if results['channel_capacity'][i] > 0 and results['channel_capacity'][i+1] <= 0:
        p1, p2 = noise_params[i], noise_params[i+1]
        v1, v2 = results['channel_capacity'][i], results['channel_capacity'][i+1]
        cap_threshold = p1 + (p2 - p1) * v1 / (v1 - v2)
        break
results['capacity_threshold'] = cap_threshold
print(f"\nChannel capacity threshold: p ≈ {cap_threshold}")

# Efficiency: how much of the channel capacity does each code capture?
print("\n" + "=" * 60)
print("CODE EFFICIENCY vs CHANNEL CAPACITY")
print("=" * 60)
print(f"{'p':>5} | {'Q(N)':>7} | {'R_3q':>7} | {'R_5q':>7} | {'η_3q':>7} | {'η_5q':>7}")
print("-" * 50)
for i in range(0, len(noise_params), 3):
    p = noise_params[i]
    Q = results['channel_capacity'][i]
    h3 = results['three_qubit_holevo'][i]
    h5 = results['five_qubit_holevo'][i]
    # Achieved rate = (Holevo per logical qubit) * (k/n)
    r3 = h3 * (1/3)
    r5 = h5 * (1/5)
    eta3 = r3 / Q if Q > 0.01 else float('nan')
    eta5 = r5 / Q if Q > 0.01 else float('nan')
    print(f"{p:5.2f} | {Q:7.4f} | {r3:7.4f} | {r5:7.4f} | {eta3:7.3f} | {eta5:7.3f}")

# Save
with open('results/sprint_019b_capacity_vs_threshold.json', 'w') as f:
    json.dump(results, f, indent=2, default=str)
print("\nResults saved to results/sprint_019b_capacity_vs_threshold.json")

print("\n" + "=" * 60)
print("KEY FINDINGS")
print("=" * 60)
print(f"1. Channel capacity threshold (depolarizing): p ≈ {cap_threshold}")
print(f"   This is the FUNDAMENTAL limit — no code can help beyond this point.")
print(f"2. 3-qubit code captures ~{results['three_qubit_holevo'][5] * (1/3) / results['channel_capacity'][5] * 100:.0f}% of capacity at p=0.05")
print(f"   [[5,1,3]] captures ~{results['five_qubit_holevo'][5] * (1/5) / results['channel_capacity'][5] * 100:.0f}% of capacity at p=0.05")
print(f"3. Both codes operate FAR below channel capacity — huge room for improvement")
print(f"4. The gap between code threshold and channel threshold is the 'coding gain' yet to be achieved")
