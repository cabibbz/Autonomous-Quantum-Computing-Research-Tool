"""Sprint 015b: Error correction performance comparison.

Compares [[5,1,3]], Steane [[7,1,3]], and Shor [[9,1,3]] under depolarizing noise.
Measures logical fidelity vs physical error rate. Finds break-even thresholds.
"""

import numpy as np
import json
import time

# Pauli matrices
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
pauli_map = {'I': I2, 'X': X, 'Y': Y, 'Z': Z}

def kron_list(mats):
    result = mats[0]
    for m in mats[1:]:
        result = np.kron(result, m)
    return result

def pauli_on_qubit(pauli, qubit, n):
    ops = [I2] * n
    ops[qubit] = pauli
    return kron_list(ops)

def pauli_string(s):
    return kron_list([pauli_map[c] for c in s])

def depolarizing_channel(rho, p, n):
    for q in range(n):
        XI = pauli_on_qubit(X, q, n)
        YI = pauli_on_qubit(Y, q, n)
        ZI = pauli_on_qubit(Z, q, n)
        rho = (1 - p) * rho + (p / 3) * (XI @ rho @ XI + YI @ rho @ YI + ZI @ rho @ ZI)
    return rho

def fidelity(rho, target_state):
    return float(np.real(target_state.conj() @ rho @ target_state))


# ============================================================
# Build encoding isometry from stabilizers
# ============================================================
def build_code(stabilizers, n, logical_x_str):
    dim = 2**n
    state_0 = np.zeros(dim, dtype=complex)
    state_0[0] = 1.0
    for stab in stabilizers:
        S = pauli_string(stab)
        proj = (np.eye(dim) + S) / 2.0
        state_0 = proj @ state_0
        norm = np.linalg.norm(state_0)
        if norm < 1e-12:
            raise ValueError(f"Stabilizer {stab} annihilated state!")
        state_0 = state_0 / norm

    X_L = pauli_string(logical_x_str)
    state_1 = X_L @ state_0
    state_1 = state_1 / np.linalg.norm(state_1)

    # Verify orthogonality
    overlap = abs(np.dot(state_0.conj(), state_1))
    if overlap > 1e-10:
        raise ValueError(f"Logical states not orthogonal! overlap={overlap:.6f}")

    V = np.zeros((dim, 2), dtype=complex)
    V[:, 0] = state_0
    V[:, 1] = state_1
    return V


# ============================================================
# Syndrome correction (general, on-the-fly projectors)
# ============================================================
def build_error_table(stabilizers, n):
    """Map syndromes to correction operators."""
    dim = 2**n
    stab_ops = [pauli_string(s) for s in stabilizers]
    num_stabs = len(stab_ops)
    table = {0: np.eye(dim, dtype=complex)}

    for q in range(n):
        for P in [X, Y, Z]:
            error = pauli_on_qubit(P, q, n)
            s = 0
            for k, S in enumerate(stab_ops):
                comm_test = S @ error - error @ S
                if np.linalg.norm(comm_test) > 1e-10:
                    s |= (1 << k)
            if s not in table:
                table[s] = error
    return table, stab_ops

def syndrome_correct(rho, stab_ops, error_table, n):
    """Apply syndrome correction, computing projectors on the fly."""
    dim = 2**n
    num_stabs = len(stab_ops)
    rho_corrected = np.zeros((dim, dim), dtype=complex)

    for syndrome in range(2**num_stabs):
        # Build projector
        proj = np.eye(dim, dtype=complex)
        for k, S in enumerate(stab_ops):
            sign = 1 if (syndrome >> k) & 1 == 0 else -1
            proj = (np.eye(dim) + sign * S) / 2.0 @ proj

        correction = error_table.get(syndrome, np.eye(dim, dtype=complex))
        projected = proj @ rho @ proj.conj().T
        rho_corrected += correction @ projected @ correction.conj().T

    return rho_corrected


# ============================================================
# Define codes
# ============================================================
print("=== Building codes ===")

# [[5,1,3]]
t0 = time.time()
V5 = build_code(['XZZXI', 'IXZZX', 'XIXZZ', 'ZXIXZ'], 5, 'XXXXX')
et5, so5 = build_error_table(['XZZXI', 'IXZZX', 'XIXZZ', 'ZXIXZ'], 5)
print(f"[[5,1,3]]: built in {time.time()-t0:.2f}s, {len(et5)} syndromes mapped")

# Steane [[7,1,3]]
t0 = time.time()
stab_7 = ['IIIXXXX', 'IXXIIXX', 'XIXIXIX', 'IIIZZZZ', 'IZZIIZZ', 'ZIZIZIZ']
V7 = build_code(stab_7, 7, 'XXXXXXX')
et7, so7 = build_error_table(stab_7, 7)
print(f"Steane [[7,1,3]]: built in {time.time()-t0:.2f}s, {len(et7)} syndromes mapped")

# Shor [[9,1,3]] — build directly (stabilizer projection gives |+>_L not |0>_L)
t0 = time.time()
stab_9 = [
    'ZZIIIIIII', 'IZZIIIIII',
    'IIIZZIIII', 'IIIIZZIII',
    'IIIIIIZZI', 'IIIIIIIZZ',
    'XXXXXXIII', 'IIIXXXXXX',
]
# Build Shor code states directly: |0>_L = |+>|+>|+>, |1>_L = |->|->|-> (block GHZ)
plus_block = np.zeros(8, dtype=complex)
plus_block[0b000] = 1/np.sqrt(2); plus_block[0b111] = 1/np.sqrt(2)
minus_block = np.zeros(8, dtype=complex)
minus_block[0b000] = 1/np.sqrt(2); minus_block[0b111] = -1/np.sqrt(2)
state_0_shor = np.kron(np.kron(plus_block, plus_block), plus_block)
state_1_shor = np.kron(np.kron(minus_block, minus_block), minus_block)
V9 = np.zeros((512, 2), dtype=complex)
V9[:, 0] = state_0_shor; V9[:, 1] = state_1_shor
# Verify orthogonality
assert abs(np.dot(state_0_shor.conj(), state_1_shor)) < 1e-10, "Shor code states not orthogonal!"
et9, so9 = build_error_table(stab_9, 9)
print(f"Shor [[9,1,3]]: built in {time.time()-t0:.2f}s, {len(et9)} syndromes mapped")


# ============================================================
# Timing test for Shor correction
# ============================================================
print("\n=== Timing test (Shor syndrome correction) ===")
logical_sv = np.array([1, 0], dtype=complex)
logical_rho = np.outer(logical_sv, logical_sv.conj())
rho_coded = V9 @ logical_rho @ V9.conj().T
rho_noisy = depolarizing_channel(rho_coded, 0.1, 9)
t0 = time.time()
rho_corrected = syndrome_correct(rho_noisy, so9, et9, 9)
t1 = time.time()
print(f"Shor correction: {t1-t0:.2f}s")


# ============================================================
# Run performance comparison
# ============================================================
p_values = np.linspace(0, 0.3, 16)  # Focus on interesting range
logical_states = {
    'zero': np.array([1, 0], dtype=complex),
    'plus': np.array([1, 1], dtype=complex) / np.sqrt(2),
}

results = {}

for code_name, V, et, so, n in [
    ('five_qubit', V5, et5, so5, 5),
    ('steane', V7, et7, so7, 7),
    ('shor', V9, et9, so9, 9),
]:
    print(f"\n=== {code_name} (n={n}) under depolarizing noise ===")
    for lname, logical_sv in logical_states.items():
        logical_rho = np.outer(logical_sv, logical_sv.conj())

        fid_uncoded = []
        fid_coded = []

        t0 = time.time()
        for p in p_values:
            # Uncoded
            rho_phys = np.array(logical_rho, dtype=complex)
            rho_phys_noisy = (1 - p) * rho_phys + (p/3) * (X @ rho_phys @ X + Y @ rho_phys @ Y + Z @ rho_phys @ Z)
            f_unc = fidelity(rho_phys_noisy, logical_sv)

            # Coded + corrected
            rho_coded = V @ logical_rho @ V.conj().T
            rho_noisy = depolarizing_channel(rho_coded, p, n)
            rho_corrected = syndrome_correct(rho_noisy, so, et, n)
            rho_decoded = V.conj().T @ rho_corrected @ V
            tr = np.trace(rho_decoded).real
            if tr > 1e-12:
                rho_decoded = rho_decoded / tr
            f_cod = fidelity(rho_decoded, logical_sv)

            fid_uncoded.append(round(float(f_unc), 6))
            fid_coded.append(round(float(f_cod), 6))

        elapsed = time.time() - t0
        results[f'{code_name}_{lname}'] = {
            'n': n,
            'p_values': [round(float(p), 4) for p in p_values],
            'fid_uncoded': fid_uncoded,
            'fid_coded': fid_coded,
        }

        # Find break-even
        breakeven = None
        for i, p in enumerate(p_values):
            if p > 0 and fid_coded[i] < fid_uncoded[i]:
                breakeven = round(float(p), 4)
                break

        print(f"  |{lname}>: break-even ~{breakeven}, time={elapsed:.1f}s")
        if len(p_values) > 5:
            i5 = 5
            print(f"    p={p_values[i5]:.3f}: uncoded={fid_uncoded[i5]:.4f}, coded={fid_coded[i5]:.4f}")
        if len(p_values) > 10:
            i10 = 10
            print(f"    p={p_values[i10]:.3f}: uncoded={fid_uncoded[i10]:.4f}, coded={fid_coded[i10]:.4f}")

# ============================================================
# Summary comparison
# ============================================================
print("\n\n=== BREAK-EVEN COMPARISON (|zero>) ===")
print(f"{'p':<8}", end='')
for code_name in ['five_qubit', 'steane', 'shor']:
    print(f"{'coded_'+code_name:<18}", end='')
print(f"{'uncoded':<12}")
print("-" * 70)

for i, p in enumerate(p_values):
    print(f"{p:.3f}   ", end='')
    for code_name in ['five_qubit', 'steane', 'shor']:
        key = f'{code_name}_zero'
        f = results[key]['fid_coded'][i]
        print(f"{f:.4f}            ", end='')
    print(f"{results['five_qubit_zero']['fid_uncoded'][i]:.4f}")

# Save
with open('results/sprint_015b_code_families_performance.json', 'w') as f:
    json.dump(results, f, indent=2)
print("\nResults saved to results/sprint_015b_code_families_performance.json")
