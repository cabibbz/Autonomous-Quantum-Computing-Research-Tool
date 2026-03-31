"""Sprint 015c: Page curves for three code families.

Entangle a reference qubit with the logical qubit, encode into physical qubits,
then measure MI(reference : k physical qubits) for k=1..n.
Tests whether the Singleton bound (k >= n-d+1) predicts the recovery transition universally.

For each code [[n,1,d]], total system is n+1 qubits (1 reference + n physical).
We keep the reference (qubit 0) + k physical qubits and trace out the rest.
"""

import numpy as np
import json
import time
from itertools import combinations

# ============================================================
# Helpers
# ============================================================

def partial_trace(rho, keep, n):
    keep = sorted(keep)
    remove = sorted(set(range(n)) - set(keep))
    rho_t = rho.reshape([2]*n + [2]*n)
    for i, q in enumerate(remove):
        ax_row = q - i
        ax_col = ax_row + (n - i)
        rho_t = np.trace(rho_t, axis1=ax_row, axis2=ax_col)
    k = len(keep)
    return rho_t.reshape(2**k, 2**k)

def von_neumann_entropy(rho):
    evals = np.linalg.eigvalsh(rho)
    evals = evals[evals > 1e-12]
    return float(-np.sum(evals * np.log2(evals)))

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
pauli_map = {'I': I2, 'X': X, 'Y': Y, 'Z': Z}

def pauli_string(s):
    result = pauli_map[s[0]]
    for c in s[1:]:
        result = np.kron(result, pauli_map[c])
    return result

def kron_list(mats):
    result = mats[0]
    for m in mats[1:]:
        result = np.kron(result, m)
    return result

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

    overlap = abs(np.dot(state_0.conj(), state_1))
    if overlap > 1e-10:
        raise ValueError(f"Logical states not orthogonal! overlap={overlap:.6f}")

    V = np.zeros((dim, 2), dtype=complex)
    V[:, 0] = state_0
    V[:, 1] = state_1
    return V


# ============================================================
# Build codes
# ============================================================
print("=== Building codes ===")

# [[5,1,3]]
t0 = time.time()
V5 = build_code(['XZZXI', 'IXZZX', 'XIXZZ', 'ZXIXZ'], 5, 'XXXXX')
print(f"[[5,1,3]]: built in {time.time()-t0:.2f}s")

# Steane [[7,1,3]]
t0 = time.time()
stab_7 = ['IIIXXXX', 'IXXIIXX', 'XIXIXIX', 'IIIZZZZ', 'IZZIIZZ', 'ZIZIZIZ']
V7 = build_code(stab_7, 7, 'XXXXXXX')
print(f"Steane [[7,1,3]]: built in {time.time()-t0:.2f}s")

# Shor [[9,1,3]] — build directly
t0 = time.time()
plus_block = np.zeros(8, dtype=complex)
plus_block[0b000] = 1/np.sqrt(2); plus_block[0b111] = 1/np.sqrt(2)
minus_block = np.zeros(8, dtype=complex)
minus_block[0b000] = 1/np.sqrt(2); minus_block[0b111] = -1/np.sqrt(2)
state_0_shor = np.kron(np.kron(plus_block, plus_block), plus_block)
state_1_shor = np.kron(np.kron(minus_block, minus_block), minus_block)
V9 = np.zeros((512, 2), dtype=complex)
V9[:, 0] = state_0_shor; V9[:, 1] = state_1_shor
print(f"Shor [[9,1,3]]: built in {time.time()-t0:.2f}s")


# ============================================================
# Build encoded Bell state: |Phi+> = (|0>_ref|0_L> + |1>_ref|1_L>) / sqrt(2)
# Total: n+1 qubits. Qubit 0 = reference, qubits 1..n = code.
# ============================================================

def build_encoded_bell(V, n_code):
    """Build (|0>|0_L> + |1>|1_L>) / sqrt(2) as (n_code+1)-qubit state."""
    dim_code = 2**n_code
    # |0>|0_L>
    state_00 = np.zeros(2 * dim_code, dtype=complex)
    state_00[:dim_code] = V[:, 0]  # |0>_ref tensor |0_L>
    # |1>|1_L>
    state_11 = np.zeros(2 * dim_code, dtype=complex)
    state_11[dim_code:] = V[:, 1]  # |1>_ref tensor |1_L>
    psi = (state_00 + state_11) / np.sqrt(2)
    return np.outer(psi, psi.conj())


# ============================================================
# Page curve: MI(reference : k physical qubits) for all k
# For each k, average over all C(n,k) subsets (or sample if too many)
# ============================================================

def page_curve(rho_full, n_total, n_code, max_subsets=50):
    """Compute MI(ref : subset) for subsets of k=1..n_code physical qubits.

    Qubit 0 is reference, qubits 1..n_code are physical.
    For each k, we keep ref (qubit 0) + k physical qubits.
    MI(ref : phys_subset) = S(ref) + S(phys_subset) - S(ref + phys_subset)
    """
    results = {}

    for k in range(1, n_code + 1):
        phys_qubits = list(range(1, n_code + 1))
        all_subsets = list(combinations(phys_qubits, k))

        # Sample if too many
        if len(all_subsets) > max_subsets:
            rng = np.random.default_rng(42)
            indices = rng.choice(len(all_subsets), max_subsets, replace=False)
            subsets = [all_subsets[i] for i in indices]
        else:
            subsets = all_subsets

        mi_values = []
        for subset in subsets:
            keep_all = [0] + list(subset)  # ref + physical subset
            keep_phys = list(subset)       # physical subset only

            rho_ref_phys = partial_trace(rho_full, keep_all, n_total)
            rho_phys = partial_trace(rho_full, keep_phys, n_total)
            rho_ref = partial_trace(rho_full, [0], n_total)

            S_ref = von_neumann_entropy(rho_ref)
            S_phys = von_neumann_entropy(rho_phys)
            S_joint = von_neumann_entropy(rho_ref_phys)

            mi = S_ref + S_phys - S_joint
            mi_values.append(round(float(mi), 6))

        results[k] = {
            'mean_MI': round(float(np.mean(mi_values)), 6),
            'min_MI': round(float(np.min(mi_values)), 6),
            'max_MI': round(float(np.max(mi_values)), 6),
            'std_MI': round(float(np.std(mi_values)), 6),
            'n_subsets': len(subsets),
        }
        print(f"    k={k}: MI = {results[k]['mean_MI']:.4f} "
              f"(min={results[k]['min_MI']:.4f}, max={results[k]['max_MI']:.4f}, "
              f"n={len(subsets)} subsets)")

    return results


# ============================================================
# Timing test with [[5,1,3]] (smallest)
# ============================================================
print("\n=== Timing test: [[5,1,3]] Page curve ===")
t0 = time.time()
rho5 = build_encoded_bell(V5, 5)
print(f"Bell state built: {time.time()-t0:.2f}s, shape={rho5.shape}")

t0 = time.time()
# Test single partial trace at largest size
rho_test = partial_trace(rho5, [0, 1, 2], 6)
print(f"Single partial_trace (6 qubits, keep 3): {time.time()-t0:.2f}s")


# ============================================================
# Run Page curves for all three codes
# ============================================================
all_results = {}

for code_name, V, n_code in [
    ('five_qubit', V5, 5),
    ('steane', V7, 7),
    ('shor', V9, 9),
]:
    n_total = n_code + 1
    print(f"\n=== {code_name} [[{n_code},1,3]]: Page curve (n_total={n_total}) ===")

    # Check feasibility: n_total for Shor = 10 (at the limit!)
    if n_total > 10:
        print(f"  SKIPPING: n_total={n_total} exceeds 10-qubit limit")
        continue

    t0 = time.time()
    rho = build_encoded_bell(V, n_code)
    print(f"  Encoded Bell state: {time.time()-t0:.2f}s")

    t0 = time.time()
    pc = page_curve(rho, n_total, n_code, max_subsets=50)
    elapsed = time.time() - t0
    print(f"  Page curve computed in {elapsed:.1f}s")

    # Singleton bound: recovery requires k >= n - d + 1 = n - 2
    singleton_k = n_code - 2

    all_results[code_name] = {
        'n_code': n_code,
        'n_total': n_total,
        'distance': 3,
        'singleton_threshold': singleton_k,
        'page_curve': pc,
        'time_s': round(elapsed, 2),
    }


# ============================================================
# Summary comparison
# ============================================================
print("\n\n=== PAGE CURVE COMPARISON ===")
print(f"{'k':<5}", end='')
for code_name in ['five_qubit', 'steane', 'shor']:
    if code_name in all_results:
        print(f"{'MI_'+code_name:<22}", end='')
print()
print("-" * 71)

max_k = max(r['n_code'] for r in all_results.values())
for k in range(1, max_k + 1):
    print(f"{k:<5}", end='')
    for code_name in ['five_qubit', 'steane', 'shor']:
        if code_name in all_results and k <= all_results[code_name]['n_code']:
            pc = all_results[code_name]['page_curve']
            mi = pc[k]['mean_MI']
            spread = pc[k]['max_MI'] - pc[k]['min_MI']
            print(f"{mi:.4f} (±{spread:.3f})       ", end='')
        else:
            print(f"{'—':<22}", end='')
    print()

print("\nSingleton bound predictions (k >= n-d+1 for full recovery):")
for code_name in ['five_qubit', 'steane', 'shor']:
    if code_name in all_results:
        r = all_results[code_name]
        print(f"  {code_name}: k >= {r['singleton_threshold']} "
              f"(actual MI at k={r['singleton_threshold']}: "
              f"{r['page_curve'][r['singleton_threshold']]['mean_MI']:.4f})")


# Save
with open('results/sprint_015c_code_families_page_curves.json', 'w') as f:
    json.dump(all_results, f, indent=2)
print("\nResults saved to results/sprint_015c_code_families_page_curves.json")
