"""
Sprint 022c: WHY does isotropy always win? Information-theoretic analysis.

Three analyses:
1. Effective distance decomposition — how does each code protect against X, Y, Z errors?
2. Per-qubit efficiency — Holevo/n_physical. Does [[5,1,3]]'s overhead pay for itself?
3. Hybrid strategy — what if you allocate 5 qubits as "3 phase-flip + 2 uncoded"
   vs "5 isotropic"? Can combining specialized codes approach isotropy?
"""

import numpy as np
import json, time

start = time.time()

I2 = np.eye(2)
X = np.array([[0, 1], [1, 0]])
Y = np.array([[0, -1j], [1j, 0]])
Z = np.array([[1, 0], [0, -1]])

def density_matrix(sv):
    sv = sv.reshape(-1, 1)
    return sv @ sv.conj().T

def von_neumann_entropy(rho):
    evals = np.linalg.eigvalsh(rho)
    evals = evals[evals > 1e-12]
    return -np.sum(evals * np.log2(evals))

def multi_pauli(ops_dict, n):
    result = np.eye(1)
    for q in range(n):
        result = np.kron(result, ops_dict.get(q, I2))
    return result

def apply_single_error(rho, error_op, p, n):
    """Apply independent single-qubit errors of type error_op with probability p."""
    result = rho.copy().astype(complex)
    for q in range(n):
        K0_coeff = np.sqrt(1 - p)
        K1_coeff = np.sqrt(p)
        K0_full = K0_coeff * np.eye(1, dtype=complex)
        K1_full = K1_coeff * np.eye(1, dtype=complex)
        for i in range(n):
            K0_full = np.kron(K0_full, I2 if True else I2)
            K1_full = np.kron(K1_full, error_op if i == q else I2)
        # Simpler: just apply identity*(1-p) + error*p per qubit
        E = np.eye(1, dtype=complex)
        for i in range(n):
            E = np.kron(E, error_op if i == q else I2)
        result = (1 - p) * result + p * (E @ result @ E.conj().T)
    return result

def apply_amplitude_damping(rho, gamma, n):
    K0 = np.array([[1, 0], [0, np.sqrt(1 - gamma)]])
    K1 = np.array([[0, np.sqrt(gamma)], [0, 0]])
    result = rho.copy().astype(complex)
    for q in range(n):
        rho_new = np.zeros_like(result)
        for K in [K0, K1]:
            K_full = np.eye(1, dtype=complex)
            for i in range(n):
                K_full = np.kron(K_full, K if i == q else I2)
            rho_new += K_full @ result @ K_full.conj().T
        result = rho_new
    return result

def apply_phase_damping(rho, lam, n):
    K0 = np.array([[1, 0], [0, np.sqrt(1 - lam)]])
    K1 = np.array([[0, 0], [0, np.sqrt(lam)]])
    result = rho.copy().astype(complex)
    for q in range(n):
        rho_new = np.zeros_like(result)
        for K in [K0, K1]:
            K_full = np.eye(1, dtype=complex)
            for i in range(n):
                K_full = np.kron(K_full, K if i == q else I2)
            rho_new += K_full @ result @ K_full.conj().T
        result = rho_new
    return result

def apply_combined_noise(rho, gamma, lam, n):
    rho = apply_amplitude_damping(rho, gamma, n)
    rho = apply_phase_damping(rho, lam, n)
    return rho

def holevo_info(rho_0, rho_1):
    rho_avg = (rho_0 + rho_1) / 2
    return von_neumann_entropy(rho_avg) - (von_neumann_entropy(rho_0) + von_neumann_entropy(rho_1)) / 2

def make_bases(sv_0z, sv_1z):
    return {
        'Z': (sv_0z, sv_1z),
        'X': ((sv_0z + sv_1z) / np.sqrt(2), (sv_0z - sv_1z) / np.sqrt(2)),
        'Y': ((sv_0z + 1j * sv_1z) / np.sqrt(2), (sv_0z - 1j * sv_1z) / np.sqrt(2)),
    }

def basis_avg_holevo_combined(states, gamma, lam, n):
    total = 0.0
    per_basis = {}
    for basis in ['Z', 'X', 'Y']:
        sv0, sv1 = states[basis]
        rho_0 = density_matrix(sv0)
        rho_1 = density_matrix(sv1)
        rho_0n = apply_combined_noise(rho_0, gamma, lam, n)
        rho_1n = apply_combined_noise(rho_1, gamma, lam, n)
        h = holevo_info(rho_0n, rho_1n)
        per_basis[basis] = h
        total += h
    return total / 3.0, per_basis

# ============================================================
# Build code states
# ============================================================

n3 = 3; dim3 = 8
sv_0z_bf = np.zeros(dim3); sv_0z_bf[0] = 1.0
sv_1z_bf = np.zeros(dim3); sv_1z_bf[7] = 1.0

plus = np.array([1, 1]) / np.sqrt(2)
minus = np.array([1, -1]) / np.sqrt(2)
sv_0z_pf = np.kron(np.kron(plus, plus), plus)
sv_1z_pf = np.kron(np.kron(minus, minus), minus)

n5 = 5; dim5 = 32
stabs_513 = [
    {0: X, 1: Z, 2: Z, 3: X},
    {1: X, 2: Z, 3: Z, 4: X},
    {0: X, 2: X, 3: Z, 4: Z},
    {0: Z, 1: X, 3: X, 4: Z},
]
proj_513 = np.eye(dim5)
for s in stabs_513:
    S = multi_pauli(s, n5)
    proj_513 = proj_513 @ (np.eye(dim5) + S) / 2
state_0 = np.zeros(dim5); state_0[0] = 1.0
sv_0z_513 = proj_513 @ state_0; sv_0z_513 /= np.linalg.norm(sv_0z_513)
X_L_513 = multi_pauli({q: X for q in range(5)}, 5)
sv_1z_513 = X_L_513 @ sv_0z_513; sv_1z_513 /= np.linalg.norm(sv_1z_513)

sv_0z_unc = np.array([1.0, 0.0])
sv_1z_unc = np.array([0.0, 1.0])

codes = {
    'uncoded': {'n': 1, 'states': make_bases(sv_0z_unc, sv_1z_unc)},
    '3q-bitflip': {'n': n3, 'states': make_bases(sv_0z_bf, sv_1z_bf)},
    '3q-phaseflip': {'n': n3, 'states': make_bases(sv_0z_pf, sv_1z_pf)},
    '[[5,1,3]]': {'n': n5, 'states': make_bases(sv_0z_513, sv_1z_513)},
}

print(f"States built in {time.time()-start:.1f}s")

# ============================================================
# Analysis 1: Per-error-type Holevo under pure X, Y, Z errors
# ============================================================
print("\n=== Analysis 1: Per-Error-Type Protection (Pauli channel) ===\n")

error_probs = [0.05, 0.10, 0.15, 0.20]
pauli_ops = {'X': X, 'Y': Y, 'Z': Z}

analysis1 = {}
for p in error_probs:
    analysis1[f"p={p}"] = {}
    print(f"\nError probability p={p}:")
    print(f"  {'Code':>14s} | {'X-err':>8s} {'Y-err':>8s} {'Z-err':>8s} | {'Min':>6s} {'Max':>6s} {'Balance':>8s}")
    print("  " + "-" * 70)
    for code_name, code_data in codes.items():
        n = code_data['n']
        err_holevo = {}
        for err_name, err_op in pauli_ops.items():
            # Apply single-type errors, measure basis-averaged Holevo
            total = 0.0
            for basis in ['Z', 'X', 'Y']:
                sv0, sv1 = code_data['states'][basis]
                rho_0 = density_matrix(sv0)
                rho_1 = density_matrix(sv1)
                rho_0n = apply_single_error(rho_0, err_op, p, n)
                rho_1n = apply_single_error(rho_1, err_op, p, n)
                total += holevo_info(rho_0n, rho_1n)
            err_holevo[err_name] = total / 3.0

        vals = list(err_holevo.values())
        balance = 1.0 - (max(vals) - min(vals)) / max(max(vals), 1e-10)
        analysis1[f"p={p}"][code_name] = {k: round(float(v), 4) for k, v in err_holevo.items()}
        analysis1[f"p={p}"][code_name]['balance'] = round(float(balance), 4)
        print(f"  {code_name:>14s} | {err_holevo['X']:8.3f} {err_holevo['Y']:8.3f} {err_holevo['Z']:8.3f} | {min(vals):6.3f} {max(vals):6.3f} {balance:8.3f}")

# ============================================================
# Analysis 2: Per-qubit efficiency under combined noise
# ============================================================
print("\n\n=== Analysis 2: Per-Qubit Efficiency (Holevo/n_physical) ===\n")

test_points = [(0.05, 0.05), (0.05, 0.20), (0.20, 0.05), (0.10, 0.10), (0.20, 0.20)]
analysis2 = {}

print(f"  {'(γ,λ)':>12s} | {'Code':>14s} {'n':>3s} {'Holevo':>8s} {'H/n':>8s}")
print("  " + "-" * 55)
for gamma, lam in test_points:
    key = f"({gamma},{lam})"
    analysis2[key] = {}
    for code_name, code_data in codes.items():
        n = code_data['n']
        avg, _ = basis_avg_holevo_combined(code_data['states'], gamma, lam, n)
        eff = avg / n
        analysis2[key][code_name] = {'holevo': round(float(avg), 4), 'n': n, 'efficiency': round(float(eff), 4)}
        print(f"  ({gamma:.2f},{lam:.2f}) | {code_name:>14s} {n:3d} {avg:8.3f} {eff:8.4f}")
    print()

# ============================================================
# Analysis 3: Hybrid strategy — 5 qubits as "3-pf + 2 uncoded"
# ============================================================
print("\n=== Analysis 3: Hybrid Strategy (5-qubit budget) ===\n")
print("Compare: [[5,1,3]] (1 logical qubit, isotropic)")
print("     vs: 3q-phaseflip (1 logical) + 1 uncoded (1 logical) = 2 logical qubits")
print("     vs: 3q-bitflip (1 logical) + 1 uncoded (1 logical) = 2 logical qubits")
print("     vs: 5 uncoded qubits = 5 logical qubits")
print()
print("Metric: total logical information preserved (sum of Holevo over logical qubits)\n")

analysis3 = {}
test_noise = [(0.05, 0.05), (0.05, 0.20), (0.05, 0.40), (0.10, 0.10), (0.20, 0.20)]

for gamma, lam in test_noise:
    key = f"({gamma},{lam})"

    # [[5,1,3]]: 1 logical qubit
    h_513, _ = basis_avg_holevo_combined(codes['[[5,1,3]]']['states'], gamma, lam, 5)

    # 3q-phaseflip + 1 uncoded: 2 logical qubits (3+1=4 physical, but let's use 5-qubit budget)
    h_pf, _ = basis_avg_holevo_combined(codes['3q-phaseflip']['states'], gamma, lam, 3)
    h_unc, _ = basis_avg_holevo_combined(codes['uncoded']['states'], gamma, lam, 1)
    h_hybrid_pf = h_pf + h_unc  # 2 logical qubits, 4 physical

    # 3q-bitflip + 1 uncoded
    h_bf, _ = basis_avg_holevo_combined(codes['3q-bitflip']['states'], gamma, lam, 3)
    h_hybrid_bf = h_bf + h_unc  # 2 logical qubits, 4 physical

    # Best hybrid: pick whichever 3q code is better
    h_hybrid_best = max(h_hybrid_pf, h_hybrid_bf)
    hybrid_name = 'pf+unc' if h_hybrid_pf >= h_hybrid_bf else 'bf+unc'

    # 5 uncoded: 5 logical qubits
    h_5unc = 5 * h_unc

    analysis3[key] = {
        '[[5,1,3]]_1logical': round(float(h_513), 4),
        'hybrid_2logical': round(float(h_hybrid_best), 4),
        'hybrid_type': hybrid_name,
        '5uncoded_5logical': round(float(h_5unc), 4),
        'uncoded_per_qubit': round(float(h_unc), 4),
        '513_per_qubit': round(float(h_513), 4),  # 1 logical / 5 physical
    }

    print(f"(γ={gamma:.2f}, λ={lam:.2f}):")
    print(f"  [[5,1,3]]:       {h_513:.3f} Holevo (1 logical qubit, 5 physical)")
    print(f"  Best hybrid:     {h_hybrid_best:.3f} Holevo ({hybrid_name}, 2 logical, 4 physical)")
    print(f"  5×uncoded:       {h_5unc:.3f} Holevo (5 logical, 5 physical)")
    print(f"  Per-physical-qubit: [[5,1,3]]={h_513/5:.3f}, hybrid={h_hybrid_best/4:.3f}, uncoded={h_unc:.3f}")
    # Quality vs quantity tradeoff
    print(f"  Quality ([[5,1,3]]/uncoded per logical): {h_513/h_unc:.2f}x")
    print(f"  Quantity (5×uncoded/[[5,1,3]] total): {h_5unc/h_513:.2f}x")
    print()

# ============================================================
# Analysis 4: The isotropy theorem — distance decomposition
# ============================================================
print("\n=== Analysis 4: Distance Decomposition ===\n")
print("Effective distance against each Pauli type:")
print("  Code             | d_X  d_Y  d_Z | Min(d) = d_code")
print("  " + "-" * 50)

# For repetition codes:
# bit-flip: d_X = 3, d_Y = 1, d_Z = 1 (detects only X errors)
# phase-flip: d_X = 1, d_Y = 1, d_Z = 3 (detects only Z errors)
# [[5,1,3]]: d_X = 3, d_Y = 3, d_Z = 3 (detects all)

distance_data = {
    '3q-bitflip': {'d_X': 3, 'd_Y': 1, 'd_Z': 1, 'min_d': 1},
    '3q-phaseflip': {'d_X': 1, 'd_Y': 1, 'd_Z': 3, 'min_d': 1},
    '[[5,1,3]]': {'d_X': 3, 'd_Y': 3, 'd_Z': 3, 'min_d': 3},
}

for code_name, dd in distance_data.items():
    print(f"  {code_name:>14s}   |  {dd['d_X']}    {dd['d_Y']}    {dd['d_Z']}  |    {dd['min_d']}")

print("\nKey insight: The CODE DISTANCE is min(d_X, d_Y, d_Z).")
print("3-qubit codes have distance 1 against their unprotected errors,")
print("making them effectively distance-1 codes overall.")
print("[[5,1,3]] with distance 3 against ALL error types is genuinely distance-3.")
print("\nIsotropy wins because: the weakest direction determines the code distance,")
print("and unknown quantum states have components in ALL directions.")

elapsed = time.time() - start
print(f"\nTotal time: {elapsed:.1f}s")

results = {
    'sprint': '022c',
    'description': 'Why isotropy always wins: distance decomposition and efficiency analysis',
    'analysis1_per_error': analysis1,
    'analysis2_efficiency': analysis2,
    'analysis3_hybrid': analysis3,
    'analysis4_distances': distance_data,
    'elapsed_seconds': round(elapsed, 1)
}

with open('results/sprint_022c_why_isotropy_wins.json', 'w') as f:
    json.dump(results, f, indent=2)
print("Results saved to results/sprint_022c_why_isotropy_wins.json")
