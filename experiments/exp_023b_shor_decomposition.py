"""
Sprint 023b: Decompose Shor's advantage — structure vs redundancy

Key question: Does Shor win because of its two-level concatenated STRUCTURE,
or just because 9 qubits > 5 qubits (more redundancy)?

Compare:
1. Shor [[9,1,3]] — two-level concatenation (phase inside bit)
2. [[5,1,3]] — isotropic, 5 qubits
3. 9-qubit bit-flip repetition — same qubit count as Shor, single-level
4. 9-qubit phase-flip repetition — same qubit count as Shor, single-level
5. Uncoded — baseline

Also analyze: Shor's per-basis breakdown to see which level handles which noise.
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

def make_bases(sv0z, sv1z):
    sv0x = (sv0z + sv1z) / np.sqrt(2)
    sv1x = (sv0z - sv1z) / np.sqrt(2)
    sv0y = (sv0z + 1j * sv1z) / np.sqrt(2)
    sv1y = (sv0z - 1j * sv1z) / np.sqrt(2)
    return {'Z': (sv0z, sv1z), 'X': (sv0x, sv1x), 'Y': (sv0y, sv1y)}

# ============================================================
# Build code states
# ============================================================

# --- Shor [[9,1,3]] ---
n9 = 9; dim9 = 2**n9
block_plus = np.zeros(8)
block_plus[0b000] = 1.0; block_plus[0b111] = 1.0
block_plus /= np.linalg.norm(block_plus)
block_minus = np.zeros(8)
block_minus[0b000] = 1.0; block_minus[0b111] = -1.0
block_minus /= np.linalg.norm(block_minus)
sv_0z_shor = np.kron(np.kron(block_plus, block_plus), block_plus)
sv_1z_shor = np.kron(np.kron(block_minus, block_minus), block_minus)

# --- 9-qubit bit-flip repetition ---
# |0>_L = |000000000>, |1>_L = |111111111>
sv_0z_9bf = np.zeros(dim9); sv_0z_9bf[0] = 1.0
sv_1z_9bf = np.zeros(dim9); sv_1z_9bf[dim9 - 1] = 1.0

# --- 9-qubit phase-flip repetition ---
# |0>_L = |+++++++++>, |1>_L = |--------->
plus = np.array([1, 1]) / np.sqrt(2)
minus = np.array([1, -1]) / np.sqrt(2)
sv_0z_9pf = plus.copy()
sv_1z_9pf = minus.copy()
for _ in range(8):
    sv_0z_9pf = np.kron(sv_0z_9pf, plus)
    sv_1z_9pf = np.kron(sv_1z_9pf, minus)

# --- [[5,1,3]] ---
n5 = 5; dim5 = 2**n5
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
sv_0z_513 = proj_513 @ state_0
sv_0z_513 /= np.linalg.norm(sv_0z_513)
X_L_513 = multi_pauli({q: X for q in range(5)}, 5)
sv_1z_513 = X_L_513 @ sv_0z_513
sv_1z_513 /= np.linalg.norm(sv_1z_513)

# --- Uncoded ---
sv_0z_unc = np.array([1.0, 0.0])
sv_1z_unc = np.array([0.0, 1.0])

print(f"Code states built in {time.time()-start:.1f}s")

# ============================================================
# Targeted comparison at key noise points
# ============================================================
# Select points where Shor won and lost in 23a
test_points = [
    (0.05, 0.05, "low balanced"),
    (0.19, 0.05, "high-γ low-λ (Shor best zone)"),
    (0.26, 0.19, "moderate (Shor wins)"),
    (0.05, 0.33, "low-γ high-λ ([[5,1,3]] best zone)"),
    (0.12, 0.26, "moderate dephasing"),
    (0.33, 0.33, "high balanced"),
]

codes_9q = {
    'shor-9':     {'n': n9, 'bases': make_bases(sv_0z_shor, sv_1z_shor)},
    '9q-bitflip': {'n': n9, 'bases': make_bases(sv_0z_9bf, sv_1z_9bf)},
    '9q-phaseflip': {'n': n9, 'bases': make_bases(sv_0z_9pf, sv_1z_9pf)},
    '[[5,1,3]]':  {'n': n5, 'bases': make_bases(sv_0z_513, sv_1z_513)},
    'uncoded':    {'n': 1,  'bases': make_bases(sv_0z_unc, sv_1z_unc)},
}

print("\n=== Structure vs Redundancy: 9-qubit codes ===\n")
print(f"{'noise point':>30s} | {'shor-9':>8s} {'9q-bf':>8s} {'9q-pf':>8s} {'[[5,1,3]]':>9s} {'uncoded':>8s} | winner")
print("-" * 95)

all_results = {}
for gamma, lam, label in test_points:
    key = f"({gamma:.2f},{lam:.2f})"
    all_results[key] = {'label': label}

    for code_name, code_data in codes_9q.items():
        n = code_data['n']
        holevos = []
        for basis in ['Z', 'X', 'Y']:
            sv0, sv1 = code_data['bases'][basis]
            rho_0 = density_matrix(sv0)
            rho_1 = density_matrix(sv1)
            rho_0n = apply_combined_noise(rho_0, gamma, lam, n)
            rho_1n = apply_combined_noise(rho_1, gamma, lam, n)
            h = holevo_info(rho_0n, rho_1n)
            holevos.append(float(h))

        avg = float(np.mean(holevos))
        asym = float(max(holevos) - min(holevos))
        all_results[key][code_name] = {
            'Z': round(holevos[0], 4),
            'X': round(holevos[1], 4),
            'Y': round(holevos[2], 4),
            'average': round(avg, 4),
            'asymmetry': round(asym, 4),
        }

    avgs = {name: all_results[key][name]['average'] for name in codes_9q}
    winner = max(avgs, key=avgs.get)
    all_results[key]['winner'] = winner

    print(f"  {label:>28s} | {avgs['shor-9']:8.3f} {avgs['9q-bitflip']:8.3f} {avgs['9q-phaseflip']:8.3f} {avgs['[[5,1,3]]']:9.3f} {avgs['uncoded']:8.3f} | {winner}")

# ============================================================
# Shor per-basis breakdown vs [[5,1,3]]
# ============================================================
print("\n\n=== Per-Basis Breakdown: Shor vs [[5,1,3]] ===\n")
for gamma, lam, label in test_points:
    key = f"({gamma:.2f},{lam:.2f})"
    print(f"\n{label} (γ={gamma}, λ={lam}):")
    for code_name in ['shor-9', '[[5,1,3]]']:
        r = all_results[key][code_name]
        print(f"  {code_name:14s}: Z={r['Z']:.3f}  X={r['X']:.3f}  Y={r['Y']:.3f}  avg={r['average']:.3f}  asym={r['asymmetry']:.3f}")

# ============================================================
# The critical test: does Shor beat 9q-repetition codes?
# If YES: structure matters. If NO: just redundancy.
# ============================================================
print("\n\n=== CRITICAL TEST: Shor vs Best 9-qubit Repetition ===")
print("If Shor > max(9q-bf, 9q-pf): concatenation STRUCTURE adds value")
print("If Shor <= max(9q-bf, 9q-pf): just redundancy, no structural advantage\n")

structure_wins = 0
for gamma, lam, label in test_points:
    key = f"({gamma:.2f},{lam:.2f})"
    shor_avg = all_results[key]['shor-9']['average']
    best_rep = max(all_results[key]['9q-bitflip']['average'],
                   all_results[key]['9q-phaseflip']['average'])
    diff = shor_avg - best_rep
    verdict = "STRUCTURE wins" if diff > 0.001 else "REDUNDANCY only" if diff > -0.001 else "REPETITION wins"
    if diff > 0.001:
        structure_wins += 1
    print(f"  {label:>28s}: Shor={shor_avg:.3f}, best_rep={best_rep:.3f}, diff={diff:+.3f} → {verdict}")

print(f"\nStructure wins: {structure_wins}/{len(test_points)} points")

# ============================================================
# Shor Z-Holevo vs X-Holevo ratio (measures inner vs outer balance)
# ============================================================
print("\n\n=== Shor Z/X Holevo Ratio (inner vs outer balance) ===")
print("Z=1, X<<1: outer (bit-flip) dominates. Z<<1, X=1: inner (phase-flip) dominates.\n")
for gamma, lam, label in test_points:
    key = f"({gamma:.2f},{lam:.2f})"
    r = all_results[key]['shor-9']
    ratio = r['Z'] / max(r['X'], 1e-6)
    print(f"  {label:>28s}: Z/X = {ratio:.2f}  (Z={r['Z']:.3f}, X={r['X']:.3f})")

elapsed = time.time() - start
print(f"\nTotal time: {elapsed:.1f}s")

# Save
results = {
    'sprint': '023b',
    'description': 'Decomposing Shor advantage: structure vs redundancy',
    'test_points': [(g, l, lab) for g, l, lab in test_points],
    'results': all_results,
    'structure_wins': structure_wins,
    'elapsed_seconds': round(elapsed, 1)
}

with open('results/sprint_023b_shor_decomposition.json', 'w') as f:
    json.dump(results, f, indent=2)
print("Results saved to results/sprint_023b_shor_decomposition.json")
