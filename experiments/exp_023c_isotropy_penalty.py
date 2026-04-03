"""
Sprint 023c: Isotropy-adjusted comparison — does Shor really win?

23a showed Shor wins 67% of grid on basis-averaged Holevo.
23b showed Shor has extreme X-axis specialization (near-perfect X, mediocre Z/Y).

Key questions:
1. Using the min-basis Holevo (worst-case) instead of average — who wins?
   (If your code can't protect Z info, it can't do arbitrary quantum computation)
2. Per-qubit efficiency: Shor uses 9 qubits vs [[5,1,3]]'s 5. Fair comparison?
3. Code Quality Score = min-basis Holevo × isotropy — which code architecture is best?

Also test: "reversed Shor" = bit-flip inside phase-flip (opposite concatenation order).
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

# Shor [[9,1,3]] (phase inside bit)
n9 = 9; dim9 = 2**n9
block_plus = np.zeros(8)
block_plus[0b000] = 1.0; block_plus[0b111] = 1.0
block_plus /= np.linalg.norm(block_plus)
block_minus = np.zeros(8)
block_minus[0b000] = 1.0; block_minus[0b111] = -1.0
block_minus /= np.linalg.norm(block_minus)
sv_0z_shor = np.kron(np.kron(block_plus, block_plus), block_plus)
sv_1z_shor = np.kron(np.kron(block_minus, block_minus), block_minus)

# "Reversed Shor" = bit-flip inside phase-flip
# Inner: 3-qubit bit-flip: |0>_L=|000>, |1>_L=|111>
# Outer: 3-qubit phase-flip on encoded blocks
# |0>_L_rev = |block_0>|block_0>|block_0> + |block_1>|block_1>|block_1> (in Hadamard frame)
# Actually: outer phase-flip means |0>_L = |+_L>|+_L>|+_L>, |1>_L = |-_L>|-_L>|-_L>
# where |+_L> = (|000>+|111>)/√2, |-_L> = (|000>-|111>)/√2
# But that's the SAME as Shor! The Shor code IS symmetric under inner/outer swap
# because (phase-flip inside bit-flip) = (bit-flip inside phase-flip) for 3-qubit codes.

# Let me verify this is NOT the same by building it differently.
# Reversed: inner = bit-flip (|0>_I=|000>, |1>_I=|111>)
# Outer = phase-flip: |0>_L = |+_I>|+_I>|+_I>, |1>_L = |-_I>|-_I>|-_I>
# |+_I> = (|000>+|111>)/sqrt(2) = block_plus
# |-_I> = (|000>-|111>)/sqrt(2) = block_minus
# So |0>_L = block_plus ⊗ block_plus ⊗ block_plus = SAME as Shor!
# Confirmed: the Shor code is self-dual under inner/outer swap.
print("Note: 'reversed Shor' = same as Shor (self-dual under inner/outer swap)")

# [[5,1,3]]
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

# Uncoded
sv_0z_unc = np.array([1.0, 0.0])
sv_1z_unc = np.array([0.0, 1.0])

print(f"Code states built in {time.time()-start:.1f}s")

codes = {
    'uncoded':   {'n': 1,  'bases': make_bases(sv_0z_unc, sv_1z_unc)},
    '[[5,1,3]]': {'n': n5, 'bases': make_bases(sv_0z_513, sv_1z_513)},
    'shor-9':    {'n': n9, 'bases': make_bases(sv_0z_shor, sv_1z_shor)},
}

# ============================================================
# 6x6 grid: compute all metrics
# ============================================================
grid_size = 6
noise_vals = np.linspace(0.05, 0.40, grid_size)

print(f"\n=== {grid_size}x{grid_size} Grid: Multiple Figures of Merit ===\n")

grid_results = {}
for gi, gamma in enumerate(noise_vals):
    for li, lam in enumerate(noise_vals):
        key = f"({gamma:.2f},{lam:.2f})"
        grid_results[key] = {}

        for code_name, code_data in codes.items():
            n = code_data['n']
            holevos = {}
            for basis in ['Z', 'X', 'Y']:
                sv0, sv1 = code_data['bases'][basis]
                rho_0 = density_matrix(sv0)
                rho_1 = density_matrix(sv1)
                rho_0n = apply_combined_noise(rho_0, gamma, lam, n)
                rho_1n = apply_combined_noise(rho_1, gamma, lam, n)
                holevos[basis] = float(holevo_info(rho_0n, rho_1n))

            avg = np.mean(list(holevos.values()))
            min_h = min(holevos.values())
            max_h = max(holevos.values())
            asym = max_h - min_h
            isotropy = 1.0 - asym  # 1 = perfect isotropy, 0 = worst
            quality = min_h  # worst-case basis = actual quantum computational ability

            grid_results[key][code_name] = {
                'Z': round(holevos['Z'], 4),
                'X': round(holevos['X'], 4),
                'Y': round(holevos['Y'], 4),
                'average': round(float(avg), 4),
                'min_basis': round(float(min_h), 4),
                'max_basis': round(float(max_h), 4),
                'asymmetry': round(float(asym), 4),
                'isotropy': round(float(isotropy), 4),
                'quality': round(float(quality), 4),
                'per_qubit_avg': round(float(avg / n), 4),
                'per_qubit_min': round(float(min_h / n), 4),
            }

elapsed_grid = time.time() - start
print(f"Grid computed in {elapsed_grid:.1f}s")

# ============================================================
# Winner maps under different metrics
# ============================================================
metrics = {
    'average': 'Basis-Averaged Holevo (favors Shor)',
    'min_basis': 'Min-Basis Holevo (worst-case, favors isotropy)',
    'quality': 'Code Quality = min-basis (quantum computation capability)',
    'per_qubit_avg': 'Per-Qubit Efficiency (avg/n)',
    'per_qubit_min': 'Per-Qubit Min-Basis Efficiency (min/n)',
}

for metric, label in metrics.items():
    print(f"\n=== Winner Map: {label} ===")
    print(f"{'':>8s}", end="")
    for lam in noise_vals:
        print(f" λ={lam:.2f}", end="")
    print()

    counts = {}
    for gi, gamma in enumerate(noise_vals):
        print(f"γ={gamma:.2f}", end="")
        for li, lam in enumerate(noise_vals):
            key = f"({gamma:.2f},{lam:.2f})"
            vals = {name: grid_results[key][name][metric] for name in codes}
            winner = max(vals, key=vals.get)
            short = {'uncoded': 'UNC', '[[5,1,3]]': '513', 'shor-9': 'SH9'}[winner]
            print(f"  {short:>5s}", end="")
            counts[winner] = counts.get(winner, 0) + 1
        print()

    total = sum(counts.values())
    print(f"  Counts: ", end="")
    for w, c in sorted(counts.items(), key=lambda x: -x[1]):
        print(f"{w}={c}/{total}({100*c/total:.0f}%) ", end="")
    print()

# ============================================================
# Direct comparison table at key points
# ============================================================
print("\n\n=== Head-to-Head: Shor vs [[5,1,3]] at Grid Corners ===\n")
corners = [(0, 0), (0, -1), (-1, 0), (-1, -1), (grid_size//2, grid_size//2)]
corner_labels = ['low-γ low-λ', 'low-γ high-λ', 'high-γ low-λ', 'high-γ high-λ', 'balanced']

print(f"{'point':>18s} | {'metric':>12s} | {'[[5,1,3]]':>10s} {'shor-9':>10s} {'diff':>8s} | better")
print("-" * 80)

for (gi, li), label in zip(corners, corner_labels):
    gamma = noise_vals[gi]; lam = noise_vals[li]
    key = f"({gamma:.2f},{lam:.2f})"
    for metric in ['average', 'min_basis', 'per_qubit_min']:
        v513 = grid_results[key]['[[5,1,3]]'][metric]
        vsh = grid_results[key]['shor-9'][metric]
        diff = vsh - v513
        better = 'Shor' if diff > 0.001 else '[[5,1,3]]' if diff < -0.001 else 'TIE'
        print(f"  {label:>16s} | {metric:>12s} | {v513:10.4f} {vsh:10.4f} {diff:+8.4f} | {better}")
    print()

# ============================================================
# Summary statistics
# ============================================================
print("\n=== SUMMARY: Which metric, which winner? ===\n")
summary = {}
for metric in ['average', 'min_basis', 'per_qubit_avg', 'per_qubit_min']:
    wins = {'uncoded': 0, '[[5,1,3]]': 0, 'shor-9': 0}
    for gi, gamma in enumerate(noise_vals):
        for li, lam in enumerate(noise_vals):
            key = f"({gamma:.2f},{lam:.2f})"
            vals = {name: grid_results[key][name][metric] for name in codes}
            winner = max(vals, key=vals.get)
            wins[winner] += 1
    summary[metric] = wins
    total = sum(wins.values())
    print(f"{metric:>16s}: ", end="")
    for w in ['[[5,1,3]]', 'shor-9', 'uncoded']:
        print(f"{w}={wins[w]}/{total}({100*wins[w]/total:.0f}%) ", end="")
    print()

elapsed = time.time() - start
print(f"\nTotal time: {elapsed:.1f}s")

# Save
results = {
    'sprint': '023c',
    'description': 'Isotropy-adjusted comparison: which figure of merit makes Shor vs [[5,1,3]] decision?',
    'grid_size': grid_size,
    'noise_values': [round(float(v), 4) for v in noise_vals],
    'grid_results': grid_results,
    'summary_wins': summary,
    'elapsed_seconds': round(elapsed, 1)
}

with open('results/sprint_023c_isotropy_penalty.json', 'w') as f:
    json.dump(results, f, indent=2)
print("Results saved to results/sprint_023c_isotropy_penalty.json")
