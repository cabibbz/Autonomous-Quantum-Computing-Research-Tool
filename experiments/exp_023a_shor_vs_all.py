"""
Sprint 023a: Shor [[9,1,3]] vs all codes across (gamma, lambda) noise landscape

The Shor code is a two-level concatenation: inner phase-flip (within each 3-qubit block)
and outer bit-flip (across 3 blocks). Does this structure create a specialization zone
that single-level codes can't access?

Compare basis-averaged Holevo for: uncoded, 3q-bitflip, 3q-phaseflip, [[5,1,3]], Shor [[9,1,3]]
across an 8x8 (gamma, lambda) grid.
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

# ============================================================
# Build code states
# ============================================================

# --- Uncoded ---
sv_0z_unc = np.array([1.0, 0.0])
sv_1z_unc = np.array([0.0, 1.0])

# --- 3-qubit bit-flip ---
n3 = 3; dim3 = 2**n3
sv_0z_bf = np.zeros(dim3); sv_0z_bf[0] = 1.0
sv_1z_bf = np.zeros(dim3); sv_1z_bf[7] = 1.0

# --- 3-qubit phase-flip ---
plus = np.array([1, 1]) / np.sqrt(2)
minus = np.array([1, -1]) / np.sqrt(2)
sv_0z_pf = np.kron(np.kron(plus, plus), plus)
sv_1z_pf = np.kron(np.kron(minus, minus), minus)

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

# --- Shor [[9,1,3]] code ---
# |0>_L = (|000>+|111>)^{⊗3} / 2√2
# |1>_L = (|000>-|111>)^{⊗3} / 2√2
n9 = 9; dim9 = 2**n9

block_plus = np.zeros(8)
block_plus[0b000] = 1.0; block_plus[0b111] = 1.0
block_plus /= np.linalg.norm(block_plus)  # (|000>+|111>)/sqrt(2)

block_minus = np.zeros(8)
block_minus[0b000] = 1.0; block_minus[0b111] = -1.0
block_minus /= np.linalg.norm(block_minus)  # (|000>-|111>)/sqrt(2)

sv_0z_shor = np.kron(np.kron(block_plus, block_plus), block_plus)  # already normalized
sv_1z_shor = np.kron(np.kron(block_minus, block_minus), block_minus)

# Verify orthogonality
assert abs(np.dot(sv_0z_shor.conj(), sv_1z_shor)) < 1e-10, "Shor logical states not orthogonal!"

print(f"Code states built in {time.time()-start:.1f}s")

# Helper to make basis states
def make_bases(sv0z, sv1z):
    sv0x = (sv0z + sv1z) / np.sqrt(2)
    sv1x = (sv0z - sv1z) / np.sqrt(2)
    sv0y = (sv0z + 1j * sv1z) / np.sqrt(2)
    sv1y = (sv0z - 1j * sv1z) / np.sqrt(2)
    return {'Z': (sv0z, sv1z), 'X': (sv0x, sv1x), 'Y': (sv0y, sv1y)}

codes = {
    'uncoded':      {'n': 1,  'bases': make_bases(sv_0z_unc, sv_1z_unc)},
    '3q-bitflip':   {'n': n3, 'bases': make_bases(sv_0z_bf, sv_1z_bf)},
    '3q-phaseflip': {'n': n3, 'bases': make_bases(sv_0z_pf, sv_1z_pf)},
    '[[5,1,3]]':    {'n': n5, 'bases': make_bases(sv_0z_513, sv_1z_513)},
    'shor-9':       {'n': n9, 'bases': make_bases(sv_0z_shor, sv_1z_shor)},
}

# ============================================================
# Time test: single Shor noise application
# ============================================================
t0 = time.time()
rho_test = density_matrix(sv_0z_shor)
rho_test_noisy = apply_combined_noise(rho_test, 0.1, 0.1, n9)
t_single = time.time() - t0
print(f"Single Shor noise application: {t_single:.2f}s")
# 8x8 grid × 3 bases × 2 states × 5 codes (but Shor is slowest)
# Shor: 64 × 3 × 2 = 384 noise applications
est_total = t_single * 384
print(f"Estimated total for Shor: {est_total:.0f}s")

if est_total > 50:
    print(f"WARNING: Too slow for 8x8. Reducing to 6x6 grid.")
    grid_size = 6
elif est_total > 30:
    print(f"Reducing to 6x6 grid for safety.")
    grid_size = 6
else:
    grid_size = 8

# ============================================================
# Sweep (gamma, lambda) grid
# ============================================================
noise_vals = np.linspace(0.05, 0.40, grid_size)

print(f"\n=== {grid_size}x{grid_size} Grid: Basis-Averaged Holevo ===\n")

grid_results = {}
winner_grid = {}

for gi, gamma in enumerate(noise_vals):
    for li, lam in enumerate(noise_vals):
        key = f"({gamma:.2f},{lam:.2f})"
        grid_results[key] = {}

        for code_name, code_data in codes.items():
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
            grid_results[key][code_name] = {
                'Z': round(holevos[0], 4),
                'X': round(holevos[1], 4),
                'Y': round(holevos[2], 4),
                'average': round(avg, 4),
                'asymmetry': round(asym, 4),
            }

        avgs = {name: grid_results[key][name]['average'] for name in codes}
        winner = max(avgs, key=avgs.get)
        winner_grid[key] = winner

elapsed_grid = time.time() - start
print(f"Grid computed in {elapsed_grid:.1f}s")

# ============================================================
# Analysis: winner map
# ============================================================
print(f"\n=== Winner Map ({grid_size}x{grid_size}) ===")
print(f"{'':>8s}", end="")
for lam in noise_vals:
    print(f" λ={lam:.2f}", end="")
print()

winner_counts = {}
for gi, gamma in enumerate(noise_vals):
    print(f"γ={gamma:.2f}", end="")
    for li, lam in enumerate(noise_vals):
        key = f"({gamma:.2f},{lam:.2f})"
        w = winner_grid[key]
        short = {'uncoded': 'UNC', '3q-bitflip': 'BF3', '3q-phaseflip': 'PF3',
                 '[[5,1,3]]': '513', 'shor-9': 'SH9'}[w]
        print(f"  {short:>5s}", end="")
        winner_counts[w] = winner_counts.get(w, 0) + 1
    print()

print(f"\nWinner counts: {winner_counts}")
total = sum(winner_counts.values())
for w, c in sorted(winner_counts.items(), key=lambda x: -x[1]):
    print(f"  {w}: {c}/{total} ({100*c/total:.0f}%)")

# ============================================================
# Key comparison: Shor vs [[5,1,3]]
# ============================================================
print(f"\n=== Shor vs [[5,1,3]]: Basis-Averaged Holevo Difference ===")
print(f"Positive = Shor better, Negative = [[5,1,3]] better\n")
print(f"{'':>8s}", end="")
for lam in noise_vals:
    print(f" λ={lam:.2f}", end="")
print()

shor_wins = 0
for gi, gamma in enumerate(noise_vals):
    print(f"γ={gamma:.2f}", end="")
    for li, lam in enumerate(noise_vals):
        key = f"({gamma:.2f},{lam:.2f})"
        diff = grid_results[key]['shor-9']['average'] - grid_results[key]['[[5,1,3]]']['average']
        sign = "+" if diff > 0 else ""
        print(f" {sign}{diff:6.3f}", end="")
        if diff > 0:
            shor_wins += 1
    print()

print(f"\nShor wins {shor_wins}/{total} grid points ({100*shor_wins/total:.0f}%)")

# ============================================================
# Asymmetry comparison
# ============================================================
print(f"\n=== Basis Asymmetry: Shor vs [[5,1,3]] ===")
print("Lower = more isotropic (better for unknown quantum states)\n")

shor_asyms = []
fivethree_asyms = []
for gi, gamma in enumerate(noise_vals):
    for li, lam in enumerate(noise_vals):
        key = f"({gamma:.2f},{lam:.2f})"
        shor_asyms.append(grid_results[key]['shor-9']['asymmetry'])
        fivethree_asyms.append(grid_results[key]['[[5,1,3]]']['asymmetry'])

print(f"Shor avg asymmetry:     {np.mean(shor_asyms):.4f} (max {np.max(shor_asyms):.4f})")
print(f"[[5,1,3]] avg asymmetry: {np.mean(fivethree_asyms):.4f} (max {np.max(fivethree_asyms):.4f})")

# ============================================================
# Per-qubit efficiency
# ============================================================
print(f"\n=== Per-Qubit Efficiency (avg Holevo / n_qubits) ===")
print("At balanced noise (γ=λ=middle of grid):\n")

mid = grid_size // 2
gamma_mid = noise_vals[mid]
lam_mid = noise_vals[mid]
key_mid = f"({gamma_mid:.2f},{lam_mid:.2f})"
for code_name in codes:
    n = codes[code_name]['n']
    avg = grid_results[key_mid][code_name]['average']
    eff = avg / n
    print(f"  {code_name:14s}: avg={avg:.3f}, n={n}, efficiency={eff:.4f}")

elapsed = time.time() - start
print(f"\nTotal time: {elapsed:.1f}s")

# Save results
results = {
    'sprint': '023a',
    'description': 'Shor [[9,1,3]] vs all codes across noise landscape - basis-averaged Holevo',
    'grid_size': grid_size,
    'noise_values': [round(float(v), 4) for v in noise_vals],
    'grid_results': grid_results,
    'winner_grid': winner_grid,
    'winner_counts': winner_counts,
    'shor_wins_vs_513': shor_wins,
    'shor_avg_asymmetry': round(float(np.mean(shor_asyms)), 4),
    'fivethree_avg_asymmetry': round(float(np.mean(fivethree_asyms)), 4),
    'elapsed_seconds': round(elapsed, 1)
}

with open('results/sprint_023a_shor_vs_all.json', 'w') as f:
    json.dump(results, f, indent=2)
print("Results saved to results/sprint_023a_shor_vs_all.json")
