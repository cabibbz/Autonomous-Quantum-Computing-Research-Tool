"""
Sprint 024b: Surface Code [[9,1,3]] vs Shor [[9,1,3]] Under Combined T1+T2 Noise

Head-to-head comparison: same n=9, same d=3, different architecture.
Which wins on basis-averaged and min-basis Holevo?
"""

import numpy as np
import json, time

start = time.time()

I2 = np.eye(2)
X = np.array([[0, 1], [1, 0]])
Z = np.array([[1, 0], [0, -1]])

def multi_pauli(ops_dict, n):
    result = np.eye(1)
    for q in range(n):
        result = np.kron(result, ops_dict.get(q, I2))
    return result

def density_matrix(sv):
    sv = sv.reshape(-1, 1)
    return sv @ sv.conj().T

def von_neumann_entropy(rho):
    evals = np.linalg.eigvalsh(rho)
    evals = evals[evals > 1e-12]
    return -np.sum(evals * np.log2(evals))

def apply_single_qubit_channel(rho, kraus_ops, qubit, n):
    """Apply single-qubit channel using tensor contraction (fast)."""
    rho_t = rho.astype(complex).reshape([2]*2*n)
    rho_new = np.zeros_like(rho_t, dtype=complex)
    for K in kraus_ops:
        contracted = np.tensordot(K, rho_t, axes=([1], [qubit]))
        contracted = np.moveaxis(contracted, 0, qubit)
        contracted = np.tensordot(K.conj(), contracted, axes=([1], [qubit + n]))
        contracted = np.moveaxis(contracted, 0, qubit + n)
        rho_new += contracted
    return rho_new.reshape(2**n, 2**n)

def apply_channel_all_qubits(rho, kraus_ops, n):
    """Apply identical single-qubit channel to all qubits."""
    for q in range(n):
        rho = apply_single_qubit_channel(rho, kraus_ops, q, n)
    return rho

def apply_combined_noise(rho, gamma, lam, n):
    # Amplitude damping
    K0_ad = np.array([[1, 0], [0, np.sqrt(1 - gamma)]], dtype=complex)
    K1_ad = np.array([[0, np.sqrt(gamma)], [0, 0]], dtype=complex)
    rho = apply_channel_all_qubits(rho, [K0_ad, K1_ad], n)
    # Phase damping
    K0_pd = np.array([[1, 0], [0, np.sqrt(1 - lam)]], dtype=complex)
    K1_pd = np.array([[0, 0], [0, np.sqrt(lam)]], dtype=complex)
    rho = apply_channel_all_qubits(rho, [K0_pd, K1_pd], n)
    return rho

def holevo_info(rho_0, rho_1):
    rho_avg = (rho_0 + rho_1) / 2
    return von_neumann_entropy(rho_avg) - (von_neumann_entropy(rho_0) + von_neumann_entropy(rho_1)) / 2

# ============================================================
# Build code states
# ============================================================
# Uncoded
sv_0z_unc = np.array([1.0, 0.0])
sv_1z_unc = np.array([0.0, 1.0])

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
state5_0 = np.zeros(dim5); state5_0[0] = 1.0
sv_0z_513 = proj_513 @ state5_0
sv_0z_513 /= np.linalg.norm(sv_0z_513)
X_L_513 = multi_pauli({q: X for q in range(5)}, 5)
sv_1z_513 = X_L_513 @ sv_0z_513
sv_1z_513 /= np.linalg.norm(sv_1z_513)

# Shor [[9,1,3]]
n9 = 9; dim9 = 2**n9
block_plus = np.zeros(8); block_plus[0] = 1.0; block_plus[7] = 1.0
block_plus /= np.linalg.norm(block_plus)
block_minus = np.zeros(8); block_minus[0] = 1.0; block_minus[7] = -1.0
block_minus /= np.linalg.norm(block_minus)
sv_0z_shor = np.kron(np.kron(block_plus, block_plus), block_plus)
sv_1z_shor = np.kron(np.kron(block_minus, block_minus), block_minus)

# Rotated Surface Code [[9,1,3]]
x_stabs = [{0: X, 1: X}, {1: X, 2: X, 4: X, 5: X}, {3: X, 4: X, 6: X, 7: X}, {7: X, 8: X}]
z_stabs = [{3: Z, 6: Z}, {0: Z, 1: Z, 3: Z, 4: Z}, {4: Z, 5: Z, 7: Z, 8: Z}, {2: Z, 5: Z}]
proj_surf = np.eye(dim9, dtype=complex)
for s in (x_stabs + z_stabs):
    S = multi_pauli(s, n9)
    proj_surf = proj_surf @ (np.eye(dim9) + S) / 2
state9_0 = np.zeros(dim9); state9_0[0] = 1.0
sv_0z_surf = proj_surf @ state9_0
sv_0z_surf /= np.linalg.norm(sv_0z_surf)
X_L_surf = multi_pauli({0: X, 3: X, 6: X}, n9)
Z_L_surf = multi_pauli({0: Z, 1: Z, 2: Z}, n9)
zl = np.real(sv_0z_surf.conj() @ Z_L_surf @ sv_0z_surf)
if zl < 0:
    sv_1z_surf = sv_0z_surf.copy()
    sv_0z_surf = X_L_surf @ sv_1z_surf
    sv_0z_surf /= np.linalg.norm(sv_0z_surf)
else:
    sv_1z_surf = X_L_surf @ sv_0z_surf
    sv_1z_surf /= np.linalg.norm(sv_1z_surf)

print(f"Code states built in {time.time()-start:.1f}s")

# Timing test
t0 = time.time()
rho_test = density_matrix(sv_0z_surf)
rho_test = apply_combined_noise(rho_test, 0.1, 0.1, n9)
t_single = time.time() - t0
print(f"Single 9q noise: {t_single:.3f}s")

# ============================================================
# Basis states
# ============================================================
def make_bases(sv0z, sv1z):
    sv0x = (sv0z + sv1z) / np.sqrt(2)
    sv1x = (sv0z - sv1z) / np.sqrt(2)
    sv0y = (sv0z + 1j * sv1z) / np.sqrt(2)
    sv1y = (sv0z - 1j * sv1z) / np.sqrt(2)
    return {'Z': (sv0z, sv1z), 'X': (sv0x, sv1x), 'Y': (sv0y, sv1y)}

codes = {
    'uncoded':    {'n': 1,  'bases': make_bases(sv_0z_unc, sv_1z_unc)},
    '[[5,1,3]]':  {'n': n5, 'bases': make_bases(sv_0z_513, sv_1z_513)},
    'shor-9':     {'n': n9, 'bases': make_bases(sv_0z_shor, sv_1z_shor)},
    'surface-9':  {'n': n9, 'bases': make_bases(sv_0z_surf, sv_1z_surf)},
}

# ============================================================
# Sweep 5x5 grid
# ============================================================
grid_size = 5
noise_vals = np.linspace(0.05, 0.35, grid_size)

print(f"\n=== {grid_size}x{grid_size} Grid ===\n")

grid_results = {}
winner_avg = {}
winner_min = {}

for gi, gamma in enumerate(noise_vals):
    for li, lam in enumerate(noise_vals):
        key = f"({gamma:.3f},{lam:.3f})"
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

            avg_h = float(np.mean(holevos))
            min_h = float(min(holevos))
            grid_results[key][code_name] = {
                'Z': round(holevos[0], 4),
                'X': round(holevos[1], 4),
                'Y': round(holevos[2], 4),
                'average': round(avg_h, 4),
                'min_basis': round(min_h, 4),
                'asymmetry': round(max(holevos) - min(holevos), 4),
            }

        avgs = {name: grid_results[key][name]['average'] for name in codes}
        winner_avg[key] = max(avgs, key=avgs.get)
        mins = {name: grid_results[key][name]['min_basis'] for name in codes}
        winner_min[key] = max(mins, key=mins.get)

    elapsed = time.time() - start
    print(f"  Row γ={gamma:.3f} done ({elapsed:.0f}s)")

# ============================================================
# Print results
# ============================================================
short = {'uncoded': 'UNC', '[[5,1,3]]': '513', 'shor-9': 'SH9', 'surface-9': 'SF9'}

print(f"\n=== Winner Map: Average Holevo ===")
count_avg = {}
for gi, gamma in enumerate(noise_vals):
    row = f"γ={gamma:.2f} "
    for li, lam in enumerate(noise_vals):
        key = f"({gamma:.3f},{lam:.3f})"
        w = winner_avg[key]
        row += f" {short[w]:>4s}"
        count_avg[w] = count_avg.get(w, 0) + 1
    print(row)

total = grid_size**2
print(f"\nAvg Holevo winners:")
for w, c in sorted(count_avg.items(), key=lambda x: -x[1]):
    print(f"  {w}: {c}/{total} ({100*c/total:.0f}%)")

print(f"\n=== Winner Map: Min-Basis Holevo ===")
count_min = {}
for gi, gamma in enumerate(noise_vals):
    row = f"γ={gamma:.2f} "
    for li, lam in enumerate(noise_vals):
        key = f"({gamma:.3f},{lam:.3f})"
        w = winner_min[key]
        row += f" {short[w]:>4s}"
        count_min[w] = count_min.get(w, 0) + 1
    print(row)

print(f"\nMin-basis Holevo winners:")
for w, c in sorted(count_min.items(), key=lambda x: -x[1]):
    print(f"  {w}: {c}/{total} ({100*c/total:.0f}%)")

# Head-to-head
print(f"\n=== Surface vs Shor (positive = Surface better) ===")
print("--- Avg Holevo ---")
s_avg_wins = 0
for gi, gamma in enumerate(noise_vals):
    row = f"γ={gamma:.2f} "
    for li, lam in enumerate(noise_vals):
        key = f"({gamma:.3f},{lam:.3f})"
        diff = grid_results[key]['surface-9']['average'] - grid_results[key]['shor-9']['average']
        row += f" {diff:+.3f}"
        if diff > 0: s_avg_wins += 1
    print(row)
print(f"Surface wins avg: {s_avg_wins}/{total}")

print("\n--- Min-Basis Holevo ---")
s_min_wins = 0
for gi, gamma in enumerate(noise_vals):
    row = f"γ={gamma:.2f} "
    for li, lam in enumerate(noise_vals):
        key = f"({gamma:.3f},{lam:.3f})"
        diff = grid_results[key]['surface-9']['min_basis'] - grid_results[key]['shor-9']['min_basis']
        row += f" {diff:+.3f}"
        if diff > 0: s_min_wins += 1
    print(row)
print(f"Surface wins min-basis: {s_min_wins}/{total}")

# Asymmetry
print(f"\n=== Basis Asymmetry ===")
for code_name in codes:
    asyms = [grid_results[key][code_name]['asymmetry'] for key in grid_results]
    print(f"  {code_name:12s}: mean={np.mean(asyms):.4f}, max={max(asyms):.4f}")

# Sample point detail
mid = grid_size // 2
gm = noise_vals[mid]; lm = noise_vals[mid]
key_mid = f"({gm:.3f},{lm:.3f})"
print(f"\n=== Detail at γ={gm:.3f}, λ={lm:.3f} ===")
for code_name in codes:
    d = grid_results[key_mid][code_name]
    print(f"  {code_name:12s}: Z={d['Z']:.4f} X={d['X']:.4f} Y={d['Y']:.4f} avg={d['average']:.4f} min={d['min_basis']:.4f} asym={d['asymmetry']:.4f}")

elapsed = time.time() - start
print(f"\nTotal time: {elapsed:.1f}s")

output = {
    'sprint': '024b',
    'description': 'Surface vs Shor [[9,1,3]] under combined T1+T2 noise',
    'grid_size': grid_size,
    'noise_values': [round(float(v), 4) for v in noise_vals],
    'grid_results': grid_results,
    'winner_avg': winner_avg,
    'winner_min': winner_min,
    'count_avg': count_avg,
    'count_min': count_min,
    'surface_avg_wins_vs_shor': s_avg_wins,
    'surface_min_wins_vs_shor': s_min_wins,
    'elapsed_seconds': round(elapsed, 1),
}

with open('results/sprint_024b_surface_vs_shor_noise.json', 'w') as f:
    json.dump(output, f, indent=2)
print("Results saved.")
