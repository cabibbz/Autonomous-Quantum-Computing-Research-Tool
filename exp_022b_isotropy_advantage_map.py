"""
Sprint 022b: Isotropy advantage map across full (γ, λ) noise landscape

22a showed [[5,1,3]] wins everywhere. This experiment maps:
1. The [[5,1,3]] advantage over the BEST specialized code at each point
2. A "bias-aware oracle" that picks the optimal 3-qubit code for each noise point
3. Whether combining both specialized codes (bit-flip AND phase-flip) closes the gap

Key metric: basis-averaged Holevo information
"""

import numpy as np
import json, time

start = time.time()

I2 = np.eye(2)
X = np.array([[0, 1], [1, 0]])
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

def basis_avg_holevo(code_states, gamma, lam, n):
    """Compute basis-averaged Holevo for a code under combined noise."""
    total = 0.0
    for basis in ['Z', 'X', 'Y']:
        sv0, sv1 = code_states[basis]
        rho_0 = density_matrix(sv0)
        rho_1 = density_matrix(sv1)
        rho_0n = apply_combined_noise(rho_0, gamma, lam, n)
        rho_1n = apply_combined_noise(rho_1, gamma, lam, n)
        total += holevo_info(rho_0n, rho_1n)
    return total / 3.0

# ============================================================
# Build code states
# ============================================================

# --- 3-qubit bit-flip ---
n3 = 3; dim3 = 8
sv_0z_bf = np.zeros(dim3); sv_0z_bf[0] = 1.0
sv_1z_bf = np.zeros(dim3); sv_1z_bf[7] = 1.0

# --- 3-qubit phase-flip ---
plus = np.array([1, 1]) / np.sqrt(2)
minus = np.array([1, -1]) / np.sqrt(2)
sv_0z_pf = np.kron(np.kron(plus, plus), plus)
sv_1z_pf = np.kron(np.kron(minus, minus), minus)

# --- [[5,1,3]] ---
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
sv_0z_513 = proj_513 @ state_0
sv_0z_513 /= np.linalg.norm(sv_0z_513)
X_L_513 = multi_pauli({q: X for q in range(5)}, 5)
sv_1z_513 = X_L_513 @ sv_0z_513
sv_1z_513 /= np.linalg.norm(sv_1z_513)

# --- Uncoded ---
sv_0z_unc = np.array([1.0, 0.0])
sv_1z_unc = np.array([0.0, 1.0])

def make_bases(sv_0z, sv_1z):
    return {
        'Z': (sv_0z, sv_1z),
        'X': ((sv_0z + sv_1z) / np.sqrt(2), (sv_0z - sv_1z) / np.sqrt(2)),
        'Y': ((sv_0z + 1j * sv_1z) / np.sqrt(2), (sv_0z - 1j * sv_1z) / np.sqrt(2)),
    }

codes = {
    'uncoded': {'n': 1, 'states': make_bases(sv_0z_unc, sv_1z_unc)},
    '3q-bitflip': {'n': n3, 'states': make_bases(sv_0z_bf, sv_1z_bf)},
    '3q-phaseflip': {'n': n3, 'states': make_bases(sv_0z_pf, sv_1z_pf)},
    '[[5,1,3]]': {'n': n5, 'states': make_bases(sv_0z_513, sv_1z_513)},
}

print(f"Code states built in {time.time()-start:.1f}s")

# ============================================================
# Full (gamma, lambda) grid
# ============================================================

grid_size = 10
gammas = np.linspace(0.01, 0.40, grid_size)
lambdas = np.linspace(0.01, 0.40, grid_size)

results_grid = {}
winners_grid = np.empty((grid_size, grid_size), dtype='<U20')
advantage_grid = np.zeros((grid_size, grid_size))
oracle_gap_grid = np.zeros((grid_size, grid_size))

print(f"\n=== Full (γ, λ) grid: {grid_size}x{grid_size} ===\n")

for i, gamma in enumerate(gammas):
    for j, lam in enumerate(lambdas):
        avgs = {}
        for name, code_data in codes.items():
            avgs[name] = basis_avg_holevo(code_data['states'], gamma, lam, code_data['n'])

        # Best specialized code (oracle picks best 3-qubit code for this noise point)
        best_specialized = max(avgs['3q-bitflip'], avgs['3q-phaseflip'])
        best_specialized_name = '3q-bitflip' if avgs['3q-bitflip'] >= avgs['3q-phaseflip'] else '3q-phaseflip'

        # [[5,1,3]] advantage over best specialized
        advantage = avgs['[[5,1,3]]'] - best_specialized
        advantage_grid[i, j] = advantage

        # Oracle gap: [[5,1,3]] vs oracle-selected best overall
        all_avgs = {k: v for k, v in avgs.items()}
        winner = max(all_avgs, key=all_avgs.get)
        winners_grid[i, j] = winner

        # Gap between [[5,1,3]] and uncoded
        oracle_gap_grid[i, j] = avgs['[[5,1,3]]'] - avgs['uncoded']

        key = f"({gamma:.2f},{lam:.2f})"
        results_grid[key] = {
            name: round(float(v), 4) for name, v in avgs.items()
        }
        results_grid[key]['best_specialized'] = best_specialized_name
        results_grid[key]['advantage_513_vs_best_spec'] = round(float(advantage), 4)
        results_grid[key]['advantage_513_vs_uncoded'] = round(float(avgs['[[5,1,3]]'] - avgs['uncoded']), 4)
        results_grid[key]['winner'] = winner

# Print results
print("\n=== Winner Map ===\n")
print(f"{'γ\\λ':>8s}", end="")
for lam in lambdas:
    print(f"  {lam:.2f}", end="")
print()
for i, gamma in enumerate(gammas):
    print(f"γ={gamma:.2f} ", end="")
    for j in range(grid_size):
        w = winners_grid[i, j]
        sym = {'[[5,1,3]]': '5', 'uncoded': 'U', '3q-bitflip': 'B', '3q-phaseflip': 'P'}[w]
        print(f"     {sym}", end="")
    print()

print("\n=== [[5,1,3]] Advantage Over Best Specialized Code ===\n")
print(f"{'γ\\λ':>8s}", end="")
for lam in lambdas:
    print(f"  {lam:.2f}", end="")
print()
for i, gamma in enumerate(gammas):
    print(f"γ={gamma:.2f} ", end="")
    for j in range(grid_size):
        a = advantage_grid[i, j]
        print(f"  {a:+.3f}" if a != 0 else f"  {a:.3f}", end="")
    print()

print("\n=== [[5,1,3]] Advantage Over Uncoded ===\n")
print(f"{'γ\\λ':>8s}", end="")
for lam in lambdas:
    print(f"  {lam:.2f}", end="")
print()
for i, gamma in enumerate(gammas):
    print(f"γ={gamma:.2f} ", end="")
    for j in range(grid_size):
        a = oracle_gap_grid[i, j]
        print(f"  {a:+.3f}" if a != 0 else f"  {a:.3f}", end="")
    print()

# Bias ratio analysis
print("\n\n=== Bias Ratio Analysis ===\n")
print("Diagonal (balanced noise γ=λ):")
for gamma in gammas:
    key = f"({gamma:.2f},{gamma:.2f})"
    if key in results_grid:
        r = results_grid[key]
        print(f"  γ=λ={gamma:.2f}: [[5,1,3]]={r['[[5,1,3]]']:.3f}  best_spec={r[r['best_specialized']]:.3f}  uncoded={r['uncoded']:.3f}  adv={r['advantage_513_vs_best_spec']:.3f}")

# Statistics
all_advantages = advantage_grid.flatten()
print(f"\nAdvantage statistics (over best specialized code):")
print(f"  Min advantage:  {np.min(all_advantages):.4f}")
print(f"  Max advantage:  {np.max(all_advantages):.4f}")
print(f"  Mean advantage: {np.mean(all_advantages):.4f}")
print(f"  [[5,1,3]] wins {np.sum(all_advantages > 0)}/{len(all_advantages)} grid points ({100*np.sum(all_advantages>0)/len(all_advantages):.0f}%)")

# Does bias-aware oracle EVER beat [[5,1,3]]?
oracle_beats = 0
for i in range(grid_size):
    for j in range(grid_size):
        if winners_grid[i, j] != '[[5,1,3]]':
            oracle_beats += 1
print(f"  Non-[[5,1,3]] winners: {oracle_beats}/{grid_size**2}")

elapsed = time.time() - start
print(f"\nTotal time: {elapsed:.1f}s")

results = {
    'sprint': '022b',
    'description': 'Isotropy advantage map across (gamma, lambda) noise landscape',
    'grid_size': grid_size,
    'gammas': [round(float(g), 4) for g in gammas],
    'lambdas': [round(float(l), 4) for l in lambdas],
    'results': results_grid,
    'advantage_stats': {
        'min': round(float(np.min(all_advantages)), 4),
        'max': round(float(np.max(all_advantages)), 4),
        'mean': round(float(np.mean(all_advantages)), 4),
        'pct_513_wins': round(float(100 * np.sum(all_advantages > 0) / len(all_advantages)), 1),
    },
    'elapsed_seconds': round(elapsed, 1)
}

with open('results/sprint_022b_isotropy_advantage_map.json', 'w') as f:
    json.dump(results, f, indent=2)
print("Results saved to results/sprint_022b_isotropy_advantage_map.json")
