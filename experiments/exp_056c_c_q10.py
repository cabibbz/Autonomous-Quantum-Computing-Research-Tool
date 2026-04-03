#!/usr/bin/env python3
"""Sprint 056c: Measure central charge c(q=10) via exact diag.
Exact diag: n=4 (10^4=10000), n=5 (10^5=100000).
n=6 (10^6=1M) — attempt if timing allows.
g_c(q=10) = 0.684 from Sprint 052.

Predictions: log=1.59, power=1.76, quad=0.10
"""
import numpy as np, json, time
from scipy.sparse import kron, eye, csr_matrix, lil_matrix
from scipy.sparse.linalg import eigsh

def potts_hamiltonian_fast(q, n, g, J=1.0):
    """Build Potts H using lil_matrix for efficiency."""
    d = q
    dim = d**n
    H = lil_matrix((dim, dim))

    # Precompute single-site operators
    # For bond term: -J * sum_a P_a(i) P_a(i+1) = -J * delta(s_i, s_j)
    # Direct construction: iterate over basis states
    for idx in range(dim):
        # Decode state
        state = []
        tmp = idx
        for _ in range(n):
            state.append(tmp % d)
            tmp //= d
        # state[0] is site 0, state[1] is site 1, ...

        # Bond terms: -J * delta(s_i, s_j)
        for i in range(n - 1):
            if state[i] == state[i + 1]:
                H[idx, idx] -= J

        # Field terms: -g * (X + X†) on each site
        for i in range(n):
            # X|s_i> = |s_i+1 mod q>
            new_state_p = list(state)
            new_state_p[i] = (state[i] + 1) % q
            new_idx_p = sum(s * d**j for j, s in enumerate(new_state_p))
            H[idx, new_idx_p] -= g

            # X†|s_i> = |s_i-1 mod q>
            new_state_m = list(state)
            new_state_m[i] = (state[i] - 1) % q
            new_idx_m = sum(s * d**j for j, s in enumerate(new_state_m))
            H[idx, new_idx_m] -= g

    return H.tocsr()

def entanglement_entropy(psi, n, d, l):
    """Entropy of first l sites."""
    dim_A = d**l
    dim_B = d**(n - l)
    rho_A = psi.reshape(dim_A, dim_B)
    s = np.linalg.svd(rho_A, compute_uv=False)
    s2 = s[s > 1e-15]**2
    return float(-np.sum(s2 * np.log(s2)))

def extract_c_profile(S_list, n):
    """Central charge from entropy profile."""
    ls = np.arange(1, n)
    l_chord = (2 * n / np.pi) * np.sin(np.pi * ls / n)
    ln_chord = np.log(l_chord)
    S = np.array(S_list)
    quarter = max(1, n // 4)
    three_quarter = min(n - 1, 3 * n // 4)
    mask = (ls >= quarter) & (ls <= three_quarter)
    if mask.sum() < 2:
        mask = np.ones(len(ls), dtype=bool)
    A = np.vstack([ln_chord[mask], np.ones(mask.sum())]).T
    slope, _ = np.linalg.lstsq(A, S[mask], rcond=None)[0]
    return float(6 * slope)

q = 10
g_c = 0.684
data = []

for n in [4, 5]:
    t0 = time.time()
    dim = q**n
    print(f"q={q}, n={n}: dim={dim}, building H...", flush=True)
    H = potts_hamiltonian_fast(q, n, g_c)
    t_build = time.time() - t0
    print(f"  H built ({t_build:.1f}s), finding ground state...", flush=True)

    if dim > 500000:
        # For large matrices, use more Lanczos vectors
        E0, psi0 = eigsh(H, k=1, which='SA', maxiter=1000)
    else:
        E0, psi0 = eigsh(H, k=1, which='SA')
    E0 = float(E0[0])
    psi = psi0[:, 0]
    t_eig = time.time() - t0
    print(f"  E0={E0:.6f} ({t_eig:.1f}s), computing entropies...", flush=True)

    S_profile = []
    for l in range(1, n):
        S_profile.append(entanglement_entropy(psi, n, q, l))

    S_half = S_profile[n // 2 - 1]
    c_prof = extract_c_profile(S_profile, n)
    dt = time.time() - t0
    print(f"  n={n}: S_half={S_half:.6f}, c_profile={c_prof:.4f}, time={dt:.1f}s", flush=True)
    data.append({
        'method': 'exact_diag', 'n': n, 'dim': dim,
        'E0': E0, 'S_half': S_half, 'S_profile': S_profile,
        'c_profile': c_prof, 'time': dt
    })
    results = {'experiment': '056c', 'model': 'Potts', 'q': q, 'g_c': g_c,
               'predictions': {'log': 1.59, 'power': 1.76, 'quadratic': 0.10},
               'data': data}
    with open('results/exp_056c.json', 'w') as f:
        json.dump(results, f, indent=2)

    if dt > 200:
        print("Skipping larger sizes")
        break

# Try n=6 (10^6 = 1M) only if n=5 was fast enough
if data[-1]['time'] < 60:
    n = 6
    dim = q**n
    print(f"\nAttempting n={n}: dim={dim} (may be slow)...", flush=True)
    t0 = time.time()
    H = potts_hamiltonian_fast(q, n, g_c)
    t_build = time.time() - t0
    print(f"  H built ({t_build:.1f}s)", flush=True)
    if t_build < 120:
        print(f"  Finding ground state...", flush=True)
        E0, psi0 = eigsh(H, k=1, which='SA', maxiter=2000)
        E0 = float(E0[0])
        psi = psi0[:, 0]
        S_profile = []
        for l in range(1, n):
            S_profile.append(entanglement_entropy(psi, n, q, l))
        S_half = S_profile[n // 2 - 1]
        c_prof = extract_c_profile(S_profile, n)
        dt = time.time() - t0
        print(f"  n={n}: S_half={S_half:.6f}, c_profile={c_prof:.4f}, time={dt:.1f}s", flush=True)
        data.append({
            'method': 'exact_diag', 'n': n, 'dim': dim,
            'E0': E0, 'S_half': S_half, 'S_profile': S_profile,
            'c_profile': c_prof, 'time': dt
        })
    else:
        print(f"  H build too slow ({t_build:.0f}s), skipping")

# FSS from exact diag
if len(data) >= 2:
    n_list = [d['n'] for d in data]
    S_list = [d['S_half'] for d in data]
    c_fss = []
    for i in range(len(n_list) - 1):
        c = 6 * (S_list[i+1] - S_list[i]) / (np.log(n_list[i+1]) - np.log(n_list[i]))
        c_fss.append(float(c))
    results['c_fss_pairs'] = c_fss
    print(f"\nFSS pairwise c: {c_fss}")

results['data'] = data
with open('results/exp_056c.json', 'w') as f:
    json.dump(results, f, indent=2)

print(f"\n=== SUMMARY q={q} ===")
for d in data:
    print(f"  n={d['n']:3d}: S_half={d['S_half']:.6f}, c_profile={d.get('c_profile', 'N/A')}")
print(f"\nPredictions: log=1.59, power=1.76, quad=0.10")
print("Saved results/exp_056c.json")
