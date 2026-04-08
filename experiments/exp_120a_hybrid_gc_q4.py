"""Sprint 120a: Find hybrid Potts-clock g_c at q=4 via gap×N crossings.

Hybrid model: H = -Σ δ(s_i,s_j) - g(X + X†)
Known hybrid g_c: q=3→0.333, q=5→0.438. Interpolation: q=4 ≈ 0.38-0.39.

Sizes: n=6 (4096), n=8 (65536), n=10 (1M) — all feasible with GPU.
Scan g ∈ [0.33, 0.43] with 30 points, find gap×N crossing.
"""
import numpy as np
import json, time, os
from scipy.sparse import coo_matrix, csr_matrix
from gpu_utils import eigsh
from db_utils import record

q = 4
results = {
    'experiment': '120a_hybrid_gc_q4',
    'sprint': 120,
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'model': 'hybrid',
    'q': q,
    'data': {},
}

def save():
    outpath = os.path.join(os.path.dirname(__file__), '..', 'results', 'sprint_120a_hybrid_gc_q4.json')
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    with open(outpath, 'w') as f:
        json.dump(results, f, indent=2, default=str)


def build_hybrid_H(n, q, g):
    """Build hybrid Potts-clock Hamiltonian: δ coupling + (X+X†) field."""
    dim = q**n
    all_idx = np.arange(dim, dtype=np.int64)
    digits = np.zeros((dim, n), dtype=np.int64)
    tmp = all_idx.copy()
    for site in range(n):
        digits[:, site] = tmp % q
        tmp //= q
    powers = q ** np.arange(n, dtype=np.int64)

    # Coupling: -δ(s_i, s_{i+1}) periodic
    diag_vals = np.zeros(dim, dtype=np.float64)
    for site in range(n):
        nxt = (site + 1) % n
        diag_vals -= (digits[:, site] == digits[:, nxt]).astype(np.float64)

    # Field: -g(X + X†) at each site
    rows_list, cols_list, vals_list = [], [], []
    for site in range(n):
        pw = powers[site]
        old_digit = digits[:, site]
        for k in [1, q - 1]:  # X and X†
            new_digit = (old_digit + k) % q
            delta = (new_digit.astype(np.int64) - old_digit.astype(np.int64)) * pw
            new_idx = all_idx + delta
            rows_list.append(new_idx)
            cols_list.append(all_idx)
            vals_list.append(np.full(dim, -g, dtype=np.float64))
    rows = np.concatenate(rows_list)
    cols = np.concatenate(cols_list)
    vals = np.concatenate(vals_list)

    H = coo_matrix((vals, (rows, cols)), shape=(dim, dim)).tocsr()
    H.setdiag(H.diagonal() + diag_vals)
    return H


# Phase 1: Coarse scan
sizes = [6, 8, 10]
g_vals = np.linspace(0.33, 0.43, 30)

print("=" * 70)
print(f"HYBRID q={q} g_c SCAN: gap×N crossings")
print(f"Sizes: {sizes}, g range: [{g_vals[0]:.3f}, {g_vals[-1]:.3f}], {len(g_vals)} points")
print("=" * 70)

gap_data = {}  # {n: [(g, gap*N), ...]}

for n in sizes:
    dim = q**n
    gap_data[n] = []
    print(f"\n  n={n} (dim={dim:,}):", flush=True)
    t0 = time.time()

    for i, g in enumerate(g_vals):
        H = build_hybrid_H(n, q, g)
        evals, _ = eigsh(H, k=4, which='SA')
        evals.sort()
        gap = evals[1] - evals[0]
        gN = gap * n
        gap_data[n].append((float(g), float(gN)))
        if i % 10 == 0:
            print(f"    g={g:.4f}: gap×N={gN:.6f}", flush=True)

    dt = time.time() - t0
    print(f"  n={n} done in {dt:.1f}s", flush=True)
    results['data'][f'n{n}_scan'] = gap_data[n]
    save()

# Find crossings between consecutive size pairs
print(f"\n{'=' * 70}")
print("CROSSING ANALYSIS")
print("=" * 70)

crossings = []
for i in range(len(sizes) - 1):
    n1, n2 = sizes[i], sizes[i+1]
    g1_arr = np.array([p[0] for p in gap_data[n1]])
    gN1_arr = np.array([p[1] for p in gap_data[n1]])
    gN2_arr = np.array([p[1] for p in gap_data[n2]])
    diff = gN1_arr - gN2_arr
    # Find sign changes
    for j in range(len(diff) - 1):
        if diff[j] * diff[j+1] < 0:
            # Linear interpolation
            frac = diff[j] / (diff[j] - diff[j+1])
            g_cross = g1_arr[j] + frac * (g1_arr[j+1] - g1_arr[j])
            crossings.append((n1, n2, float(g_cross)))
            print(f"  ({n1},{n2}) crossing at g = {g_cross:.6f}")

if crossings:
    gc_est = np.mean([c[2] for c in crossings])
    gc_spread = np.max([c[2] for c in crossings]) - np.min([c[2] for c in crossings])
    print(f"\n  g_c estimate: {gc_est:.6f} (spread: {gc_spread:.6f})")
    results['gc_estimate'] = float(gc_est)
    results['gc_spread'] = float(gc_spread)
    results['crossings'] = crossings

    # Phase 2: Fine scan around crossing
    print(f"\n{'=' * 70}")
    print(f"FINE SCAN around g={gc_est:.4f}")
    print("=" * 70)

    g_fine = np.linspace(gc_est - 0.01, gc_est + 0.01, 20)
    fine_data = {}
    for n in sizes:
        dim = q**n
        fine_data[n] = []
        print(f"\n  n={n}:", flush=True)
        for g in g_fine:
            H = build_hybrid_H(n, q, g)
            evals, _ = eigsh(H, k=4, which='SA')
            evals.sort()
            gap = evals[1] - evals[0]
            gN = gap * n
            fine_data[n].append((float(g), float(gN)))
        results['data'][f'n{n}_fine'] = fine_data[n]

    # Fine crossings
    fine_crossings = []
    for i in range(len(sizes) - 1):
        n1, n2 = sizes[i], sizes[i+1]
        g1_arr = np.array([p[0] for p in fine_data[n1]])
        gN1_arr = np.array([p[1] for p in fine_data[n1]])
        gN2_arr = np.array([p[1] for p in fine_data[n2]])
        diff = gN1_arr - gN2_arr
        for j in range(len(diff) - 1):
            if diff[j] * diff[j+1] < 0:
                frac = diff[j] / (diff[j] - diff[j+1])
                g_cross = g1_arr[j] + frac * (g1_arr[j+1] - g1_arr[j])
                fine_crossings.append((n1, n2, float(g_cross)))
                print(f"  Fine ({n1},{n2}) crossing at g = {g_cross:.6f}")

    if fine_crossings:
        gc_fine = np.mean([c[2] for c in fine_crossings])
        print(f"\n  REFINED g_c = {gc_fine:.6f}")
        results['gc_fine'] = float(gc_fine)
        results['fine_crossings'] = fine_crossings
        record(sprint=120, model='hybrid', q=q, n=0,
               quantity='g_c', value=float(gc_fine), method='gap_crossing')

save()
print(f"\nResults saved. g_c(hybrid, q=4) determined.")
