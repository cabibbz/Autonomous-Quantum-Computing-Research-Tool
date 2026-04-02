#!/usr/bin/env python3
"""Sprint 068c: First-order vs continuous diagnostic for 2D hybrid at q=5,7,10

Three independent tests:
1. dE₀/dg kink: First-order → discontinuous slope. Continuous → smooth.
   Fine scan near g_c for q=5 L=3 (10 pts × 25s = 250s).
2. q-dependence at L=2: If transition becomes first-order above some q*,
   gap*L behavior should qualitatively change.
3. Gap scaling: In 2D, first-order → gap ~ exp(-α·L^d). Continuous → gap ~ L^{-z}.
   Compare gap(g_c) between L=2 and L=3 at q=5.
"""
import numpy as np, json, time
from scipy.sparse import csr_matrix, kron as sp_kron, eye as sp_eye, diags
from scipy.sparse.linalg import eigsh

def build_H_2d(Lx, Ly, q, g):
    n = Lx*Ly; dim = q**n
    potts_2site = np.zeros(q**2)
    for a in range(q): potts_2site[a*q+a] = 1.0
    potts_op = diags(potts_2site, 0, shape=(q**2, q**2), format='csr')
    X = np.zeros((q,q))
    for s in range(q): X[(s+1)%q, s] = 1.0
    XpXd = csr_matrix(X + X.T)
    H = csr_matrix((dim, dim))
    bonds = []
    for y in range(Ly):
        for x in range(Lx):
            site = y*Lx+x
            nbr_h = y*Lx+(x+1)%Lx
            if nbr_h != site: bonds.append((min(site,nbr_h), max(site,nbr_h)))
            nbr_v = ((y+1)%Ly)*Lx+x
            if nbr_v != site: bonds.append((min(site,nbr_v), max(site,nbr_v)))
    bonds = list(set(bonds))
    for (i,j) in bonds:
        if j == i+1:
            left = q**i; right = q**(n-i-2); op = potts_op
            if left > 1: op = sp_kron(sp_eye(left), op, format='csr')
            if right > 1: op = sp_kron(op, sp_eye(right), format='csr')
        else:
            dv = np.zeros(dim)
            for idx in range(dim):
                si = (idx//(q**i))%q; sj = (idx//(q**j))%q
                dv[idx] = 1.0 if si == sj else 0.0
            op = diags(dv, 0, shape=(dim,dim), format='csr')
        H = H - op
    for i in range(n):
        left = q**i; right = q**(n-i-1); op = XpXd.copy()
        if left > 1: op = sp_kron(sp_eye(left), op, format='csr')
        if right > 1: op = sp_kron(op, sp_eye(right), format='csr')
        H = H - g * op
    return H

def get_evals(L, q, g, k=4):
    n = L*L; dim = q**n
    H = build_H_2d(L, L, q, g)
    if dim < 500:
        return np.linalg.eigvalsh(H.toarray())[:k]
    return np.sort(eigsh(H, k=k, which='SA')[0])

results = {'experiment': '068c_first_order_test', 'timestamp': time.strftime('%Y-%m-%d %H:%M:%S')}

# === TEST 1: dE₀/dg kink test at q=5, L=3 ===
print("=" * 60)
print("TEST 1: Energy slope (dE₀/dg) at q=5, L=3")
print("Fine scan near g_c ≈ 1.59")
print("=" * 60)

# Fine scan: 10 points near g_c
g_fine = np.linspace(1.1, 1.9, 10)
q5_L3_fine = []
for gi, g in enumerate(g_fine):
    t0 = time.time()
    evals = get_evals(3, 5, g)
    dt = time.time() - t0
    E0_per_site = evals[0] / 9
    gap = evals[1] - evals[0]
    q5_L3_fine.append({
        'g': round(g, 5),
        'E0_per_site': round(float(E0_per_site), 8),
        'gap': round(float(gap), 8),
        'gap_x_L': round(float(gap * 3), 6),
        'time_s': round(dt, 1)
    })
    print(f"  g={g:.3f}: E0/N={E0_per_site:.6f}, gap={gap:.6f}, gap*L={gap*3:.4f} ({dt:.1f}s)")

# Compute dE/dg numerically
print("\n  dE₀/dg (numerical):")
for i in range(1, len(q5_L3_fine)):
    dg = q5_L3_fine[i]['g'] - q5_L3_fine[i-1]['g']
    dE = q5_L3_fine[i]['E0_per_site'] - q5_L3_fine[i-1]['E0_per_site']
    slope = dE / dg
    g_mid = (q5_L3_fine[i]['g'] + q5_L3_fine[i-1]['g']) / 2
    print(f"    g={g_mid:.3f}: dE/dg = {slope:.4f}")

results['q5_L3_fine'] = q5_L3_fine

# === TEST 2: q-dependence at L=2 (dim = q^4) ===
print("\n" + "=" * 60)
print("TEST 2: q-dependence at L=2 — q=2,3,5,7,10")
print("=" * 60)

# Expected 2D g_c from 1D scaling (ratio ~3.6)
gc_1d = {2: 0.250, 3: 0.333, 5: 0.441, 7: 0.535, 10: 0.684}
gc_2d_est = {q: gc_1d[q] * 3.6 for q in gc_1d}

L2_results = {}
for q in [2, 3, 5, 7, 10]:
    dim = q**4
    gc_est = gc_2d_est[q]
    g_scan = np.linspace(gc_est * 0.3, gc_est * 2.5, 30)
    print(f"\n  q={q} (dim={dim}, g_c estimate ≈ {gc_est:.2f})")
    data = []
    t0 = time.time()
    for g in g_scan:
        evals = get_evals(2, q, g)
        gap = evals[1] - evals[0]
        data.append({
            'g': round(g, 5),
            'gap': round(float(gap), 10),
            'gap_x_L': round(float(gap * 2), 8),
            'E0_per_site': round(float(evals[0] / 4), 8)
        })
    dt = time.time() - t0
    L2_results[q] = data
    print(f"    Time: {dt:.1f}s")

    # Print gap*L at a few key g values
    for d in data:
        if abs(d['g'] - gc_est) < (g_scan[1] - g_scan[0]):
            print(f"    gap*L at g≈{d['g']:.2f}: {d['gap_x_L']:.4f}")

results['L2_all_q'] = {str(q): v for q, v in L2_results.items()}

# === TEST 3: Gap scaling at g_c ===
print("\n" + "=" * 60)
print("TEST 3: Gap at g_c — L=2 vs L=3 comparison")
print("=" * 60)

# Use the crossing values from 068b
gc_2d = {2: 0.771, 3: 1.267, 5: 1.588}

for q in [2, 3, 5]:
    gc = gc_2d[q]
    print(f"\n  q={q}, g_c={gc:.3f}:")

    # L=2
    evals_L2 = get_evals(2, q, gc)
    gap_L2 = evals_L2[1] - evals_L2[0]

    # L=3 (q=5 will be slow)
    t0 = time.time()
    evals_L3 = get_evals(3, q, gc)
    gap_L3 = evals_L3[1] - evals_L3[0]
    dt = time.time() - t0

    ratio = gap_L3 / gap_L2
    # Continuous (z=1): gap ~ 1/L, ratio = 2/3 = 0.667
    # First-order: gap ~ exp(-α·L), ratio << 0.667
    print(f"    L=2: gap = {gap_L2:.6f}, gap*L = {gap_L2*2:.4f}")
    print(f"    L=3: gap = {gap_L3:.6f}, gap*L = {gap_L3*3:.4f} ({dt:.1f}s)")
    print(f"    Ratio gap(L=3)/gap(L=2) = {ratio:.4f}")
    print(f"    Continuous (z=1) prediction: {2/3:.4f}")
    print(f"    If z≠1: gap(L=3)/gap(L=2) = (2/3)^z → z = {np.log(ratio)/np.log(2/3):.2f}")

    results[f'gap_scaling_q{q}'] = {
        'gc': gc,
        'gap_L2': float(gap_L2),
        'gap_L3': float(gap_L3),
        'ratio': float(ratio),
        'z_eff': float(np.log(ratio)/np.log(2/3)) if ratio > 0 else None
    }

# === q=5 L=3 at q=5: also do q=7 L=2 vs theoretical ===
print("\n" + "=" * 60)
print("TEST 4: Extended q scan at L=2 — crossing with hypothetical L=3")
print("=" * 60)

# For q=7,10 we can't do L=3. But at L=2, check gap*L behavior
for q in [7, 10]:
    dim_L2 = q**4
    gc_est = gc_1d[q] * 3.6
    g_vals = [gc_est * 0.5, gc_est * 0.75, gc_est, gc_est * 1.25, gc_est * 1.5]
    print(f"\n  q={q} (L=2, dim={dim_L2}):")
    for g in g_vals:
        evals = get_evals(2, q, g)
        gap = evals[1] - evals[0]
        gap2 = evals[2] - evals[0]
        print(f"    g={g:.3f}: gap*L={gap*2:.4f}, gap2/gap1={gap2/gap:.3f}")

# === SUMMARY ===
print("\n" + "=" * 60)
print("SUMMARY")
print("=" * 60)

print("\nGap scaling ratios gap(L=3)/gap(L=2) at g_c:")
for q in [2, 3, 5]:
    r = results[f'gap_scaling_q{q}']
    print(f"  q={q}: ratio = {r['ratio']:.4f}, z_eff = {r['z_eff']:.2f} (continuous→z=1)")

print("\nVerdict:")
print("  Continuous: ratio ≈ 0.667 (z=1)")
print("  First-order: ratio << 0.667 (exponential gap)")

with open("results/sprint_068c_first_order.json", "w") as f:
    json.dump(results, f, indent=2, default=str)
print("\nSaved results/sprint_068c_first_order.json")
