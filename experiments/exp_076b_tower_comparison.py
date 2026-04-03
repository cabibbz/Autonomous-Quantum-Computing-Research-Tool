#!/usr/bin/env python3
"""Sprint 076b: Conformal tower comparison — S_q Potts vs Hybrid at q=5.

At a genuine CFT critical point, the low-energy spectrum should organize
into conformal towers with specific gap ratios:
  x_k = (E_k - E_0) / (E_1 - E_0)  →  rational or irrational CFT values

For first-order transitions, the spectrum should show level repulsion
and NO clean tower structure (pseudo-critical point ≠ true CFT).

Compare at respective g_c values:
  S_q Potts: g_c ≈ 0.200
  Hybrid: g_c ≈ 0.438
"""
import numpy as np
from scipy.sparse import csr_matrix, kron as sp_kron, eye as sp_eye, diags
from gpu_utils import eigsh
import json, time

q = 5
gc_sq = 0.200
gc_hy = 0.438


def sq_potts_hamiltonian_periodic(n, q, g):
    """S_q Potts: Potts coupling + full S_q symmetric field."""
    dim = q ** n
    potts_2site = np.zeros(q ** 2)
    for a in range(q):
        potts_2site[a * q + a] = 1.0
    potts_op = diags(potts_2site, 0, shape=(q ** 2, q ** 2), format='csr')
    field_1site = np.ones((q, q)) - np.eye(q)
    field_op = csr_matrix(field_1site)
    H = csr_matrix((dim, dim))
    for i in range(n - 1):
        left = q ** i; right = q ** (n - i - 2)
        op = potts_op
        if left > 1: op = sp_kron(sp_eye(left), op, format='csr')
        if right > 1: op = sp_kron(op, sp_eye(right), format='csr')
        H = H - op
    diag_boundary = np.zeros(dim)
    for idx in range(dim):
        s0 = idx % q; sn = (idx // (q ** (n - 1))) % q
        diag_boundary[idx] = 1.0 if s0 == sn else 0.0
    H = H - diags(diag_boundary, 0, shape=(dim, dim), format='csr')
    for i in range(n):
        left = q ** i; right = q ** (n - i - 1)
        op = field_op.copy()
        if left > 1: op = sp_kron(sp_eye(left), op, format='csr')
        if right > 1: op = sp_kron(op, sp_eye(right), format='csr')
        H = H - g * op
    return H


def hybrid_hamiltonian_periodic(n, q, g):
    """Hybrid: Potts coupling + Z_q clock field (X + X†)."""
    dim = q ** n
    potts_2site = np.zeros(q ** 2)
    for a in range(q):
        potts_2site[a * q + a] = 1.0
    potts_op = diags(potts_2site, 0, shape=(q ** 2, q ** 2), format='csr')
    X = np.zeros((q, q))
    for s in range(q):
        X[(s + 1) % q, s] = 1.0
    XpXd = csr_matrix(X + X.T)
    H = csr_matrix((dim, dim))
    for i in range(n - 1):
        left = q ** i; right = q ** (n - i - 2)
        op = potts_op
        if left > 1: op = sp_kron(sp_eye(left), op, format='csr')
        if right > 1: op = sp_kron(op, sp_eye(right), format='csr')
        H = H - op
    diag_boundary = np.zeros(dim)
    for idx in range(dim):
        s0 = idx % q; sn = (idx // (q ** (n - 1))) % q
        diag_boundary[idx] = 1.0 if s0 == sn else 0.0
    H = H - diags(diag_boundary, 0, shape=(dim, dim), format='csr')
    for i in range(n):
        left = q ** i; right = q ** (n - i - 1)
        op = XpXd
        if left > 1: op = sp_kron(sp_eye(left), op, format='csr')
        if right > 1: op = sp_kron(op, sp_eye(right), format='csr')
        H = H - g * op
    return H


def get_spectrum(n, q, g, model, k=12):
    """Get k lowest eigenvalues."""
    if model == 'sq_potts':
        H = sq_potts_hamiltonian_periodic(n, q, g)
    else:
        H = hybrid_hamiltonian_periodic(n, q, g)
    dim = H.shape[0]
    if dim < 500:
        evals = np.sort(np.linalg.eigvalsh(H.toarray()))[:k]
    else:
        k_use = min(k, dim - 1)
        evals, _ = eigsh(H, k=k_use, which='SA')
        evals = np.sort(evals)
    return evals


# ---- MAIN ----
print("=" * 60, flush=True)
print("Sprint 076b: Conformal tower comparison q=5", flush=True)
print(f"S_q Potts g_c = {gc_sq}, Hybrid g_c = {gc_hy}", flush=True)
print("=" * 60, flush=True)

results = {
    'experiment': '076b_tower_comparison',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'q': q, 'gc_sq': gc_sq, 'gc_hy': gc_hy,
}

k_eig = 12
sizes = [4, 6, 8]

for model_name, gc in [('sq_potts', gc_sq), ('hybrid', gc_hy)]:
    print(f"\n{'='*60}", flush=True)
    print(f"Model: {model_name}, g_c = {gc}", flush=True)
    print(f"{'='*60}", flush=True)

    model_data = {}
    for n in sizes:
        dim = q ** n
        print(f"\nn={n} (dim={dim:,}):", flush=True)
        t0 = time.time()
        evals = get_spectrum(n, q, gc, model_name, k=k_eig)
        dt = time.time() - t0

        E0 = evals[0]
        gap1 = evals[1] - E0

        # Scaled gaps: x_k = (E_k - E_0) / gap1
        x = [(evals[i] - E0) / gap1 for i in range(len(evals))]

        # Degeneracies (levels within 1% of each other)
        degen_groups = []
        used = set()
        for i in range(1, len(evals)):
            if i in used:
                continue
            group = [i]
            for j in range(i + 1, len(evals)):
                if j not in used and abs(evals[j] - evals[i]) < 0.01 * gap1:
                    group.append(j)
            for j in group:
                used.add(j)
            degen_groups.append(group)

        print(f"  E0 = {E0:.8f}, gap1 = {gap1:.8f} ({dt:.1f}s)", flush=True)
        print(f"  Scaled spectrum x_k:", flush=True)
        for i in range(min(len(x), 10)):
            print(f"    x_{i} = {x[i]:.4f}", flush=True)

        print(f"  Degeneracy groups:", flush=True)
        for g_idx, group in enumerate(degen_groups[:6]):
            level_vals = [f"x_{k}={x[k]:.4f}" for k in group]
            print(f"    [{len(group)}] {', '.join(level_vals)}", flush=True)

        # Harmonic ratios R_k = x_k / k² (CFT prediction for free boson: R_k=1)
        print(f"  Harmonic ratios R_k = x_k / k²:", flush=True)
        for k in range(2, min(6, len(x))):
            r_k = x[k] / (k ** 2) if k > 0 else 0
            print(f"    R_{k} = {r_k:.4f}", flush=True)

        model_data[str(n)] = {
            'dim': dim, 'E0': float(E0), 'gap1': float(gap1),
            'evals': [float(e) for e in evals],
            'scaled_gaps': [float(xi) for xi in x],
            'degeneracies': [[int(j) for j in g] for g in degen_groups],
            'time_s': dt,
        }

    results[model_name] = model_data

# Cross-model comparison
print(f"\n{'='*60}", flush=True)
print("COMPARISON: Tower structure at n=8", flush=True)
print(f"{'='*60}", flush=True)

for n_str in ['8']:
    if n_str in results['sq_potts'] and n_str in results['hybrid']:
        sq = results['sq_potts'][n_str]
        hy = results['hybrid'][n_str]

        print(f"\n{'Level':>6} {'S_q x_k':>10} {'Hybrid x_k':>12} {'Diff%':>8}", flush=True)
        print("-" * 40, flush=True)
        n_levels = min(len(sq['scaled_gaps']), len(hy['scaled_gaps']), 10)
        for i in range(n_levels):
            xs = sq['scaled_gaps'][i]
            xh = hy['scaled_gaps'][i]
            if xs > 0:
                diff_pct = (xh - xs) / xs * 100
                print(f"  x_{i:>2} {xs:>10.4f} {xh:>12.4f} {diff_pct:>+7.1f}%", flush=True)
            else:
                print(f"  x_{i:>2} {xs:>10.4f} {xh:>12.4f}", flush=True)

        # Degeneracy comparison
        print(f"\n  S_q degeneracies: {[len(g) for g in sq['degeneracies'][:6]]}", flush=True)
        print(f"  Hybrid degeneracies: {[len(g) for g in hy['degeneracies'][:6]]}", flush=True)

        # Gap1 × N comparison
        sq_gN = sq['gap1'] * int(n_str)
        hy_gN = hy['gap1'] * int(n_str)
        print(f"\n  gap×N: S_q = {sq_gN:.6f}, Hybrid = {hy_gN:.6f}", flush=True)

# n-dependence of tower ratios (convergence test)
print(f"\n{'='*60}", flush=True)
print("CONVERGENCE: x_2 across sizes", flush=True)
print(f"{'='*60}", flush=True)
for model_name in ['sq_potts', 'hybrid']:
    x2_vals = []
    for n_str in ['4', '6', '8']:
        if n_str in results[model_name]:
            x2 = results[model_name][n_str]['scaled_gaps'][2] if len(results[model_name][n_str]['scaled_gaps']) > 2 else None
            if x2 is not None:
                x2_vals.append((int(n_str), x2))
    print(f"  {model_name}: {[(n, f'{x2:.4f}') for n, x2 in x2_vals]}", flush=True)
    if len(x2_vals) >= 2:
        # Check if converging or drifting
        drift = abs(x2_vals[-1][1] - x2_vals[-2][1]) / x2_vals[-2][1] * 100
        print(f"    n={x2_vals[-2][0]}→{x2_vals[-1][0]} drift: {drift:.2f}%", flush=True)

# Save
with open("results/sprint_076b_tower_comparison.json", "w") as f:
    json.dump(results, f, indent=2, default=str)
print("\nSaved to results/sprint_076b_tower_comparison.json", flush=True)

from db_utils import record
for model_name, gc in [('sq_potts', gc_sq), ('hybrid', gc_hy)]:
    if '8' in results[model_name]:
        data = results[model_name]['8']
        if len(data['scaled_gaps']) > 2:
            record(sprint=76, model=model_name, q=q, n=8,
                   quantity='x2_at_gc', value=data['scaled_gaps'][2],
                   method='exact_diag_periodic',
                   notes=f'conformal tower x_2 at g_c={gc}')
        record(sprint=76, model=model_name, q=q, n=8,
               quantity='gap1', value=data['gap1'],
               method='exact_diag_periodic',
               notes=f'gap at g_c={gc}')
print("Recorded to DB.", flush=True)
