#!/usr/bin/env python3
"""Sprint 077b: Conformal tower comparison S_q Potts vs hybrid at q=7.

Sprint 076b found S_q q=5 has 4-fold first excited state (S_5 symmetry)
vs hybrid 2-fold (Z_5 conjugate pairs, x_3=2.41 split).

q=7: S_7 has (q-1)=6 non-trivial permutations, so S_q field produces
6-fold degenerate first excited state. Hybrid has Z_7 conjugate pairs
(3 pairs for charge 1-3, plus self-conjugate charge 0 if q even).

Compare conformal tower ratios R_n = (E_n - E_0) / (E_1 - E_0).
"""
import numpy as np
from scipy.sparse import csr_matrix, kron as sp_kron, eye as sp_eye, diags
from gpu_utils import eigsh
import json, time

q = 7


def build_hamiltonian(n, q, g, model='sq_potts'):
    """Build Hamiltonian with periodic BC."""
    dim = q ** n
    potts_2site = np.zeros(q ** 2)
    for a in range(q):
        potts_2site[a * q + a] = 1.0
    potts_op = diags(potts_2site, 0, shape=(q ** 2, q ** 2), format='csr')

    if model == 'sq_potts':
        field_1site = np.ones((q, q)) - np.eye(q)
    else:
        X = np.zeros((q, q))
        for s in range(q):
            X[(s + 1) % q, s] = 1.0
        field_1site = X + X.T
    field_op = csr_matrix(field_1site)

    H = csr_matrix((dim, dim))
    for i in range(n - 1):
        left = q ** i
        right = q ** (n - i - 2)
        op = potts_op
        if left > 1:
            op = sp_kron(sp_eye(left), op, format='csr')
        if right > 1:
            op = sp_kron(op, sp_eye(right), format='csr')
        H = H - op

    diag_boundary = np.zeros(dim)
    for idx in range(dim):
        s0 = idx % q
        sn = (idx // (q ** (n - 1))) % q
        diag_boundary[idx] = 1.0 if s0 == sn else 0.0
    H = H - diags(diag_boundary, 0, shape=(dim, dim), format='csr')

    for i in range(n):
        left = q ** i
        right = q ** (n - i - 1)
        op = field_op.copy()
        if left > 1:
            op = sp_kron(sp_eye(left), op, format='csr')
        if right > 1:
            op = sp_kron(op, sp_eye(right), format='csr')
        H = H - g * op

    return H


# ---- MAIN ----
print("=" * 60, flush=True)
print("Sprint 077b: Conformal tower at q=7 — S_q vs hybrid", flush=True)
print("=" * 60, flush=True)

results = {
    'experiment': '077b_tower_q7',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'q': q,
}

# Critical points
gc_sq = 0.144   # From 077a
gc_hy = 0.535   # Known from prior sprints

k_evals = 16  # Get 16 lowest eigenvalues

for model, gc, label in [('sq_potts', gc_sq, 'S_q Potts'), ('hybrid', gc_hy, 'Hybrid')]:
    print(f"\n--- {label} at g_c={gc} ---", flush=True)
    model_results = {}

    for n in [4, 6]:
        dim = q ** n
        print(f"  n={n} (dim={dim:,})...", end=" ", flush=True)
        t0 = time.time()
        H = build_hamiltonian(n, q, gc, model)
        k = min(k_evals, dim - 1)
        if dim < 500:
            evals_all = np.linalg.eigvalsh(H.toarray())
            evals = evals_all[:k]
        else:
            evals, _ = eigsh(H, k=k, which='SA')
            evals = np.sort(evals)
        dt = time.time() - t0

        E0 = evals[0]
        gap1 = evals[1] - E0
        ratios = [(evals[i] - E0) / gap1 for i in range(1, len(evals))]

        print(f"{dt:.1f}s", flush=True)
        print(f"    E0={E0:.6f}, gap={gap1:.6f}, gap*N={gap1*n:.4f}", flush=True)
        print(f"    Ratios: {[f'{r:.4f}' for r in ratios[:12]]}", flush=True)

        # Check degeneracies: ratios within 1% of each other
        print(f"    Degeneracy pattern:", flush=True)
        i = 0
        deg_pattern = []
        while i < len(ratios):
            degen = 1
            while i + degen < len(ratios) and abs(ratios[i + degen] - ratios[i]) / max(abs(ratios[i]), 1e-10) < 0.02:
                degen += 1
            print(f"      R={ratios[i]:.4f} (×{degen})", flush=True)
            deg_pattern.append({'R': float(ratios[i]), 'degeneracy': degen})
            i += degen

        model_results[n] = {
            'evals': [float(e) for e in evals],
            'gap': float(gap1),
            'gap_N': float(gap1 * n),
            'ratios': [float(r) for r in ratios],
            'deg_pattern': deg_pattern,
            'time_s': float(dt),
        }

    results[model] = model_results

# Comparison
print(f"\n{'='*60}", flush=True)
print("COMPARISON at n=6", flush=True)
print(f"{'='*60}", flush=True)
for model, label in [('sq_potts', 'S_q Potts'), ('hybrid', 'Hybrid')]:
    r = results[model][6]
    print(f"\n{label}:", flush=True)
    print(f"  gap×N = {r['gap_N']:.4f}", flush=True)
    print(f"  First 10 ratios: {[f'{x:.3f}' for x in r['ratios'][:10]]}", flush=True)
    degs = r['deg_pattern']
    deg_strs = [f"{d['R']:.3f}x{d['degeneracy']}" for d in degs[:6]]
    print(f"  Degeneracy: {' | '.join(deg_strs)}", flush=True)

# Key comparison: first excited state degeneracy
sq_deg1 = results['sq_potts'][6]['deg_pattern'][0]['degeneracy']
hy_deg1 = results['hybrid'][6]['deg_pattern'][0]['degeneracy']
print(f"\n1st excited degeneracy: S_q={sq_deg1}, hybrid={hy_deg1}", flush=True)
print(f"Expected: S_q=6 (S_7 permutation), hybrid=2 (Z_7 conjugate pair)", flush=True)

# Save
with open("results/sprint_077b_tower_q7.json", "w") as f:
    json.dump(results, f, indent=2, default=str)
print("\nSaved to results/sprint_077b_tower_q7.json", flush=True)

from db_utils import record
for model, n_use in [('sq_potts', 6), ('hybrid', 6)]:
    r = results[model][n_use]
    if len(r['ratios']) > 2:
        record(sprint=77, model=model, q=q, n=n_use,
               quantity='x3_ratio', value=r['ratios'][2],
               method='conformal_tower', notes=f'R_3 = (E3-E0)/(E1-E0)')
    record(sprint=77, model=model, q=q, n=n_use,
           quantity='deg1', value=r['deg_pattern'][0]['degeneracy'],
           method='conformal_tower', notes='1st excited state degeneracy')
print("Recorded to DB.", flush=True)
