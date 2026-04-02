#!/usr/bin/env python3
"""Sprint 076c: Δ·N scaling comparison — S_q Potts vs Hybrid at q=5.

First-order signature: Δ(N) ~ exp(-αN), so Δ·N GROWS with N.
Continuous signature: Δ(N) ~ 1/N, so Δ·N CONVERGES.

Also test at q=3 (both models equivalent, known continuous) as control.
And extend to q=7 for S_q Potts (expected more strongly first-order).
"""
import numpy as np
from scipy.sparse import csr_matrix, kron as sp_kron, eye as sp_eye, diags
from gpu_utils import eigsh
import json, time


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


def get_gap_at_gc(n, q, gc, model):
    """Get gap at critical point."""
    if model == 'sq_potts':
        H = sq_potts_hamiltonian_periodic(n, q, gc)
    else:
        H = hybrid_hamiltonian_periodic(n, q, gc)
    dim = H.shape[0]
    if dim < 500:
        evals = np.sort(np.linalg.eigvalsh(H.toarray()))
        return evals[1] - evals[0], evals[:min(6, len(evals))]
    k = min(6, dim - 1)
    evals, _ = eigsh(H, k=k, which='SA')
    evals = np.sort(evals)
    return evals[1] - evals[0], evals


# ---- MAIN ----
print("=" * 60, flush=True)
print("Sprint 076c: Δ·N scaling — S_q Potts vs Hybrid", flush=True)
print("=" * 60, flush=True)

results = {
    'experiment': '076c_gap_scaling_comparison',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
}

# Test configurations
# For q=3: both models are equivalent (control)
# For q=5: different (main test)
# For q=7: S_q Potts should be even more first-order
configs = [
    # (q, model, gc, sizes)
    (3, 'sq_potts', 0.333, [4, 6, 8, 10, 12]),
    (3, 'hybrid', 0.333, [4, 6, 8, 10, 12]),
    (5, 'sq_potts', 0.200, [4, 6, 8]),
    (5, 'hybrid', 0.438, [4, 6, 8]),
    (7, 'sq_potts', 0.120, [4, 6]),  # estimate g_c from field strength ratio
]

for q, model, gc, sizes in configs:
    key = f"q{q}_{model}"
    print(f"\n--- q={q}, {model}, g_c={gc} ---", flush=True)

    data = {'q': q, 'model': model, 'gc': gc, 'sizes': {}}

    for n in sizes:
        dim = q ** n
        if dim > 500000:
            # Skip if too large (>500k without GPU timing test)
            # But GPU handles up to ~10M
            if dim > 10000000:
                print(f"  n={n}: dim={dim:,} — skipping (too large)", flush=True)
                continue

        t0 = time.time()
        try:
            gap, evals = get_gap_at_gc(n, q, gc, model)
            dt = time.time() - t0
            if dt > 250:
                print(f"  n={n}: {dt:.0f}s — too slow, stopping", flush=True)
                data['sizes'][str(n)] = {'gap': float(gap), 'gap_N': float(gap * n),
                                          'dim': dim, 'time_s': dt}
                break

            gap_N = gap * n
            # Also compute gap2 for ratio
            gap2 = evals[2] - evals[0] if len(evals) > 2 else None
            ratio21 = (evals[2] - evals[0]) / gap if gap > 1e-15 and len(evals) > 2 else None

            print(f"  n={n} (dim={dim:,}): Δ₁={gap:.8f}, Δ₁·N={gap_N:.6f}", end="", flush=True)
            if ratio21 is not None:
                print(f", Δ₂/Δ₁={ratio21:.4f}", end="", flush=True)
            print(f" ({dt:.1f}s)", flush=True)

            data['sizes'][str(n)] = {
                'gap': float(gap), 'gap_N': float(gap_N),
                'gap2': float(gap2) if gap2 is not None else None,
                'ratio21': float(ratio21) if ratio21 is not None else None,
                'dim': dim, 'time_s': dt,
                'evals': [float(e) for e in evals[:6]],
            }
        except Exception as e:
            dt = time.time() - t0
            print(f"  n={n}: FAILED ({e}) ({dt:.1f}s)", flush=True)

    results[key] = data

# ============================================================
# ANALYSIS
# ============================================================
print(f"\n{'='*60}", flush=True)
print("ANALYSIS: Δ·N trends", flush=True)
print(f"{'='*60}", flush=True)

for key, data in results.items():
    if key in ['experiment', 'timestamp']:
        continue
    q = data['q']
    model = data['model']
    gc = data['gc']

    ns = []
    gap_Ns = []
    for n_str in sorted(data['sizes'].keys(), key=int):
        s = data['sizes'][n_str]
        ns.append(int(n_str))
        gap_Ns.append(s['gap_N'])

    if len(ns) < 2:
        continue

    print(f"\nq={q} {model} (g_c={gc}):", flush=True)
    print(f"  {'n':>4} {'Δ·N':>12} {'Δ(Δ·N)':>12} {'%change':>10}", flush=True)
    print(f"  {'-'*40}", flush=True)
    for i, (n, gN) in enumerate(zip(ns, gap_Ns)):
        if i == 0:
            print(f"  {n:>4} {gN:>12.6f}", flush=True)
        else:
            d = gN - gap_Ns[i-1]
            pct = d / gap_Ns[i-1] * 100
            print(f"  {n:>4} {gN:>12.6f} {d:>+12.6f} {pct:>+9.2f}%", flush=True)

    # Fit log(Δ) vs N: first-order → linear (exponential gap), continuous → -log(N)
    if len(ns) >= 3:
        gaps = [data['sizes'][str(n)]['gap'] for n in ns]
        log_gaps = np.log(gaps)
        ns_arr = np.array(ns, dtype=float)

        # Linear fit: log(Δ) = a - b·N (exponential: first-order)
        coeffs_exp = np.polyfit(ns_arr, log_gaps, 1)
        residuals_exp = log_gaps - np.polyval(coeffs_exp, ns_arr)
        rms_exp = np.sqrt(np.mean(residuals_exp ** 2))

        # Power-law fit: log(Δ) = a - z·log(N) (continuous: z=1 for CFT)
        log_ns = np.log(ns_arr)
        coeffs_pow = np.polyfit(log_ns, log_gaps, 1)
        residuals_pow = log_gaps - np.polyval(coeffs_pow, log_ns)
        rms_pow = np.sqrt(np.mean(residuals_pow ** 2))

        z_eff = -coeffs_pow[0]
        alpha_eff = -coeffs_exp[0]

        print(f"\n  Exponential fit: Δ ~ exp(-{alpha_eff:.4f}·N), RMS={rms_exp:.6f}", flush=True)
        print(f"  Power-law fit:   Δ ~ N^(-{z_eff:.4f}), RMS={rms_pow:.6f}", flush=True)

        if rms_exp < rms_pow * 0.5:
            print(f"  → EXPONENTIAL (first-order) fits BETTER by {rms_pow/rms_exp:.1f}x", flush=True)
        elif rms_pow < rms_exp * 0.5:
            print(f"  → POWER-LAW (continuous) fits BETTER by {rms_exp/rms_pow:.1f}x", flush=True)
        else:
            print(f"  → Inconclusive (comparable fits)", flush=True)

        data['fit_exp'] = {'alpha': float(alpha_eff), 'rms': float(rms_exp)}
        data['fit_pow'] = {'z': float(z_eff), 'rms': float(rms_pow)}

# Cross-model comparison at q=5
print(f"\n{'='*60}", flush=True)
print("DIRECT COMPARISON: q=5 S_q vs Hybrid", flush=True)
print(f"{'='*60}", flush=True)

if 'q5_sq_potts' in results and 'q5_hybrid' in results:
    sq = results['q5_sq_potts']
    hy = results['q5_hybrid']

    print(f"\n{'n':>4} {'S_q Δ·N':>12} {'Hybrid Δ·N':>12} {'S_q Δ₂/Δ₁':>12} {'Hyb Δ₂/Δ₁':>12}", flush=True)
    print("-" * 55, flush=True)
    for n_str in sorted(set(sq['sizes'].keys()) & set(hy['sizes'].keys()), key=int):
        s = sq['sizes'][n_str]
        h = hy['sizes'][n_str]
        sq_r = f"{s.get('ratio21', 0):.4f}" if s.get('ratio21') else "—"
        hy_r = f"{h.get('ratio21', 0):.4f}" if h.get('ratio21') else "—"
        print(f"{n_str:>4} {s['gap_N']:>12.6f} {h['gap_N']:>12.6f} {sq_r:>12} {hy_r:>12}", flush=True)

    # Verdict
    sq_gNs = [sq['sizes'][str(n)]['gap_N'] for n in [4, 6, 8] if str(n) in sq['sizes']]
    hy_gNs = [hy['sizes'][str(n)]['gap_N'] for n in [4, 6, 8] if str(n) in hy['sizes']]

    if len(sq_gNs) >= 3:
        sq_d1 = sq_gNs[1] - sq_gNs[0]
        sq_d2 = sq_gNs[2] - sq_gNs[1]
        sq_accel = sq_d2 - sq_d1
        print(f"\n  S_q Potts: diffs = {sq_d1:+.6f}, {sq_d2:+.6f}, accel = {sq_accel:+.6f}", flush=True)
        if sq_accel > 0:
            print(f"    → ACCELERATING (first-order consistent)", flush=True)
        else:
            print(f"    → DECELERATING (continuous consistent)", flush=True)

    if len(hy_gNs) >= 3:
        hy_d1 = hy_gNs[1] - hy_gNs[0]
        hy_d2 = hy_gNs[2] - hy_gNs[1]
        hy_accel = hy_d2 - hy_d1
        print(f"  Hybrid:   diffs = {hy_d1:+.6f}, {hy_d2:+.6f}, accel = {hy_accel:+.6f}", flush=True)
        if hy_accel > 0:
            print(f"    → ACCELERATING (first-order consistent)", flush=True)
        else:
            print(f"    → DECELERATING (continuous consistent)", flush=True)

# Save
with open("results/sprint_076c_gap_scaling_comparison.json", "w") as f:
    json.dump(results, f, indent=2, default=str)
print(f"\nSaved to results/sprint_076c_gap_scaling_comparison.json", flush=True)

from db_utils import record
for key, data in results.items():
    if key in ['experiment', 'timestamp']:
        continue
    q = data['q']
    model = data['model']
    for n_str, s in data['sizes'].items():
        n = int(n_str)
        record(sprint=76, model=model, q=q, n=n,
               quantity='gap_N_at_gc', value=s['gap_N'],
               method='exact_diag_periodic',
               notes=f'Δ·N at g_c={data["gc"]}')
    if 'fit_pow' in data:
        record(sprint=76, model=model, q=q, n=0,
               quantity='z_eff', value=data['fit_pow']['z'],
               method='power_law_fit',
               notes=f'Δ~N^(-z), RMS={data["fit_pow"]["rms"]:.6f}')
print("Recorded to DB.", flush=True)
