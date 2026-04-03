"""Sprint 106a: Spectral decomposition of fidelity susceptibility at g_c.

chi_F = (1/N) * sum_{n>0} |<n|H_field|0>|^2 / (E_n - E_0)^2

Decompose into individual excited-state contributions to understand
WHY walking gives alpha > 2 (super-first-order scaling).

For each q and size, compute:
- Full energy spectrum (lowest ~20 states)
- Matrix elements <n|H_field|0>
- Individual contributions chi_F^(n) = |<n|H_field|0>|^2 / (E_n - E_0)^2
- Cumulative fraction: how many states capture 90%, 99% of chi_F
"""
import numpy as np
import json, time, os
from scipy.sparse import lil_matrix, csr_matrix
from gpu_utils import eigsh
from db_utils import record

results = {
    'experiment': '106a_chi_f_spectral',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'data': {},
}

def save():
    outpath = os.path.join(os.path.dirname(__file__), 'results', 'sprint_106a_chi_f_spectral.json')
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    with open(outpath, 'w') as f:
        json.dump(results, f, indent=2, default=str)

def build_sq_potts_parts(n, q):
    """Build coupling and field parts separately: H(g) = H_coup + g * H_field.
    S_q Potts: field = sum of all cyclic permutations X^k for k=1..q-1."""
    dim = q**n
    H_coup = lil_matrix((dim, dim), dtype=float)
    H_field = lil_matrix((dim, dim), dtype=float)

    for idx in range(dim):
        c = []
        tmp = idx
        for _ in range(n):
            c.append(tmp % q)
            tmp //= q

        # Diagonal: Potts coupling -delta(s_i, s_j)
        diag = 0.0
        for site in range(n):
            nxt = (site + 1) % n
            if c[site] == c[nxt]:
                diag -= 1.0
        H_coup[idx, idx] = diag

        # Off-diagonal: S_q transverse field (all cyclic shifts)
        for site in range(n):
            for k in range(1, q):
                c2 = c.copy()
                c2[site] = (c[site] + k) % q
                idx2 = 0
                for i in range(n - 1, -1, -1):
                    idx2 = idx2 * q + c2[i]
                H_field[idx, idx2] += -1.0

    return csr_matrix(H_coup), csr_matrix(H_field)

# Configuration
configs = [
    # (q, n, k_states) — k_states = number of eigenstates to compute
    (2, 6, 20),
    (2, 8, 20),
    (2, 10, 20),
    (3, 6, 20),
    (3, 8, 20),
    (5, 6, 20),
    (5, 8, 20),
    (7, 6, 20),
]

print("=" * 70)
print("Sprint 106a: chi_F spectral decomposition at g_c = 1/q")
print("=" * 70)

for q, n, k_states in configs:
    g_c = 1.0 / q
    dim = q**n
    k_use = min(k_states, dim - 2)  # can't exceed dim

    print(f"\nq={q}, n={n} (dim={dim}), k={k_use} states...", flush=True)
    t0 = time.time()

    # Build Hamiltonian parts
    H_coup, H_field = build_sq_potts_parts(n, q)
    H = H_coup + g_c * H_field
    t_build = time.time() - t0
    print(f"  Built in {t_build:.1f}s", flush=True)

    # Get k lowest eigenstates
    evals, evecs = eigsh(H, k=k_use, which='SA')
    # Sort by energy (eigsh doesn't guarantee order)
    order = np.argsort(evals)
    evals = evals[order]
    evecs = evecs[:, order]

    E0 = evals[0]
    psi0 = evecs[:, 0]

    t_diag = time.time() - t0
    print(f"  Diagonalized in {t_diag:.1f}s", flush=True)

    # Compute matrix elements <n|H_field|0>
    H_field_psi0 = H_field.dot(psi0)

    contributions = []
    for i in range(1, k_use):
        psi_n = evecs[:, i]
        gap = evals[i] - E0
        if gap < 1e-12:
            continue
        me = np.dot(psi_n, H_field_psi0)  # <n|H_field|0>
        chi_n = me**2 / gap**2  # contribution to N*chi_F
        contributions.append({
            'level': i,
            'gap': float(gap),
            'matrix_element': float(me),
            'matrix_element_sq': float(me**2),
            'chi_contribution': float(chi_n),
        })

    # Sort by contribution (largest first)
    contributions.sort(key=lambda x: -x['chi_contribution'])

    # Total chi_F from these states
    chi_total = sum(c['chi_contribution'] for c in contributions)
    chi_F = chi_total / n  # per-site normalization

    # Cumulative analysis
    cumulative = 0.0
    n_90 = n_99 = len(contributions)
    for i, c in enumerate(contributions):
        cumulative += c['chi_contribution']
        frac = cumulative / chi_total if chi_total > 0 else 0
        c['cumulative_fraction'] = float(frac)
        if frac >= 0.90 and n_90 == len(contributions):
            n_90 = i + 1
        if frac >= 0.99 and n_99 == len(contributions):
            n_99 = i + 1

    # Top contributor analysis
    top1_frac = contributions[0]['chi_contribution'] / chi_total if contributions else 0
    top3_frac = sum(c['chi_contribution'] for c in contributions[:3]) / chi_total if contributions else 0

    dt = time.time() - t0

    # Store
    key = f"q{q}_n{n}"
    results['data'][key] = {
        'q': q, 'n': n, 'dim': dim, 'g_c': g_c,
        'chi_F': float(chi_F),
        'chi_total_Nchi': float(chi_total),
        'n_states_computed': k_use,
        'n_contributing': len(contributions),
        'top1_fraction': float(top1_frac),
        'top3_fraction': float(top3_frac),
        'n_for_90pct': n_90,
        'n_for_99pct': n_99,
        'E0': float(E0),
        'gap1': float(evals[1] - E0) if k_use > 1 else None,
        'contributions_top10': contributions[:10],
        'time_s': dt,
    }

    record(sprint=106, model='S_q_Potts', q=q, n=n,
           quantity='chi_F_spectral', value=chi_F, method='spectral_decomp')
    record(sprint=106, model='S_q_Potts', q=q, n=n,
           quantity='chi_F_top1_frac', value=top1_frac, method='spectral_decomp')
    record(sprint=106, model='S_q_Potts', q=q, n=n,
           quantity='chi_F_n_for_90pct', value=float(n_90), method='spectral_decomp')

    # Print summary
    print(f"  chi_F = {chi_F:.4f}")
    print(f"  Top contributor: level {contributions[0]['level']}, "
          f"gap={contributions[0]['gap']:.6f}, "
          f"|me|²={contributions[0]['matrix_element_sq']:.6f}, "
          f"frac={top1_frac:.4f}")
    print(f"  Top 3 capture {top3_frac*100:.1f}%, "
          f"90% needs {n_90} states, 99% needs {n_99} states")
    print(f"  Time: {dt:.1f}s")

    save()

    if dt > 50:
        print("  WARNING: approaching time limit, may skip larger configs")

# Summary table
print(f"\n{'=' * 70}")
print("SUMMARY: chi_F spectral decomposition")
print(f"{'q':>3} {'n':>4} {'chi_F':>10} {'gap1':>10} {'|me1|^2':>10} "
      f"{'top1%':>7} {'top3%':>7} {'n_90%':>6} {'n_99%':>6}")
print("-" * 70)
for key in sorted(results['data'].keys()):
    d = results['data'][key]
    top1_pct = d['top1_fraction'] * 100
    top3_pct = d['top3_fraction'] * 100
    c0 = d['contributions_top10'][0] if d['contributions_top10'] else {}
    me_sq = c0.get('matrix_element_sq', 0)
    print(f"{d['q']:>3} {d['n']:>4} {d['chi_F']:>10.4f} {d['gap1']:>10.6f} "
          f"{me_sq:>10.4f} {top1_pct:>6.1f}% {top3_pct:>6.1f}% "
          f"{d['n_for_90pct']:>6} {d['n_for_99pct']:>6}")

save()
print(f"\nResults saved.")
