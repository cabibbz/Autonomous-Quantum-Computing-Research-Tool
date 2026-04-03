"""Sprint 109b: Extend q=4 to n=11 (dim=4.2M, GPU) and q=2 to n=18 (dim=262k).

q=4 BKT: add one more size to test log correction form.
q=2 Ising: extend to n=18 for a third test case (real CFT, exact nu=1, alpha=1).
"""
import numpy as np
import json, time, os
from scipy.sparse import lil_matrix, csr_matrix
from gpu_utils import eigsh
from db_utils import record

results = {
    'experiment': '109b_q4_q2_extend',
    'sprint': 109,
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'data': {},
}

def save():
    outpath = os.path.join(os.path.dirname(__file__), 'results', 'sprint_109b_q4_q2_extend.json')
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    with open(outpath, 'w') as f:
        json.dump(results, f, indent=2, default=str)

def build_sq_potts_parts(n, q):
    """Build coupling and field parts separately for S_q Potts."""
    dim = q**n
    H_coup = lil_matrix((dim, dim), dtype=float)
    H_field = lil_matrix((dim, dim), dtype=float)
    for idx in range(dim):
        c = []
        tmp = idx
        for _ in range(n):
            c.append(tmp % q)
            tmp //= q
        diag = 0.0
        for site in range(n):
            nxt = (site + 1) % n
            if c[site] == c[nxt]:
                diag -= 1.0
        H_coup[idx, idx] = diag
        for site in range(n):
            for k in range(1, q):
                c2 = c.copy()
                c2[site] = (c[site] + k) % q
                idx2 = 0
                for i in range(n - 1, -1, -1):
                    idx2 = idx2 * q + c2[i]
                H_field[idx, idx2] += -1.0
    return csr_matrix(H_coup), csr_matrix(H_field)

def run_spectral(q, sizes, k_need, label):
    g_c = 1.0 / q
    print(f"\n{'=' * 70}")
    print(f"q={q} ({label}): chi_F spectral decomposition")
    print("=" * 70)

    for n in sizes:
        dim = q**n
        print(f"\n  n={n} (dim={dim:,})...", end="", flush=True)
        t0 = time.time()

        H_coup, H_field = build_sq_potts_parts(n, q)
        H = H_coup + g_c * H_field
        k_use = min(k_need, dim - 2)
        evals, evecs = eigsh(H, k=k_use, which='SA')
        order = np.argsort(evals)
        evals = evals[order]
        evecs = evecs[:, order]

        E0 = evals[0]
        psi0 = evecs[:, 0]
        H_field_psi0 = H_field.dot(psi0)

        best_chi = 0; best_gap = 0; best_me_sq = 0; total_chi = 0
        spectral_gap = None
        for i in range(1, k_use):
            gap = evals[i] - E0
            if gap < 1e-12:
                continue
            if spectral_gap is None:
                spectral_gap = gap
            me = np.dot(evecs[:, i], H_field_psi0)
            chi_n = me**2 / gap**2
            total_chi += chi_n
            if chi_n > best_chi:
                best_chi = chi_n
                best_gap = gap
                best_me_sq = float(me**2)

        chi_F = total_chi / n
        dominant_frac = best_chi / total_chi if total_chi > 0 else 0
        dt = time.time() - t0
        print(f" gap_m={best_gap:.6f}, |me|²={best_me_sq:.4f}, chi_F={chi_F:.4f}, "
              f"frac={dominant_frac:.4f} ({dt:.1f}s)")

        key = f"q{q}_n{n}"
        results['data'][key] = {
            'n': n, 'dim': dim, 'q': q,
            'gap_multiplet': float(best_gap),
            'gap_spectral': float(spectral_gap) if spectral_gap else None,
            'me_sq': best_me_sq,
            'chi_F': float(chi_F),
            'dominant_fraction': float(dominant_frac),
            'E0': float(E0),
            'time_s': round(dt, 1),
        }
        save()

        record(sprint=109, model='sq_potts', q=q, n=n,
               quantity='chi_F', value=float(chi_F), method='spectral_exact_gpu')
        record(sprint=109, model='sq_potts', q=q, n=n,
               quantity='gap_multiplet', value=float(best_gap), method='spectral_exact_gpu')
        record(sprint=109, model='sq_potts', q=q, n=n,
               quantity='me_sq', value=best_me_sq, method='spectral_exact_gpu')

# q=2: full range n=4-18 (2^18 = 262144, CPU/GPU trivial)
# Need more eigenstates for q=2 since dominant fraction < 100%
run_spectral(2, [4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18], k_need=10, label="Ising, real CFT")

# q=4: existing n=4-10 + new n=11 (4^11 = 4.2M, GPU)
run_spectral(4, [4, 5, 6, 7, 8, 9, 10, 11], k_need=8, label="BKT boundary")

save()

# Pairwise analysis for both
print(f"\n{'=' * 70}")
print("PAIRWISE EXPONENTS")
print("=" * 70)

for q_val in [2, 4]:
    Ns = []; gaps = []; mes = []; chis = []
    for key, d in sorted(results['data'].items()):
        if d['q'] == q_val:
            Ns.append(d['n']); gaps.append(d['gap_multiplet'])
            mes.append(d['me_sq']); chis.append(d['chi_F'])
    Ns = np.array(Ns, dtype=float)
    gaps = np.array(gaps); mes = np.array(mes); chis = np.array(chis)
    log_N = np.log(Ns)
    log_chi = np.log(chis)

    pw_a = []
    for i in range(len(Ns)-1):
        dln = log_N[i+1] - log_N[i]
        pw_a.append((log_chi[i+1] - log_chi[i]) / dln)

    print(f"\nq={q_val}:")
    print(f"  Sizes: {[int(n) for n in Ns]}")
    print(f"  Pairwise alpha: {['%.4f' % a for a in pw_a]}")
    if len(Ns) >= 3:
        p = np.polyfit(log_N[-5:], log_chi[-5:], 1)
        print(f"  Global alpha (last 5): {p[0]:.4f}")

    results[f'pairwise_q{q_val}'] = {
        'Ns': [int(n) for n in Ns],
        'alpha': [float(a) for a in pw_a],
    }

save()
print("\nResults saved.")
