"""Sprint 108c: Multiplet structure comparison at q=3, 4, 5.

At q=4, the dominant multiplet is (q-1)=3-fold degenerate. But q=4 is special:
it has both S_4 permutation symmetry AND is at the BKT boundary.

Questions:
1. How does the multiplet gap degeneracy pattern compare across q=3,4,5?
2. Is the (q-1)-fold degeneracy exact or split at q=4?
3. How does the spectral gap (level 1) compare to the multiplet gap?
4. What is the gap ratio gap_multiplet/gap_spectral vs N?

This reveals whether q=4 multiplet structure is BKT-like or walking-like.
"""
import numpy as np
import json, time, os
from scipy.sparse import lil_matrix, csr_matrix
from gpu_utils import eigsh
from db_utils import record

results = {
    'experiment': '108c_multiplet_structure',
    'sprint': 108,
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'data': {},
}

def save():
    outpath = os.path.join(os.path.dirname(__file__), 'results', 'sprint_108c_multiplet_structure.json')
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    with open(outpath, 'w') as f:
        json.dump(results, f, indent=2, default=str)

def build_sq_potts_parts(n, q):
    """Build coupling and field parts separately."""
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

print("=" * 70)
print("Sprint 108c: Multiplet structure at q=3,4,5")
print("=" * 70)

total_time = 0

for q in [3, 4, 5]:
    g_c = 1.0 / q
    # Use enough states to see the full low-energy structure
    k_states = min(3*q, 20)
    sizes = [6, 8]
    if q <= 4:
        sizes.append(10)  # q=3: dim=59049, q=4: dim=1M
    if q == 3:
        sizes.append(12)  # dim=531441

    results['data'][f'q{q}'] = {'q': q, 'g_c': g_c, 'sizes': {}}

    print(f"\n{'='*50}")
    print(f"q={q}, g_c={g_c:.4f}, k_states={k_states}")
    print(f"{'='*50}")

    for n in sizes:
        dim = q**n
        if dim > 2_000_000:
            print(f"\n  n={n}: dim={dim:,} — skipping (too large for k={k_states})")
            continue

        print(f"\n  n={n} (dim={dim:,})...", flush=True)
        t0 = time.time()

        H_coup, H_field = build_sq_potts_parts(n, q)
        H = H_coup + g_c * H_field

        k_use = min(k_states, dim - 2)
        evals, evecs = eigsh(H, k=k_use, which='SA')
        order = np.argsort(evals)
        evals = evals[order]
        evecs = evecs[:, order]

        E0 = evals[0]
        psi0 = evecs[:, 0]
        H_field_psi0 = H_field.dot(psi0)

        # Gap structure
        gaps = evals[1:] - E0
        print(f"  Energy gaps (first {min(k_use-1, 2*q)} levels):")
        for i, g in enumerate(gaps[:2*q]):
            # Matrix element with field
            me = np.dot(evecs[:, i+1], H_field_psi0)
            me_sq = me**2
            print(f"    level {i+1}: gap={g:.8f}  |me|²={me_sq:.6f}")

        # Identify degeneracies (gaps within 1e-6 of each other)
        degen_groups = []
        used = set()
        for i in range(len(gaps)):
            if i in used:
                continue
            group = [i]
            used.add(i)
            for j in range(i+1, len(gaps)):
                if j in used:
                    continue
                if abs(gaps[j] - gaps[i]) / max(abs(gaps[i]), 1e-15) < 1e-4:
                    group.append(j)
                    used.add(j)
            degen_groups.append(group)

        print(f"\n  Degeneracy pattern (relative tol=1e-4):")
        for grp in degen_groups[:6]:
            deg = len(grp)
            gap_val = gaps[grp[0]]
            # Matrix elements for this group
            mes = []
            for idx in grp:
                me = np.dot(evecs[:, idx+1], H_field_psi0)
                mes.append(me**2)
            total_me = sum(mes)
            print(f"    gap={gap_val:.8f}  degeneracy={deg}  Σ|me|²={total_me:.6f}")

        # Gap ratio: multiplet/spectral
        gap_spectral = gaps[0]
        # Find the dominant multiplet (largest total |me|^2)
        best_group = None
        best_me_total = 0
        for grp in degen_groups:
            total_me = sum(np.dot(evecs[:, idx+1], H_field_psi0)**2 for idx in grp)
            if total_me > best_me_total:
                best_me_total = total_me
                best_group = grp

        gap_multiplet = gaps[best_group[0]] if best_group else gaps[0]
        gap_ratio = gap_multiplet / gap_spectral if gap_spectral > 1e-15 else float('inf')
        multiplet_deg = len(best_group) if best_group else 0

        dt = time.time() - t0
        total_time += dt

        print(f"\n  Gap ratio (multiplet/spectral): {gap_ratio:.4f}")
        print(f"  Multiplet degeneracy: {multiplet_deg}")
        print(f"  Time: {dt:.1f}s")

        results['data'][f'q{q}']['sizes'][str(n)] = {
            'n': n, 'dim': dim,
            'gaps': [float(g) for g in gaps[:2*q]],
            'degeneracy_groups': [[int(i) for i in grp] for grp in degen_groups[:8]],
            'gap_spectral': float(gap_spectral),
            'gap_multiplet': float(gap_multiplet),
            'gap_ratio': float(gap_ratio),
            'multiplet_degeneracy': multiplet_deg,
            'multiplet_me_sq_total': float(best_me_total),
            'time_s': dt,
        }

        record(sprint=108, model='S_q_Potts', q=q, n=n,
               quantity='gap_ratio', value=gap_ratio, method='spectral_structure')

        save()

        if total_time > 200:
            print("\n  Total time approaching limit, stopping")
            break

    if total_time > 200:
        break

# Summary
print(f"\n{'=' * 70}")
print("SUMMARY: Multiplet structure comparison")
print(f"{'q':>3} {'n':>4} {'gap_spec':>10} {'gap_mult':>10} {'ratio':>8} {'deg':>5}")
print("-" * 45)
for q in [3, 4, 5]:
    qkey = f'q{q}'
    if qkey in results['data']:
        for nk in sorted(results['data'][qkey]['sizes'].keys(), key=int):
            d = results['data'][qkey]['sizes'][nk]
            print(f"{q:>3} {d['n']:>4} {d['gap_spectral']:>10.6f} "
                  f"{d['gap_multiplet']:>10.6f} {d['gap_ratio']:>8.4f} "
                  f"{d['multiplet_degeneracy']:>5}")

save()
print(f"\nTotal time: {total_time:.1f}s")
print("Results saved.")
