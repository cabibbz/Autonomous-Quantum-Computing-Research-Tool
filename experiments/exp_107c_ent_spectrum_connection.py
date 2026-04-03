"""Sprint 107c: Connect chi_F multiplet gap to entanglement spectrum multiplet.

Sprint 106: chi_F dominated by (q-1)-fold energy multiplet with gap ~ N^{-z_m}
Sprint 084: Entanglement spectrum has (q-1)-fold multiplet with gap Delta_xi

Question: Are the energy multiplet gap and entanglement gap related?
- Does Delta_xi scale the same as Delta_m with N?
- Does the entropy concentration in level 1 correlate with beta_me?

Compute both at same sizes for q=2,3,5,7.
"""
import numpy as np
import json, time, os
from scipy.sparse import csr_matrix, coo_matrix
from scipy.linalg import eigh
from gpu_utils import eigsh
from db_utils import record

results = {
    'experiment': '107c_ent_spectrum_connection',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'data': {},
}

def save():
    outpath = os.path.join(os.path.dirname(__file__), 'results', 'sprint_107c_ent_spectrum_connection.json')
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    with open(outpath, 'w') as f:
        json.dump(results, f, indent=2, default=str)

def build_sq_potts_parts_fast(n, q):
    """Vectorized S_q Potts Hamiltonian builder."""
    dim = q**n
    all_idx = np.arange(dim, dtype=np.int64)
    digits = np.zeros((dim, n), dtype=np.int64)
    tmp = all_idx.copy()
    for site in range(n):
        digits[:, site] = tmp % q
        tmp //= q
    powers = q ** np.arange(n, dtype=np.int64)

    diag_vals = np.zeros(dim, dtype=np.float64)
    for site in range(n):
        nxt = (site + 1) % n
        diag_vals -= (digits[:, site] == digits[:, nxt]).astype(np.float64)
    H_coup = csr_matrix((diag_vals, (all_idx, all_idx)), shape=(dim, dim))

    rows_list, cols_list, vals_list = [], [], []
    for site in range(n):
        pw = powers[site]
        old_digit = digits[:, site]
        for k in range(1, q):
            new_digit = (old_digit + k) % q
            delta = (new_digit.astype(np.int64) - old_digit.astype(np.int64)) * pw
            new_idx = all_idx + delta
            rows_list.append(new_idx)
            cols_list.append(all_idx)
            vals_list.append(np.full(dim, -1.0, dtype=np.float64))
    rows = np.concatenate(rows_list)
    cols = np.concatenate(cols_list)
    vals = np.concatenate(vals_list)
    H_field = csr_matrix(coo_matrix((vals, (rows, cols)), shape=(dim, dim)))
    return H_coup, H_field

def get_entanglement_spectrum(psi, n, q, nA):
    """Get entanglement spectrum for half-chain bipartition."""
    dim_A = q**nA
    dim_B = q**(n - nA)
    psi_mat = psi.reshape(dim_B, dim_A)  # reshape as |B>|A>
    # Reduced density matrix for A
    rho_A = psi_mat.T.conj() @ psi_mat
    eigenvalues = np.linalg.eigvalsh(rho_A)
    eigenvalues = eigenvalues[::-1]  # descending
    eigenvalues = eigenvalues[eigenvalues > 1e-14]
    return eigenvalues

configs = [
    (2, [6, 8, 10, 12]),
    (3, [6, 8, 10]),
    (5, [6, 7, 8]),
    (7, [6, 7]),
]

print("=" * 70)
print("Sprint 107c: Energy multiplet gap vs entanglement gap")
print("=" * 70)

for q, sizes in configs:
    g_c = 1.0 / q
    k_need = min(q + 3, 12)
    results['data'][f'q{q}'] = {'q': q, 'sizes': {}}

    energy_gaps = []
    ent_gaps = []
    me_sqs = []
    ent_S1_fracs = []
    Ns = []

    for n in sizes:
        dim = q**n
        nA = n // 2
        print(f"\nq={q}, n={n} (dim={dim:,}, nA={nA})...", flush=True)
        t0 = time.time()

        H_coup, H_field = build_sq_potts_parts_fast(n, q)
        H = H_coup + g_c * H_field

        k_use = min(k_need, dim - 2)
        evals, evecs = eigsh(H, k=k_use, which='SA')
        order = np.argsort(evals)
        evals = evals[order]
        evecs = evecs[:, order]

        E0 = evals[0]
        psi0 = evecs[:, 0]
        H_field_psi0 = H_field.dot(psi0)

        # Energy multiplet: find dominant chi_F contributor
        best_gap = 0
        best_me_sq = 0
        for i in range(1, k_use):
            gap = evals[i] - E0
            if gap < 1e-12:
                continue
            me = np.dot(evecs[:, i], H_field_psi0)
            chi_n = me**2 / gap**2
            if chi_n > best_me_sq / (best_gap**2 if best_gap > 0 else 1):
                best_gap = gap
                best_me_sq = me**2

        # Entanglement spectrum
        ent_evals = get_entanglement_spectrum(psi0, n, q, nA)
        ent_energies = -np.log(ent_evals)
        ent_gap = ent_energies[1] - ent_energies[0] if len(ent_energies) > 1 else 0

        # Entropy decomposition by level
        total_S = -np.sum(ent_evals * np.log(ent_evals + 1e-30))
        # Level 0 = ground state eigenvalue
        s0 = -ent_evals[0] * np.log(ent_evals[0] + 1e-30)
        # Level 1 = (q-1)-fold multiplet (next q-1 eigenvalues)
        n_lev1 = min(q - 1, len(ent_evals) - 1)
        s1 = 0
        for j in range(1, 1 + n_lev1):
            if j < len(ent_evals):
                s1 += -ent_evals[j] * np.log(ent_evals[j] + 1e-30)
        s1_frac = s1 / total_S if total_S > 0 else 0

        dt = time.time() - t0

        energy_gaps.append(best_gap)
        ent_gaps.append(ent_gap)
        me_sqs.append(best_me_sq)
        ent_S1_fracs.append(s1_frac)
        Ns.append(n)

        results['data'][f'q{q}']['sizes'][str(n)] = {
            'dim': dim, 'nA': nA,
            'energy_multiplet_gap': float(best_gap),
            'entanglement_gap': float(ent_gap),
            'me_sq': float(best_me_sq),
            'S_total': float(total_S),
            'S_lev0_frac': float(s0 / total_S if total_S > 0 else 0),
            'S_lev1_frac': float(s1_frac),
            'time_s': dt,
        }

        record(sprint=107, model='S_q_Potts', q=q, n=n,
               quantity='ent_gap', value=ent_gap, method='half_chain')

        print(f"  Energy gap_m = {best_gap:.6f}  |me|² = {best_me_sq:.4f}")
        print(f"  Ent gap Δξ  = {ent_gap:.6f}")
        print(f"  S_total = {total_S:.4f}  S_lev1_frac = {s1_frac:.4f}")
        print(f"  ({dt:.1f}s)")

        save()

    # Scaling comparison
    if len(Ns) >= 3:
        log_N = np.log(np.array(Ns, dtype=float))
        z_energy = -np.polyfit(log_N, np.log(np.array(energy_gaps)), 1)[0]
        z_ent = -np.polyfit(log_N, np.log(np.array(ent_gaps)), 1)[0] if all(g > 0 for g in ent_gaps) else float('nan')
        beta_me = np.polyfit(log_N, np.log(np.array(me_sqs)), 1)[0]

        # Ratio of gaps
        ratios = [eg / mg for eg, mg in zip(ent_gaps, energy_gaps)]

        print(f"\n  q={q} SCALING:")
        print(f"    z_energy (gap_m ~ N^{{-z}}) = {z_energy:.4f}")
        print(f"    z_ent (Δξ ~ N^{{-z}})       = {z_ent:.4f}")
        print(f"    beta_me                     = {beta_me:.4f}")
        print(f"    Δξ/gap_m ratios: {['%.3f' % r for r in ratios]}")
        print(f"    S_lev1 fracs:    {['%.4f' % f for f in ent_S1_fracs]}")

        results['data'][f'q{q}']['z_energy'] = float(z_energy)
        results['data'][f'q{q}']['z_ent'] = float(z_ent)
        results['data'][f'q{q}']['beta_me'] = float(beta_me)
        results['data'][f'q{q}']['gap_ratios'] = [float(r) for r in ratios]

        save()

# Summary table
print(f"\n{'=' * 70}")
print("SUMMARY: Energy vs Entanglement Multiplet Gap Scaling")
print(f"{'q':>3} {'#pts':>5} {'z_energy':>10} {'z_ent':>10} {'ratio':>8}")
print("-" * 45)
for q, _ in configs:
    d = results['data'].get(f'q{q}', {})
    if 'z_energy' in d:
        npts = len(d.get('sizes', {}))
        print(f"{q:>3} {npts:>5} {d['z_energy']:>10.4f} {d['z_ent']:>10.4f} {d['z_ent']/d['z_energy']:>8.3f}")

save()
print("\nResults saved.")
