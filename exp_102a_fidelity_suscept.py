#!/usr/bin/env python3
"""Sprint 102a: Fidelity susceptibility scan near g_c = 1/q for S_q Potts chains.

chi_F(g) = (2/N) * (1 - |<psi(g)|psi(g+dg)>|) / dg^2

Optimization: Build coupling and field matrices separately, then H(g) = H_coup + g*H_field.
This avoids rebuilding the full Hamiltonian for each g value.
"""
import numpy as np
import json, time
from scipy.sparse import lil_matrix, csr_matrix
from gpu_utils import eigsh

results = {
    'experiment': '102a_fidelity_suscept',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'data': {},
}

def save():
    with open("results/sprint_102a_fidelity_suscept.json", "w") as f:
        json.dump(results, f, indent=2, default=str)

def build_sq_potts_parts(n, q):
    """Build coupling and field parts separately: H(g) = H_coup + g * H_field."""
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

def ground_state(H_coup, H_field, g):
    """Get ground state for H = H_coup + g * H_field."""
    H = H_coup + g * H_field
    _, vecs = eigsh(H, k=1, which='SA')
    return vecs[:, 0]

# Configuration: sizes and g-scan ranges
# Keep sizes moderate to fit in time budget
configs = {
    2: {'sizes': [6, 8, 10, 12], 'g_c': 0.5, 'g_range': 0.15, 'ng': 21},
    3: {'sizes': [6, 8, 10], 'g_c': 1/3, 'g_range': 0.10, 'ng': 19},
    5: {'sizes': [6, 8], 'g_c': 0.2, 'g_range': 0.06, 'ng': 17},
    7: {'sizes': [6], 'g_c': 1/7, 'g_range': 0.04, 'ng': 15},
}

dg = 1e-4  # finite difference step

print("=" * 60)
print("Sprint 102a: Fidelity susceptibility scan")
print("=" * 60)

for q in [2, 3, 5, 7]:
    cfg = configs[q]
    g_c = cfg['g_c']
    g_vals = np.linspace(g_c - cfg['g_range'], g_c + cfg['g_range'], cfg['ng'])
    results['data'][str(q)] = {'g_c': g_c, 'sizes': {}}

    for n in cfg['sizes']:
        dim = q**n
        print(f"\nq={q}, n={n} (dim={dim})...", flush=True)
        t0 = time.time()

        # Build once
        H_coup, H_field = build_sq_potts_parts(n, q)
        t_build = time.time() - t0
        print(f"  Built H parts in {t_build:.1f}s", flush=True)

        chi_vals = []
        for g in g_vals:
            psi1 = ground_state(H_coup, H_field, g)
            psi2 = ground_state(H_coup, H_field, g + dg)
            overlap = abs(np.vdot(psi1, psi2))
            overlap = min(overlap, 1.0 - 1e-15)
            chi_F = 2.0 * (1.0 - overlap) / (dg**2 * n)
            chi_vals.append(float(chi_F))

        dt = time.time() - t0

        # Find peak
        i_peak = np.argmax(chi_vals)
        g_peak = g_vals[i_peak]
        chi_peak = chi_vals[i_peak]

        # Estimate FWHM
        half_max = chi_peak / 2.0
        above = np.array(chi_vals) > half_max
        if above.any():
            left = np.where(above)[0][0]
            right = np.where(above)[0][-1]
            width = g_vals[right] - g_vals[left]
        else:
            width = float('nan')

        results['data'][str(q)]['sizes'][str(n)] = {
            'dim': dim,
            'g_vals': g_vals.tolist(),
            'chi_F': chi_vals,
            'g_peak': float(g_peak),
            'chi_peak': float(chi_peak),
            'width_FWHM': float(width),
            'time_s': dt,
        }

        print(f"  Total time: {dt:.1f}s")
        print(f"  g_peak={g_peak:.4f} (g_c={g_c:.4f}, shift={g_peak-g_c:+.4f})")
        print(f"  chi_peak={chi_peak:.2f}, FWHM={width:.4f}")

        save()

        # Time guard: skip if next size would exceed budget
        if dt > 50:
            print(f"  Skipping larger sizes (last took {dt:.0f}s)")
            break

# Summary table
print("\n" + "=" * 60)
print("SUMMARY: Peak fidelity susceptibility")
print(f"{'q':>3} {'n':>4} {'g_peak':>8} {'g_c':>8} {'shift':>8} {'chi_max':>10} {'FWHM':>8}")
print("-" * 60)
for q in [2, 3, 5, 7]:
    g_c = configs[q]['g_c']
    for n_str, d in results['data'][str(q)]['sizes'].items():
        n = int(n_str)
        print(f"{q:>3} {n:>4} {d['g_peak']:>8.4f} {g_c:>8.4f} {d['g_peak']-g_c:>+8.4f} {d['chi_peak']:>10.2f} {d['width_FWHM']:>8.4f}")

save()
print("\nResults saved to results/sprint_102a_fidelity_suscept.json")
