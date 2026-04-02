#!/usr/bin/env python3
"""Sprint 082b: Exact diag spin-spin correlator on PERIODIC chain for q=5,7.

Key advantage over 082a (DMRG open BC): on periodic chain, conformal prediction is:
  G(r) = C * [(N/pi) sin(pi*r/N)]^{-2*x_sigma}

This gives clean x_sigma extraction without boundary conformal factors.
Compare x_sigma from correlator to x_1 from gap (same computation).

Also check for oscillatory residuals from complex scaling dimensions.
"""
import numpy as np
import json, time
from scipy.sparse import lil_matrix, csr_matrix
from gpu_utils import eigsh

results = {
    'experiment': '082b_periodic_corr_q5q7',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'data': [],
}

def save():
    with open("results/sprint_082b_periodic_corr_q5q7.json", "w") as f:
        json.dump(results, f, indent=2, default=str)

def build_sq_potts_periodic(n, q, g):
    """Build S_q Potts Hamiltonian on periodic chain of n sites.
    H = -J sum_<ij> delta(s_i, s_j) - g sum_i sum_{k=1}^{q-1} X^k_i
    """
    dim = q**n
    H = lil_matrix((dim, dim), dtype=complex)

    # Helper: state index <-> configuration
    def idx_to_config(idx):
        c = []
        for _ in range(n):
            c.append(idx % q)
            idx //= q
        return c

    def config_to_idx(c):
        idx = 0
        for i in range(n - 1, -1, -1):
            idx = idx * q + c[i]
        return idx

    # Coupling: -J * delta(s_i, s_{i+1}) including periodic bond (n-1, 0)
    for idx in range(dim):
        c = idx_to_config(idx)
        diag = 0.0
        for site in range(n):
            nxt = (site + 1) % n
            if c[site] == c[nxt]:
                diag -= 1.0
        H[idx, idx] += diag

    # Field: -g * sum_i sum_{k=1}^{q-1} X^k_i
    # X|s> = |s+1 mod q>, so X^k|s> = |s+k mod q>
    for idx in range(dim):
        c = idx_to_config(idx)
        for site in range(n):
            for k in range(1, q):
                c2 = c.copy()
                c2[site] = (c[site] + k) % q
                idx2 = config_to_idx(c2)
                H[idx, idx2] += -g

    return csr_matrix(H)

def compute_correlator(psi, n, q):
    """Compute <delta(s_0, s_r)> for all r from ground state vector psi.

    delta(s_i, s_j) = sum_a P_a(i) P_a(j) where P_a = |a><a|
    <delta(s_i, s_j)> = sum_a sum_{configs with s_i=a and s_j=a} |psi(config)|^2
    """
    probs = np.abs(psi)**2  # |psi(config)|^2

    def idx_to_config(idx):
        c = []
        for _ in range(n):
            c.append(idx % q)
            idx //= q
        return c

    # Pre-compute all configurations
    dim = q**n
    configs = np.zeros((dim, n), dtype=int)
    for idx in range(dim):
        c = idx_to_config(idx)
        for site in range(n):
            configs[idx, site] = c[site]

    # For each distance r, compute G(r) = <delta(s_0, s_r)> - 1/q
    # Average over all reference sites (translation invariance on periodic chain)
    G = np.zeros(n)
    counts = np.zeros(n)

    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            r = min(abs(i - j), n - abs(i - j))  # periodic distance
            # delta(s_i, s_j) for each config, weighted by |psi|^2
            same = (configs[:, i] == configs[:, j]).astype(float)
            delta_avg = np.sum(probs * same)
            G[r] += delta_avg
            counts[r] += 1

    # Average
    G_conn = {}
    for r in range(1, n // 2 + 1):
        if counts[r] > 0:
            G_conn[r] = G[r] / counts[r] - 1.0 / q

    return G_conn

def fit_conformal_correlator(G_conn, n):
    """Fit G(r) = C * [(N/pi)*sin(pi*r/N)]^{-eta} to extract eta = 2*x_sigma.
    Uses chord distance on periodic chain.
    """
    r_vals = sorted(G_conn.keys())
    G_vals = [G_conn[r] for r in r_vals]

    # Filter to r >= 2 (short-distance lattice artifacts)
    mask = [r >= 2 for r in r_vals]
    r_fit = np.array([r for r, m in zip(r_vals, mask) if m])
    G_fit = np.array([g for g, m in zip(G_vals, mask) if m])

    if len(r_fit) < 2 or np.any(G_fit <= 0):
        return None

    # Chord distance
    chord = (n / np.pi) * np.sin(np.pi * r_fit / n)

    log_chord = np.log(chord)
    log_G = np.log(G_fit)

    A = np.vstack([log_chord, np.ones_like(log_chord)]).T
    (slope, intercept), _, _, _ = np.linalg.lstsq(A, log_G, rcond=None)

    eta = -slope
    x_sigma = eta / 2

    # R² and residuals
    G_pred = np.exp(intercept) * chord**(-eta)
    ss_res = np.sum((G_fit - G_pred)**2)
    ss_tot = np.sum((G_fit - np.mean(G_fit))**2)
    R2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0

    frac_resid = (G_fit - G_pred) / G_pred
    osc_amp = float(np.std(frac_resid))

    return {
        'eta': float(eta),
        'x_sigma': float(x_sigma),
        'R2': float(R2),
        'C': float(np.exp(intercept)),
        'oscillation_amplitude': osc_amp,
        'max_abs_resid': float(np.max(np.abs(frac_resid))),
        'r_fit': r_fit.tolist(),
        'G_fit': G_fit.tolist(),
        'G_pred': G_pred.tolist(),
        'frac_resid': frac_resid.tolist(),
    }

# Run for q=5 and q=7
cases = [
    (5, 8, "q=5 n=8 (dim=390k, walking)"),
    (7, 6, "q=7 n=6 (dim=118k, breaking)"),
    (7, 7, "q=7 n=7 (dim=824k, breaking)"),
]

print("Sprint 082b: Periodic chain correlator for S_q Potts")
print("=" * 60, flush=True)

for q, n, label in cases:
    gc = 1.0 / q
    dim = q**n
    print(f"\n--- {label}, g_c=1/{q}={gc:.4f}, dim={dim:,} ---", flush=True)

    # Build Hamiltonian
    t0 = time.time()
    H = build_sq_potts_periodic(n, q, gc)
    t_build = time.time() - t0
    print(f"  H built in {t_build:.1f}s", flush=True)

    # Get ground state + first excited (for gap comparison)
    t0 = time.time()
    evals, evecs = eigsh(H, k=4, which='SA')
    t_eig = time.time() - t0

    order = np.argsort(evals)
    evals = evals[order]
    evecs = evecs[:, order]

    E0 = float(evals[0])
    gap = float(evals[1] - evals[0])
    gap_N = gap * n
    print(f"  eigsh in {t_eig:.1f}s, E0={E0:.6f}, gap={gap:.6f}, gap×N={gap_N:.4f}", flush=True)

    # Compute correlator from ground state
    t0 = time.time()
    psi = evecs[:, 0].real  # Ground state should be real (up to phase)
    # Check if ground state has imaginary part
    if np.max(np.abs(evecs[:, 0].imag)) > 1e-10:
        print(f"  WARNING: ground state has significant imaginary part", flush=True)
        psi = evecs[:, 0]  # Use complex

    G_conn = compute_correlator(psi, n, q)
    t_corr = time.time() - t0
    print(f"  Correlator in {t_corr:.1f}s", flush=True)

    # Show correlator values
    print(f"  r   G_conn(r)    chord_dist", flush=True)
    for r in sorted(G_conn.keys()):
        chord = (n / np.pi) * np.sin(np.pi * r / n)
        print(f"  {r:2d}   {G_conn[r]:12.6f}   {chord:.3f}", flush=True)

    # Fit conformal form
    fit = fit_conformal_correlator(G_conn, n)
    if fit:
        print(f"  Conformal fit: eta={fit['eta']:.4f}, x_sigma={fit['x_sigma']:.4f}, R²={fit['R2']:.6f}", flush=True)
        print(f"  Oscillation amp: {fit['oscillation_amplitude']:.4f}, max|resid|: {fit['max_abs_resid']:.4f}", flush=True)

        # Compare to gap
        x1_gap = gap_N / (2 * np.pi)  # Assuming v=1; we'll correct below
        print(f"  x_1(gap, v=1): {x1_gap:.4f}", flush=True)

        # Get velocity from Casimir energy: E0/N = e_inf - pi*v*c/(6*N^2)
        # Need two sizes for v extraction; for now just report x_sigma

    # First excited degeneracy
    degen_gap = float(evals[2] - evals[1])
    print(f"  E1 degeneracy check: E2-E1={degen_gap:.6f} (should be ~0 if degenerate)", flush=True)

    entry = {
        'q': q, 'n': n, 'gc': gc, 'dim': dim,
        'E0': E0, 'gap': gap, 'gap_N': gap_N,
        'G_conn': {str(k): v for k, v in G_conn.items()},
        'conformal_fit': fit,
        'time_build': t_build, 'time_eig': t_eig, 'time_corr': t_corr,
    }
    results['data'].append(entry)
    save()
    print(f"  Saved.", flush=True)

    total_time = t_build + t_eig + t_corr
    if total_time > 200:
        print(f"  Total {total_time:.0f}s — approaching limit", flush=True)

# Summary comparison
print(f"\n{'='*60}")
print("SUMMARY: Conformal correlator analysis")
print(f"{'='*60}")
print(f"{'q':>3} {'n':>3} {'eta':>8} {'x_sigma':>8} {'R²':>10} {'osc_amp':>10} {'gap×N':>8}")
for d in results['data']:
    f = d['conformal_fit']
    if f:
        print(f"{d['q']:3d} {d['n']:3d} {f['eta']:8.4f} {f['x_sigma']:8.4f} {f['R2']:10.6f} {f['oscillation_amplitude']:10.4f} {d['gap_N']:8.4f}")
    else:
        print(f"{d['q']:3d} {d['n']:3d}      —        —          —          — {d['gap_N']:8.4f}")

# Compare q=5 vs q=7
q5_data = [d for d in results['data'] if d['q'] == 5]
q7_data = [d for d in results['data'] if d['q'] == 7]
if q5_data and q7_data:
    f5 = q5_data[0]['conformal_fit']
    f7 = q7_data[-1]['conformal_fit']
    if f5 and f7:
        print(f"\nq=5 vs q=7 comparison:")
        print(f"  x_sigma: q=5={f5['x_sigma']:.4f}, q=7={f7['x_sigma']:.4f}")
        print(f"  R²: q=5={f5['R2']:.6f}, q=7={f7['R2']:.6f}")
        print(f"  Oscillation: q=5={f5['oscillation_amplitude']:.4f}, q=7={f7['oscillation_amplitude']:.4f}")

save()
print(f"\nFinal save to results/sprint_082b_periodic_corr_q5q7.json")

from db_utils import record
for d in results['data']:
    f = d['conformal_fit']
    if f and f['x_sigma'] is not None:
        record(sprint=82, model='sq_potts', q=d['q'], n=d['n'],
               quantity='x_sigma_corr', value=f['x_sigma'],
               method='periodic_correlator',
               notes=f'gc=1/{d["q"]}, eta={f["eta"]:.4f}, R2={f["R2"]:.4f}')
        record(sprint=82, model='sq_potts', q=d['q'], n=d['n'],
               quantity='corr_osc_amp', value=f['oscillation_amplitude'],
               method='periodic_correlator',
               notes=f'fractional residual std from power law fit')
print("Recorded to DB.")
