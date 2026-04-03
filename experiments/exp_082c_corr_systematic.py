#!/usr/bin/env python3
"""Sprint 082c: Systematic correlator across q=2-8 on largest feasible periodic chains.

Goal: More r-points to detect oscillatory corrections from complex exponents.
Also extract x_sigma(q) curve and compare to gap-derived x_1(q).

Sizes chosen for ~60s per case:
  q=2 n=14 (dim=16k), q=3 n=10 (dim=59k), q=4 n=8 (dim=65k),
  q=5 n=8 (dim=390k), q=6 n=7 (dim=280k), q=7 n=7 (dim=824k),
  q=8 n=6 (dim=262k)
"""
import numpy as np
import json, time
from scipy.sparse import lil_matrix, csr_matrix
from gpu_utils import eigsh

results = {
    'experiment': '082c_corr_systematic',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'data': [],
}

def save():
    with open("results/sprint_082c_corr_systematic.json", "w") as f:
        json.dump(results, f, indent=2, default=str)

def build_sq_potts_periodic(n, q, g):
    """Build S_q Potts Hamiltonian on periodic chain."""
    dim = q**n
    H = lil_matrix((dim, dim), dtype=complex)

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

    for idx in range(dim):
        c = idx_to_config(idx)
        diag = 0.0
        for site in range(n):
            nxt = (site + 1) % n
            if c[site] == c[nxt]:
                diag -= 1.0
        H[idx, idx] += diag

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
    """Compute connected correlator G(r) = <delta(s_0,s_r)> - 1/q, averaged over sites."""
    probs = np.abs(psi)**2
    dim = q**n

    configs = np.zeros((dim, n), dtype=int)
    for idx in range(dim):
        tmp = idx
        for site in range(n):
            configs[idx, site] = tmp % q
            tmp //= q

    G = np.zeros(n)
    counts = np.zeros(n)

    for i in range(n):
        for j in range(i + 1, n):
            r = min(j - i, n - (j - i))
            same = (configs[:, i] == configs[:, j]).astype(float)
            delta_avg = np.sum(probs * same)
            G[r] += delta_avg
            counts[r] += 1

    G_conn = {}
    for r in range(1, n // 2 + 1):
        if counts[r] > 0:
            G_conn[r] = G[r] / counts[r] - 1.0 / q

    return G_conn

def fit_conformal(G_conn, n, min_r=2):
    """Fit G(r) = C * [(N/pi)*sin(pi*r/N)]^{-eta} on periodic chain."""
    r_vals = sorted([r for r in G_conn.keys() if r >= min_r])
    G_vals = np.array([G_conn[r] for r in r_vals])

    if len(r_vals) < 2 or np.any(G_vals <= 0):
        return None

    chord = (n / np.pi) * np.sin(np.pi * np.array(r_vals, float) / n)
    log_chord = np.log(chord)
    log_G = np.log(G_vals)

    A = np.vstack([log_chord, np.ones_like(log_chord)]).T
    (slope, intercept), _, _, _ = np.linalg.lstsq(A, log_G, rcond=None)

    eta = -slope
    x_sigma = eta / 2

    G_pred = np.exp(intercept) * chord**(-eta)
    ss_res = np.sum((G_vals - G_pred)**2)
    ss_tot = np.sum((G_vals - np.mean(G_vals))**2)
    R2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0

    frac_resid = (G_vals - G_pred) / G_pred
    osc_amp = float(np.std(frac_resid))

    return {
        'eta': float(eta),
        'x_sigma': float(x_sigma),
        'R2': float(R2),
        'oscillation_amplitude': osc_amp,
        'max_abs_resid': float(np.max(np.abs(frac_resid))),
        'n_fit_points': len(r_vals),
        'r_vals': [int(r) for r in r_vals],
        'G_vals': G_vals.tolist(),
        'G_pred': G_pred.tolist(),
        'frac_resid': frac_resid.tolist(),
    }

# Cases: (q, n) chosen for largest feasible periodic chain
cases = [
    (2, 14, "Ising reference"),
    (3, 10, "3-state Potts"),
    (4, 8,  "Ashkin-Teller marginal"),
    (5, 8,  "walking sweet spot"),
    (6, 7,  "marginal breaking"),
    (7, 7,  "walking broken"),
    (8, 6,  "walking broken"),
]

# Complex CFT predictions
def complex_cft_x_sigma(q):
    """Compute Re(x_sigma) from complex CFT Coulomb gas."""
    if q <= 4:
        return None  # Real CFT, use exact values
    sqrt_q = np.sqrt(q)
    alpha = np.arccosh(sqrt_q / 2)
    p_complex = 1j * np.pi / alpha  # complex p
    # x_sigma from Kac table: x_{1,2} = ((p-1)(p+2))/(8p) for spin field
    # Actually for Potts, x_sigma = x_{(1,0)} or similar. Let me use the standard:
    # For q-state Potts, the spin field dimension is:
    # h_sigma = (p-1)/(2p) for q>4 complex
    # Actually, from Gorbenko et al., the magnetic dimension is complex.
    # Re(h) ≈ (5p_R^2 + p_I^2 - 6p_R)/(8(p_R^2 + p_I^2)) ... too complex.
    # Just return None and compare empirically.
    return None

print("Sprint 082c: Systematic correlator across q=2-8 (S_q Potts, periodic)")
print("=" * 70, flush=True)

# Known exact values for comparison
exact_x_sigma = {2: 0.125, 3: 2.0/15.0}  # Ising: 1/8, 3-state: 2/15

for q, n, label in cases:
    gc = 1.0 / q
    dim = q**n
    print(f"\n--- q={q} n={n} ({label}), dim={dim:,} ---", flush=True)

    t0 = time.time()
    H = build_sq_potts_periodic(n, q, gc)
    t_build = time.time() - t0

    t0 = time.time()
    evals, evecs = eigsh(H, k=4, which='SA')
    t_eig = time.time() - t0

    order = np.argsort(evals)
    evals = evals[order]
    evecs = evecs[:, order]

    E0 = float(evals[0])
    gap = float(evals[1] - evals[0])
    gap_N = gap * n
    print(f"  Build {t_build:.1f}s, eigsh {t_eig:.1f}s, gap×N={gap_N:.4f}", flush=True)

    # Degeneracy of first excited
    degen = float(evals[2] - evals[1])
    print(f"  1st excited (q-1)-fold degen: E2-E1={degen:.2e}", flush=True)

    t0 = time.time()
    psi = evecs[:, 0]
    G_conn = compute_correlator(psi, n, q)
    t_corr = time.time() - t0

    # Show correlator
    print(f"  r    G(r)         chord", flush=True)
    for r in sorted(G_conn.keys()):
        chord = (n / np.pi) * np.sin(np.pi * r / n)
        print(f"  {r:2d}   {G_conn[r]:12.6f}   {chord:.3f}", flush=True)

    # Fit conformal form
    fit = fit_conformal(G_conn, n, min_r=2)
    if fit:
        exact_str = ""
        if q in exact_x_sigma:
            exact = exact_x_sigma[q]
            dev_pct = 100 * (fit['x_sigma'] - exact) / exact
            exact_str = f" (exact={exact:.4f}, dev={dev_pct:+.1f}%)"
        print(f"  Fit: eta={fit['eta']:.4f}, x_sigma={fit['x_sigma']:.4f}{exact_str}", flush=True)
        print(f"  R²={fit['R2']:.6f}, osc_amp={fit['oscillation_amplitude']:.4f}, "
              f"max|resid|={fit['max_abs_resid']:.4f}, n_pts={fit['n_fit_points']}", flush=True)

        # Show per-point residuals
        if fit['n_fit_points'] > 2:
            print(f"  Per-point fractional residuals:", flush=True)
            for r, fr in zip(fit['r_vals'], fit['frac_resid']):
                print(f"    r={r}: {fr:+.6f}", flush=True)

    entry = {
        'q': q, 'n': n, 'gc': gc, 'dim': dim,
        'E0': E0, 'gap': gap, 'gap_N': gap_N,
        'first_excited_degen': degen,
        'G_conn': {str(k): v for k, v in G_conn.items()},
        'conformal_fit': fit,
        'time_total': t_build + t_eig + t_corr,
    }
    results['data'].append(entry)
    save()

    total = t_build + t_eig + t_corr
    if total > 200:
        print(f"  WARNING: {total:.0f}s total — may need to reduce sizes for remaining q", flush=True)

# Grand summary
print(f"\n{'='*70}")
print("GRAND SUMMARY: x_sigma(q) from periodic chain correlator")
print(f"{'='*70}")
print(f"{'q':>3} {'n':>3} {'x_sigma':>8} {'exact':>8} {'dev%':>7} {'R²':>10} {'osc_amp':>8} {'n_pts':>5} {'gap×N':>7}")
for d in results['data']:
    f = d['conformal_fit']
    q = d['q']
    if f:
        exact = exact_x_sigma.get(q, None)
        exact_str = f"{exact:.4f}" if exact else "—"
        dev_str = f"{100*(f['x_sigma']-exact)/exact:+.1f}" if exact else "—"
        print(f"{q:3d} {d['n']:3d} {f['x_sigma']:8.4f} {exact_str:>8} {dev_str:>7} "
              f"{f['R2']:10.6f} {f['oscillation_amplitude']:8.4f} {f['n_fit_points']:5d} {d['gap_N']:7.4f}")

# x_sigma trend analysis
print(f"\nx_sigma(q) trend:")
xs = [(d['q'], d['conformal_fit']['x_sigma']) for d in results['data'] if d['conformal_fit']]
for q, x in xs:
    bar = '#' * int(x * 200)
    print(f"  q={q}: {x:.4f} {bar}")

# Check for walking signature: does oscillation amplitude increase with q?
print(f"\nOscillation amplitude vs q (walking signature):")
for d in results['data']:
    f = d['conformal_fit']
    if f and f['n_fit_points'] > 2:
        print(f"  q={d['q']}: osc_amp={f['oscillation_amplitude']:.6f}, max|resid|={f['max_abs_resid']:.6f}")

save()
print(f"\nFinal save to results/sprint_082c_corr_systematic.json")

from db_utils import record
for d in results['data']:
    f = d['conformal_fit']
    if f:
        record(sprint=82, model='sq_potts', q=d['q'], n=d['n'],
               quantity='x_sigma_corr', value=f['x_sigma'],
               method='periodic_correlator',
               notes=f'gc=1/{d["q"]}, eta={f["eta"]:.4f}, R2={f["R2"]:.4f}, n_pts={f["n_fit_points"]}')
print("Recorded to DB.")
