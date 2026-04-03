#!/usr/bin/env python3
"""Sprint 049b: Correlation length ξ(g) for q=2,3,4,5 Potts + FSS collapse for ν.

Extract ξ from exponential decay of ⟨P0_i P0_j⟩_conn for Potts (q≥3),
or ⟨σx_i σx_j⟩_conn for TFIM (q=2). Multiple system sizes enable
ξ/n FSS collapse: ξ(g,n)/n = f((g-g_c)·n^{1/ν}).

For BKT (q=4?): ξ ~ exp(c/|g-g_c|^{1/2}) — log(ξ) vs log|g-g_c| curves upward.
For power-law: ξ ~ |g-g_c|^{-ν} — log(ξ) vs log|g-g_c| is a straight line.
"""
import numpy as np
import numpy.linalg as la
import json, time
from scipy.optimize import minimize_scalar
import warnings
warnings.filterwarnings('ignore')

from tenpy.models.tf_ising import TFIChain
from tenpy.models.model import CouplingMPOModel, NearestNeighborModel
from tenpy.networks.site import Site
from tenpy.linalg import np_conserved as npc
from tenpy.networks.mps import MPS
from tenpy.algorithms import dmrg


class PottsSite(Site):
    def __init__(self, q):
        leg = npc.LegCharge.from_trivial(q)
        Site.__init__(self, leg, [str(a) for a in range(q)], sort_charge=False)
        for a in range(q):
            P = np.zeros((q, q), dtype=complex); P[a, a] = 1.0
            self.add_op(f'P{a}', P)
        X = np.zeros((q, q), dtype=complex)
        for a in range(q):
            X[(a + 1) % q, a] = 1.0
        self.add_op('X', X, hc='Xhc')
        self.add_op('Xhc', X.conj().T, hc='X')
        self.add_op('Xphc', X + X.conj().T, hc='Xphc')


class PottsChain(CouplingMPOModel, NearestNeighborModel):
    def init_sites(self, model_params):
        return PottsSite(model_params.get('q', 3))
    def init_terms(self, model_params):
        J = model_params.get('J', 1.0)
        g = model_params.get('g', 1.0)
        q = model_params.get('q', 3)
        for a in range(q):
            self.add_coupling(-J, 0, f'P{a}', 0, f'P{a}', 1)
        self.add_onsite(-g, 0, 'Xphc')


def extract_xi(psi, n, op='Sigmax'):
    """Extract ξ from bulk connected correlator decay."""
    i_start = n // 4
    i_end = 3 * n // 4
    sites2 = list(range(i_start, i_end))

    C = psi.correlation_function(op, op, sites1=[i_start], sites2=sites2)
    C = np.array(C).flatten()

    exp_vals = psi.expectation_value(op)
    exp_O = np.array([exp_vals[i] for i in sites2])
    C_conn = C - exp_O[0] * exp_O

    r_vals = np.arange(1, len(C_conn))
    C_abs = np.abs(C_conn[1:]).flatten()

    mask = C_abs > 1e-14
    if np.sum(mask) < 3:
        return float('inf'), 0.0

    r_fit = r_vals[mask]
    C_fit = C_abs[mask]
    log_C = np.log(C_fit)

    coeffs = np.polyfit(r_fit, log_C, 1)
    xi = -1.0 / coeffs[0] if coeffs[0] < 0 else float('inf')

    residuals = log_C - np.polyval(coeffs, r_fit)
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((log_C - np.mean(log_C))**2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0

    return xi, r2


def run_point(q, n, g, chi_max=40):
    """Run DMRG and extract ξ for one (q, n, g) point."""
    t0 = time.time()
    if q == 2:
        model = TFIChain({'L': n, 'J': 1.0, 'g': g, 'bc_MPS': 'finite'})
        op = 'Sigmax'
    else:
        model = PottsChain({'L': n, 'q': q, 'J': 1.0, 'g': g, 'bc_MPS': 'finite'})
        op = 'P0'

    psi = MPS.from_product_state(model.lat.mps_sites(), [0]*n, bc='finite')
    eng = dmrg.TwoSiteDMRGEngine(psi, model, {
        'mixer': True, 'max_E_err': 1e-8,
        'trunc_params': {'chi_max': chi_max, 'svd_min': 1e-10},
        'max_sweeps': 20,
    })
    E0, _ = eng.run()
    t_dmrg = time.time() - t0

    xi, r2 = extract_xi(psi, n, op=op)
    S_half = float(psi.entanglement_entropy()[n//2 - 1])

    return float(E0), xi, r2, S_half, t_dmrg


def fss_collapse_quality(nu, g_c, g_arr, xi_arr, n_arr):
    """Compute FSS collapse quality for ξ/n = f((g-g_c)·n^{1/ν}).

    Returns sum of squared differences between interpolated curves.
    """
    # Group by n
    sizes = sorted(set(n_arr))
    if len(sizes) < 2:
        return 1e10

    # For each size, compute scaled x = (g-g_c) * n^(1/nu) and y = xi/n
    curves = {}
    for s in sizes:
        mask = (n_arr == s)
        x = (g_arr[mask] - g_c) * s**(1.0/nu)
        y = xi_arr[mask] / s
        order = np.argsort(x)
        curves[s] = (x[order], y[order])

    # Compare each pair of curves in their overlapping x-range
    total_err = 0.0
    count = 0
    for i, s1 in enumerate(sizes):
        for s2 in sizes[i+1:]:
            x1, y1 = curves[s1]
            x2, y2 = curves[s2]

            # Overlapping range
            x_min = max(x1.min(), x2.min())
            x_max = min(x1.max(), x2.max())
            if x_min >= x_max:
                continue

            # Interpolate both curves onto common x grid
            x_grid = np.linspace(x_min, x_max, 50)
            y1_interp = np.interp(x_grid, x1, y1)
            y2_interp = np.interp(x_grid, x2, y2)

            # Normalized difference
            diff = (y1_interp - y2_interp)**2
            norm = (y1_interp**2 + y2_interp**2) / 2
            total_err += np.sum(diff / (norm + 1e-15))
            count += len(x_grid)

    return total_err / max(count, 1)


# ================================================================
# Configuration: q values with their g_c estimates and g sweep ranges
# ================================================================
configs = {
    2: {'g_c': 1.0, 'sizes': [20, 40, 80], 'chi': 60,
        'g_range': np.concatenate([
            np.arange(0.80, 1.01, 0.04), np.arange(1.02, 1.21, 0.04),
            np.array([1.30, 1.50])
        ])},
    3: {'g_c': 1.0, 'sizes': [20, 40, 80], 'chi': 40,
        'g_range': np.concatenate([
            np.arange(0.80, 1.01, 0.04), np.arange(1.02, 1.21, 0.04),
            np.array([1.30, 1.50])
        ])},
    4: {'g_c': 0.89, 'sizes': [16, 24, 32], 'chi': 40,
        'g_range': np.concatenate([
            np.arange(0.60, 0.91, 0.05), np.arange(0.92, 1.11, 0.04),
            np.array([1.20, 1.40])
        ])},
    5: {'g_c': 0.45, 'sizes': [12, 16, 24], 'chi': 40,
        'g_range': np.concatenate([
            np.arange(0.20, 0.46, 0.05), np.arange(0.48, 0.71, 0.04),
            np.array([0.80, 1.00])
        ])},
}

# === Timing test ===
print("=== Timing test: q=3, n=40, g=1.5, chi=40 ===")
t0 = time.time()
E0, xi, r2, S, dt_dmrg = run_point(3, 40, 1.5, chi_max=40)
print(f"  E0={E0:.6f}, xi={xi:.3f}, R²={r2:.3f}, S={S:.4f}, time={time.time()-t0:.1f}s")

print("\n=== Timing test: q=5, n=24, g=0.80, chi=40 ===")
t0 = time.time()
E0, xi, r2, S, dt_dmrg = run_point(5, 24, 0.80, chi_max=40)
print(f"  E0={E0:.6f}, xi={xi:.3f}, R²={r2:.3f}, S={S:.4f}, time={time.time()-t0:.1f}s")

# === Main sweeps ===
all_results = {}
for q, cfg in configs.items():
    print(f"\n{'='*60}")
    print(f"q={q} Potts: g_c≈{cfg['g_c']}, sizes={cfg['sizes']}, chi={cfg['chi']}")
    print(f"{'='*60}")

    q_results = {}
    for n in cfg['sizes']:
        print(f"\n  --- n={n} ---")
        n_results = []
        t_start = time.time()

        for g in cfg['g_range']:
            if time.time() - t_start > 300:
                print(f"  Time limit for n={n}, stopping")
                break
            t0 = time.time()
            try:
                E0, xi, r2, S, dt_dmrg = run_point(q, n, float(g), chi_max=cfg['chi'])
                dt = time.time() - t0
                xi_str = f"{xi:.3f}" if xi < 1e6 else f"{xi:.1e}"
                print(f"    g={g:.2f}: xi={xi_str}, R²={r2:.3f}, S={S:.4f}, t={dt:.1f}s")
                n_results.append({
                    'g': float(g), 'E0': E0, 'xi': min(xi, 1e10),
                    'R2': r2, 'S_half': S, 'time': dt, 'status': 'ok'
                })
            except Exception as e:
                print(f"    g={g:.2f}: FAILED ({str(e)[:60]})")
                n_results.append({'g': float(g), 'status': 'failed', 'error': str(e)[:200]})

        q_results[n] = n_results

        # Save incrementally
        all_results[q] = q_results
        save_data = {}
        for qq, qr in all_results.items():
            save_data[str(qq)] = {str(nn): nr for nn, nr in qr.items()}
        with open('results/sprint_049b_potts_xi.json', 'w') as f:
            json.dump({'sprint': '049b', 'method': 'correlation_function_decay',
                       'configs': {str(k): {'g_c': v['g_c'], 'sizes': v['sizes'], 'chi': v['chi']}
                                   for k, v in configs.items()},
                       'data': save_data}, f, indent=2)

    # === FSS collapse for this q ===
    print(f"\n  === FSS collapse for q={q} ===")
    g_all, xi_all, n_all = [], [], []
    for n, ndata in q_results.items():
        for d in ndata:
            if d['status'] == 'ok' and d['R2'] > 0.5 and d['xi'] < n and d['xi'] > 0.1:
                g_all.append(d['g'])
                xi_all.append(d['xi'])
                n_all.append(n)

    if len(g_all) < 6:
        print(f"  Insufficient data ({len(g_all)} points)")
        continue

    g_all = np.array(g_all)
    xi_all = np.array(xi_all)
    n_all = np.array(n_all)
    g_c = cfg['g_c']

    # Scan ν from 0.3 to 4.0
    nu_scan = np.linspace(0.3, 4.0, 75)
    qualities = []
    for nu in nu_scan:
        q_val = fss_collapse_quality(nu, g_c, g_all, xi_all, n_all)
        qualities.append(q_val)
    qualities = np.array(qualities)
    best_idx = np.argmin(qualities)
    nu_best = nu_scan[best_idx]
    q_best = qualities[best_idx]

    # Also check BKT (very large ν)
    q_bkt = fss_collapse_quality(10.0, g_c, g_all, xi_all, n_all)

    print(f"  Best ν = {nu_best:.2f} (quality = {q_best:.4f})")
    print(f"  ν=1.0 quality: {fss_collapse_quality(1.0, g_c, g_all, xi_all, n_all):.4f}")
    if q >= 3:
        print(f"  ν=5/6 quality: {fss_collapse_quality(5/6, g_c, g_all, xi_all, n_all):.4f}")
    if q >= 4:
        print(f"  ν=2.0 quality: {fss_collapse_quality(2.0, g_c, g_all, xi_all, n_all):.4f}")
    print(f"  ν=10 (BKT proxy) quality: {q_bkt:.4f}")

    # ξ ratio between sizes at each g (diagnostic)
    sizes = sorted(set(n_all))
    if len(sizes) >= 2:
        print(f"  ξ ratios ({sizes[-1]}/{sizes[0]}):")
        for g in sorted(set(g_all)):
            xi_small = [xi_all[i] for i in range(len(g_all))
                        if abs(g_all[i]-g)<0.001 and n_all[i]==sizes[0]]
            xi_large = [xi_all[i] for i in range(len(g_all))
                        if abs(g_all[i]-g)<0.001 and n_all[i]==sizes[-1]]
            if xi_small and xi_large:
                ratio = xi_large[0] / xi_small[0]
                n_ratio = sizes[-1] / sizes[0]
                print(f"    g={g:.2f}: ξ ratio={ratio:.2f} (n ratio={n_ratio:.1f})")

# Final save
save_data = {}
for qq, qr in all_results.items():
    save_data[str(qq)] = {str(nn): nr for nn, nr in qr.items()}
with open('results/sprint_049b_potts_xi.json', 'w') as f:
    json.dump({'sprint': '049b', 'method': 'correlation_function_decay',
               'configs': {str(k): {'g_c': v['g_c'], 'sizes': v['sizes'], 'chi': v['chi']}
                           for k, v in configs.items()},
               'data': save_data}, f, indent=2)

print("\nResults saved to results/sprint_049b_potts_xi.json")
