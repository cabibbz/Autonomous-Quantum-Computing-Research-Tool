"""Sprint 124b: iDMRG chi_F with warm-start — fix convergence issues from 124a.

Key fix: start the g+dg run from the converged psi(g), NOT from random state.
This ensures both states are in the same "basin" and the overlap is meaningful.

Also extract xi from correlation_length() with proper convergence checks.
"""
import numpy as np
import json, time, os, sys, copy
from scipy.optimize import curve_fit
from db_utils import record

results = {
    'experiment': '124b_idmrg_chif_warmstart',
    'sprint': 124,
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'q2': {}, 'q4': {}, 'q5': {},
}

def save():
    outpath = os.path.join(os.path.dirname(__file__), '..', 'results',
                           'sprint_124b_idmrg_chif_warmstart.json')
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    with open(outpath, 'w') as f:
        json.dump(results, f, indent=2, default=str)

# === TeNPy setup ===
from tenpy.models.model import CouplingMPOModel
from tenpy.networks.site import Site
from tenpy.linalg import np_conserved as npc
from tenpy.networks.mps import MPS
from tenpy.algorithms import dmrg

import warnings
warnings.filterwarnings('ignore')

class SqPottsSite(Site):
    def __init__(self, q_val):
        leg = npc.LegCharge.from_trivial(q_val)
        Site.__init__(self, leg, [str(a) for a in range(q_val)], sort_charge=False)
        for a in range(q_val):
            P = np.zeros((q_val, q_val), dtype=complex); P[a, a] = 1.0
            self.add_op(f'P{a}', P)
        Sq_field = np.ones((q_val, q_val), dtype=complex) - np.eye(q_val, dtype=complex)
        self.add_op('SqField', Sq_field, hc='SqField')

class SqPottsInfinite(CouplingMPOModel):
    def init_sites(self, model_params):
        return SqPottsSite(model_params.get('q', 4))
    def init_terms(self, model_params):
        J = model_params.get('J', 1.0)
        g = model_params.get('g', 0.25)
        q_val = model_params.get('q', 4)
        for a in range(q_val):
            self.add_coupling(-J, 0, f'P{a}', 0, f'P{a}', 1)
        self.add_onsite(-g, 0, 'SqField')

def get_correlation_length(psi):
    """Extract correlation length from transfer matrix."""
    try:
        xi = psi.correlation_length()
        if xi > 0 and np.isfinite(xi):
            return float(xi)
    except Exception:
        pass
    # Fallback: entropy-based estimate
    S = float(psi.entanglement_entropy()[0])
    return float(np.exp(6.0 * S))  # assume c~1

def idmrg_chif_warmstart(q_val, g_c, dg, chi_max, max_sweeps=80):
    """Compute chi_F density via warm-started iDMRG.

    Key: run g+dg from the converged psi(g) to ensure same basin.
    """
    t0 = time.time()

    # Step 1: converge ground state at g_c
    model0 = SqPottsInfinite({
        'L': 2, 'q': q_val, 'J': 1.0, 'g': g_c,
        'bc_MPS': 'infinite', 'bc_x': 'periodic',
    })
    np.random.seed(42 + q_val + chi_max)
    init = [np.random.randint(q_val) for _ in range(2)]
    psi0 = MPS.from_product_state(model0.lat.mps_sites(), init, bc='infinite')
    eng0 = dmrg.TwoSiteDMRGEngine(psi0, model0, {
        'mixer': True, 'max_E_err': 1e-14,
        'trunc_params': {'chi_max': chi_max, 'svd_min': 1e-15},
        'max_sweeps': max_sweeps,
    })
    E0, _ = eng0.run()
    eps0 = float(E0) / 2.0
    S0 = float(psi0.entanglement_entropy()[0])
    xi0 = get_correlation_length(psi0)
    chi_used = max(psi0.chi)

    t_mid = time.time()

    # Step 2: warm-start g+dg from converged psi0
    model1 = SqPottsInfinite({
        'L': 2, 'q': q_val, 'J': 1.0, 'g': g_c + dg,
        'bc_MPS': 'infinite', 'bc_x': 'periodic',
    })
    # Deep copy the converged state as initial guess for g+dg
    psi1 = psi0.copy()
    eng1 = dmrg.TwoSiteDMRGEngine(psi1, model1, {
        'mixer': True, 'max_E_err': 1e-14,
        'trunc_params': {'chi_max': chi_max, 'svd_min': 1e-15},
        'max_sweeps': max_sweeps,
    })
    E1, _ = eng1.run()
    eps1 = float(E1) / 2.0

    dt = time.time() - t0

    # Step 3: compute overlap
    overlap_raw = psi0.overlap(psi1)
    overlap_abs = abs(overlap_raw)

    # Per-site fidelity: |<psi0|psi1>|^{1/L_uc} where L_uc=2
    fidelity_per_site = overlap_abs ** 0.5  # sqrt for 2-site unit cell
    F2 = fidelity_per_site ** 2

    # chi_F density from fidelity
    if F2 > 1.0 - 1e-12:
        chi_F = -np.log(max(F2, 1e-300)) / dg**2
    else:
        chi_F = (1.0 - F2) / dg**2

    # Also compute energy-based susceptibility (d²E/dg²)
    # Need E at g-dg too for central difference
    model_m = SqPottsInfinite({
        'L': 2, 'q': q_val, 'J': 1.0, 'g': g_c - dg,
        'bc_MPS': 'infinite', 'bc_x': 'periodic',
    })
    psi_m = psi0.copy()
    eng_m = dmrg.TwoSiteDMRGEngine(psi_m, model_m, {
        'mixer': True, 'max_E_err': 1e-14,
        'trunc_params': {'chi_max': chi_max, 'svd_min': 1e-15},
        'max_sweeps': max_sweeps,
    })
    E_m, _ = eng_m.run()
    eps_m = float(E_m) / 2.0

    d2E_dg2 = (eps1 + eps_m - 2 * eps0) / dg**2

    dt_total = time.time() - t0

    return {
        'chi_F': float(chi_F),
        'd2E_dg2': float(d2E_dg2),
        'fidelity': float(fidelity_per_site),
        'overlap_abs': float(overlap_abs),
        'xi': float(xi0),
        'eps0': float(eps0), 'eps1': float(eps1), 'eps_m': float(eps_m),
        'S': float(S0),
        'chi_max': chi_max,
        'chi_used': int(chi_used),
        'time_s': round(dt_total, 1),
        'time_g_s': round(t_mid - t0, 1),
    }

# ============================================================
print("=" * 70)
print("Sprint 124b: iDMRG chi_F with warm-start")
print("=" * 70)

dg = 1e-3

# === Timing test ===
print("\n--- Timing test: q=4, chi_max=20 ---")
t0 = time.time()
test = idmrg_chif_warmstart(4, 0.25, dg, 20, max_sweeps=40)
dt = time.time() - t0
print(f"  chi_F={test['chi_F']:.4f}, d2E={test['d2E_dg2']:.4f}, "
      f"xi={test['xi']:.1f}, F={test['fidelity']:.8f}, t={dt:.1f}s")

# Estimate schedule based on timing
if dt > 40:
    chi_schedule = [20, 40, 60, 80]
    print(f"  Slow ({dt:.0f}s) — using short schedule")
elif dt > 20:
    chi_schedule = [20, 40, 60, 80, 100]
    print(f"  Moderate ({dt:.0f}s) — using medium schedule")
else:
    chi_schedule = [20, 40, 60, 80, 100, 150]
    print(f"  Fast ({dt:.0f}s) — using full schedule")

total_time = 0

# === q=2 cross-check ===
q = 2; g_c = 0.5
print(f"\n{'='*60}")
print(f"q={q}, g_c={g_c} — cross-check (exact alpha_N=1.0)")
print(f"{'='*60}")

for chi_max in chi_schedule:
    if total_time > 350:
        print(f"  Total time {total_time:.0f}s, stopping to leave room")
        break
    r = idmrg_chif_warmstart(q, g_c, dg, chi_max)
    total_time += r['time_s']
    print(f"  chi={chi_max:3d}: chi_F={r['chi_F']:10.4f}, d2E={r['d2E_dg2']:8.4f}, "
          f"xi={r['xi']:8.1f}, F={r['fidelity']:.8f}, t={r['time_s']:.0f}s")
    results['q2'][str(chi_max)] = r
    save()

# === q=4 main target ===
q = 4; g_c = 0.25
print(f"\n{'='*60}")
print(f"q={q}, g_c={g_c} — main target (periodic alpha=1.77)")
print(f"{'='*60}")

for chi_max in chi_schedule:
    if total_time > 350:
        print(f"  Total time {total_time:.0f}s, stopping")
        break
    r = idmrg_chif_warmstart(q, g_c, dg, chi_max)
    total_time += r['time_s']
    print(f"  chi={chi_max:3d}: chi_F={r['chi_F']:10.4f}, d2E={r['d2E_dg2']:8.4f}, "
          f"xi={r['xi']:8.1f}, F={r['fidelity']:.8f}, t={r['time_s']:.0f}s")
    results['q4'][str(chi_max)] = r
    save()

# === q=5 if time permits ===
q = 5; g_c = 0.2
if total_time < 300:
    print(f"\n{'='*60}")
    print(f"q={q}, g_c={g_c} — walking reference")
    print(f"{'='*60}")

    for chi_max in [c for c in chi_schedule if c <= 80]:
        if total_time > 450:
            break
        r = idmrg_chif_warmstart(q, g_c, dg, chi_max)
        total_time += r['time_s']
        print(f"  chi={chi_max:3d}: chi_F={r['chi_F']:10.4f}, d2E={r['d2E_dg2']:8.4f}, "
              f"xi={r['xi']:8.1f}, F={r['fidelity']:.8f}, t={r['time_s']:.0f}s")
        results['q5'][str(chi_max)] = r
        save()

# === Analysis ===
print(f"\n{'='*70}")
print("SCALING ANALYSIS")
print("=" * 70)

for qlabel, qkey in [("q=2", "q2"), ("q=4", "q4"), ("q=5", "q5")]:
    data = results[qkey]
    if len(data) < 2:
        print(f"\n  {qlabel}: insufficient data ({len(data)} points)")
        continue

    chi_vals = sorted([int(k) for k in data.keys()])
    xis = np.array([data[str(c)]['xi'] for c in chi_vals])
    chifs = np.array([data[str(c)]['chi_F'] for c in chi_vals])
    d2Es = np.array([data[str(c)]['d2E_dg2'] for c in chi_vals])

    mask = (xis > 0) & (chifs > 0) & np.isfinite(xis) & np.isfinite(chifs)
    if mask.sum() < 2:
        print(f"\n  {qlabel}: insufficient valid data")
        continue

    xis_f = xis[mask]; chifs_f = chifs[mask]; d2Es_f = d2Es[mask]
    chi_f = [chi_vals[i] for i in range(len(chi_vals)) if mask[i]]

    log_xi = np.log(xis_f)
    log_chif = np.log(chifs_f)

    print(f"\n  {qlabel}: {len(xis_f)} points")
    print(f"    xi range: [{xis_f.min():.1f}, {xis_f.max():.1f}]")
    print(f"    chi_F range: [{chifs_f.min():.4f}, {chifs_f.max():.4f}]")
    print(f"    d2E range: [{d2Es_f.min():.4f}, {d2Es_f.max():.4f}]")

    # Check monotonicity
    mono = all(chifs_f[i] <= chifs_f[i+1] for i in range(len(chifs_f)-1))
    print(f"    chi_F monotonic: {mono}")

    if len(xis_f) >= 2:
        for i in range(len(xis_f) - 1):
            a = (log_chif[i+1] - log_chif[i]) / (log_xi[i+1] - log_xi[i])
            print(f"    xi=({xis_f[i]:.1f},{xis_f[i+1]:.1f}): alpha = {a:.4f}")

        if len(xis_f) >= 3:
            p = np.polyfit(log_xi, log_chif, 1)
            print(f"    Global alpha (chi_F vs xi): {p[0]:.4f}")

    # d2E/dg2 scaling
    if len(xis_f) >= 2:
        log_d2E = np.log(np.abs(d2Es_f))
        for i in range(len(xis_f) - 1):
            a = (log_d2E[i+1] - log_d2E[i]) / (log_xi[i+1] - log_xi[i])
            print(f"    xi=({xis_f[i]:.1f},{xis_f[i+1]:.1f}): d2E alpha = {a:.4f}")

    results[f'{qkey}_analysis'] = {
        'chi_vals': chi_f,
        'xis': [float(x) for x in xis_f],
        'chi_Fs': [float(c) for c in chifs_f],
        'd2Es': [float(d) for d in d2Es_f],
    }

save()
print(f"\nTotal wall time: {total_time:.0f}s")
print("Results saved.")
