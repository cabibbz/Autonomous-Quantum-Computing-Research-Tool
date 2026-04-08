"""Sprint 124a: iDMRG chi_F scaling for S_q q=4 — probing log corrections in the thermodynamic limit.

Key question: periodic BC exact diag gives alpha=1.77 at n=4-11. Literature predicts
alpha=2 with log corrections: chi_F ~ N^2(ln N)^{-p}. At n<=11, alpha=2 is the
WORST fit (Sprint 122). Can iDMRG reach the asymptotic regime?

Method: iDMRG ground state at g_c and g_c+dg for multiple bond dimensions chi.
- Correlation length xi(chi) from transfer matrix eigenvalues
- Per-site fidelity from infinite MPS overlap
- chi_F density = (1 - fidelity) / dg^2
- Fit chi_F(xi) to extract bulk scaling exponent

Cross-checks: q=2 (exact alpha=1.0) and q=5 (walking, alpha~2.09).
"""
import numpy as np
import json, time, os, sys
from scipy.optimize import curve_fit
from db_utils import record

results = {
    'experiment': '124a_idmrg_chif_scaling',
    'sprint': 124,
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'q2': {}, 'q4': {}, 'q5': {},
}

def save():
    outpath = os.path.join(os.path.dirname(__file__), '..', 'results',
                           'sprint_124a_idmrg_chif_scaling.json')
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

def run_idmrg(q_val, g, chi_max, max_sweeps=60, seed=42):
    """Run iDMRG. Returns (energy_density, psi, chi_used, entropy)."""
    model = SqPottsInfinite({
        'L': 2, 'q': q_val, 'J': 1.0, 'g': g,
        'bc_MPS': 'infinite', 'bc_x': 'periodic',
    })
    np.random.seed(seed + q_val + chi_max + int(g * 10000))
    init = [np.random.randint(q_val) for _ in range(2)]
    psi = MPS.from_product_state(model.lat.mps_sites(), init, bc='infinite')
    eng = dmrg.TwoSiteDMRGEngine(psi, model, {
        'mixer': True, 'max_E_err': 1e-14,
        'trunc_params': {'chi_max': chi_max, 'svd_min': 1e-15},
        'max_sweeps': max_sweeps,
    })
    E0, _ = eng.run()
    eps = float(E0) / 2.0  # energy per site (unit cell = 2 sites)
    chi_used = max(psi.chi)
    S = float(psi.entanglement_entropy()[0])
    return eps, psi, chi_used, S

def get_correlation_length(psi):
    """Extract correlation length from transfer matrix of infinite MPS."""
    try:
        xi = psi.correlation_length()
        return float(xi)
    except Exception:
        # Fallback: use entanglement entropy as proxy
        S = float(psi.entanglement_entropy()[0])
        # For CFT: S = (c/6) ln(xi) => xi = exp(6S/c), approximate c~1
        return float(np.exp(6.0 * S))

def idmrg_chif(q_val, g_c, dg, chi_max, max_sweeps=60):
    """Compute chi_F density from iDMRG at given bond dimension."""
    t0 = time.time()
    eps0, psi0, chi0, S0 = run_idmrg(q_val, g_c, chi_max, max_sweeps)
    eps1, psi1, chi1, S1 = run_idmrg(q_val, g_c + dg, chi_max, max_sweeps, seed=137)
    dt = time.time() - t0

    # Per-site overlap from infinite MPS transfer matrix
    overlap = psi0.overlap(psi1)
    # For infinite MPS, overlap returns the dominant eigenvalue of mixed TM per unit cell
    # Fidelity per site = |overlap|^{1/L_uc} where L_uc=2
    fidelity_per_site = abs(overlap) ** (1.0 / 2.0)  # unit cell = 2

    # chi_F density = (1 - F^2) / dg^2 where F = fidelity per site
    # Actually for small dg: -2 ln(F) / dg^2 is more numerically stable
    F2 = fidelity_per_site ** 2
    if F2 > 0.999999:
        # Use -ln(F^2)/dg^2 ≈ (1-F^2)/dg^2 for numerical stability
        chi_F = -np.log(F2) / dg**2
    else:
        chi_F = (1.0 - F2) / dg**2

    xi = get_correlation_length(psi0)

    return {
        'chi_F': float(chi_F),
        'fidelity': float(fidelity_per_site),
        'overlap_raw': float(abs(overlap)),
        'xi': float(xi),
        'eps0': float(eps0),
        'eps1': float(eps1),
        'S': float(S0),
        'chi_max': chi_max,
        'chi_used': int(chi0),
        'time_s': round(dt, 1),
    }


# ============================================================
print("=" * 70)
print("Sprint 124a: iDMRG chi_F scaling — thermodynamic limit")
print("=" * 70)

dg = 1e-3
chi_schedule = [20, 40, 60, 80, 100, 150, 200]

# Timing test first
print("\n--- Timing test: q=4, chi_max=20 ---")
t0 = time.time()
test = idmrg_chif(4, 0.25, dg, 20, max_sweeps=30)
dt = time.time() - t0
print(f"  chi_max=20: chi_F={test['chi_F']:.4f}, xi={test['xi']:.1f}, t={dt:.1f}s")
if dt > 30:
    print(f"  WARNING: {dt:.0f}s for chi=20 — reducing schedule")
    chi_schedule = [20, 40, 60, 80, 100]

# === q=2 cross-check (exact alpha=1.0) ===
q = 2; g_c = 0.5
print(f"\n{'='*60}")
print(f"q={q}, g_c={g_c} — cross-check (exact alpha=1.0)")
print(f"{'='*60}")

for chi_max in chi_schedule:
    t0 = time.time()
    r = idmrg_chif(q, g_c, dg, chi_max)
    print(f"  chi={chi_max:3d}: chi_F={r['chi_F']:10.4f}, xi={r['xi']:8.1f}, "
          f"S={r['S']:.4f}, chi_used={r['chi_used']}, t={r['time_s']:.0f}s")
    results['q2'][str(chi_max)] = r
    save()

    record(sprint=124, model='sq_potts', q=q, n=chi_max,
           quantity='chi_F_idmrg', value=r['chi_F'], method='idmrg_overlap',
           notes=f'xi={r["xi"]:.1f}')

    if r['time_s'] > 90:
        print(f"  Time limit approaching, stopping q=2")
        break

# === q=4 main target ===
q = 4; g_c = 0.25
print(f"\n{'='*60}")
print(f"q={q}, g_c={g_c} — main target (periodic alpha=1.77, log corrections?)")
print(f"{'='*60}")

for chi_max in chi_schedule:
    t0 = time.time()
    r = idmrg_chif(q, g_c, dg, chi_max)
    print(f"  chi={chi_max:3d}: chi_F={r['chi_F']:10.4f}, xi={r['xi']:8.1f}, "
          f"S={r['S']:.4f}, chi_used={r['chi_used']}, t={r['time_s']:.0f}s")
    results['q4'][str(chi_max)] = r
    save()

    record(sprint=124, model='sq_potts', q=q, n=chi_max,
           quantity='chi_F_idmrg', value=r['chi_F'], method='idmrg_overlap',
           notes=f'xi={r["xi"]:.1f}')

    if r['time_s'] > 90:
        print(f"  Time limit approaching, stopping q=4")
        break

# === q=5 comparison (walking) ===
q = 5; g_c = 0.2
print(f"\n{'='*60}")
print(f"q={q}, g_c={g_c} — walking reference (periodic alpha=2.09)")
print(f"{'='*60}")

# Shorter schedule for q=5 (larger local dim = slower)
chi_q5 = [c for c in chi_schedule if c <= 100]

for chi_max in chi_q5:
    t0 = time.time()
    r = idmrg_chif(q, g_c, dg, chi_max)
    print(f"  chi={chi_max:3d}: chi_F={r['chi_F']:10.4f}, xi={r['xi']:8.1f}, "
          f"S={r['S']:.4f}, chi_used={r['chi_used']}, t={r['time_s']:.0f}s")
    results['q5'][str(chi_max)] = r
    save()

    record(sprint=124, model='sq_potts', q=q, n=chi_max,
           quantity='chi_F_idmrg', value=r['chi_F'], method='idmrg_overlap',
           notes=f'xi={r["xi"]:.1f}')

    if r['time_s'] > 90:
        print(f"  Time limit approaching, stopping q=5")
        break

# === Analysis: chi_F(xi) scaling ===
print(f"\n{'='*70}")
print("SCALING ANALYSIS: chi_F vs xi")
print("=" * 70)

for qlabel, qkey in [("q=2", "q2"), ("q=4", "q4"), ("q=5", "q5")]:
    data = results[qkey]
    if len(data) < 3:
        print(f"\n  {qlabel}: insufficient data ({len(data)} points)")
        continue

    chi_vals = sorted([int(k) for k in data.keys()])
    xis = np.array([data[str(c)]['xi'] for c in chi_vals])
    chifs = np.array([data[str(c)]['chi_F'] for c in chi_vals])

    # Filter out any zero/nan
    mask = (xis > 0) & (chifs > 0) & np.isfinite(xis) & np.isfinite(chifs)
    xis, chifs, chi_vals_f = xis[mask], chifs[mask], [chi_vals[i] for i in range(len(chi_vals)) if mask[i]]

    if len(xis) < 3:
        print(f"\n  {qlabel}: insufficient valid data after filtering")
        continue

    log_xi = np.log(xis)
    log_chif = np.log(chifs)

    print(f"\n  {qlabel}: {len(xis)} points, xi range [{xis.min():.1f}, {xis.max():.1f}]")

    # Pairwise alpha(xi)
    pw_alpha = []
    for i in range(len(xis) - 1):
        a = (log_chif[i+1] - log_chif[i]) / (log_xi[i+1] - log_xi[i])
        pw_alpha.append(a)
        print(f"    xi=({xis[i]:.1f},{xis[i+1]:.1f}): alpha = {a:.4f}")

    # Global power law fit: chi_F = A * xi^alpha
    p = np.polyfit(log_xi, log_chif, 1)
    alpha_global = p[0]
    print(f"    Global alpha (power law): {alpha_global:.4f}")

    # Log-corrected fit: chi_F = A * xi^2 * (ln xi)^{-p}
    if len(xis) >= 4:
        try:
            def log_corrected(log_xi, A, p_corr):
                return A + 2.0 * log_xi + p_corr * np.log(np.maximum(log_xi, 0.1))
            popt, pcov = curve_fit(log_corrected, log_xi, log_chif, p0=[0, -1])
            A_lc, p_lc = popt
            residuals_lc = log_chif - log_corrected(log_xi, *popt)
            ss_res_lc = np.sum(residuals_lc**2)
            ss_tot = np.sum((log_chif - np.mean(log_chif))**2)
            r2_lc = 1 - ss_res_lc / ss_tot if ss_tot > 0 else 0
            print(f"    Log-corrected (alpha=2): p = {p_lc:.4f}, R² = {r2_lc:.6f}")
        except Exception as e:
            print(f"    Log-corrected fit failed: {e}")
            p_lc, r2_lc = None, None
    else:
        p_lc, r2_lc = None, None

    # Power law R²
    residuals_pl = log_chif - np.polyval(p, log_xi)
    ss_res_pl = np.sum(residuals_pl**2)
    ss_tot = np.sum((log_chif - np.mean(log_chif))**2)
    r2_pl = 1 - ss_res_pl / ss_tot if ss_tot > 0 else 0
    print(f"    Power law R² = {r2_pl:.6f}")

    results[f'{qkey}_analysis'] = {
        'chi_vals': chi_vals_f,
        'xis': [float(x) for x in xis],
        'chi_Fs': [float(c) for c in chifs],
        'pw_alpha': [float(a) for a in pw_alpha],
        'alpha_global': float(alpha_global),
        'r2_power_law': float(r2_pl),
        'p_log_corrected': float(p_lc) if p_lc is not None else None,
        'r2_log_corrected': float(r2_lc) if r2_lc is not None else None,
    }

save()
print(f"\nAll results saved.")
