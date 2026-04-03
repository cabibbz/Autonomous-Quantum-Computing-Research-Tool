#!/usr/bin/env python3
"""Sprint 049a: Correlation length from correlation function decay — TFIM validation.

Extract ξ by fitting ⟨σx_i σx_j⟩_conn ~ exp(-|i-j|/ξ) from DMRG ground state.
Use bulk sites to avoid boundary effects. For TFIM, ν=1 exactly, g_c=1.0.

Strategy: use TeNPy's correlation_function method at n=40,80 with χ=60.
Fit ξ from exponential decay of connected correlator in the bulk.
Then fit ν from ξ(g) ~ |g - g_c|^(-ν) on the disordered side.
"""
import numpy as np
import numpy.linalg as la
import json, time
from scipy.optimize import curve_fit

from tenpy.models.tf_ising import TFIChain
from tenpy.networks.mps import MPS
from tenpy.algorithms import dmrg
import warnings
warnings.filterwarnings('ignore')


def extract_xi_from_correlations(psi, n, op1='Sigmax', op2='Sigmax'):
    """Extract correlation length from connected correlator decay in the bulk.

    Compute C(r) = ⟨O_i O_{i+r}⟩ - ⟨O_i⟩⟨O_{i+r}⟩ for bulk sites,
    then fit |C(r)| ~ A * exp(-r/ξ).
    """
    # Use middle quarter of chain to avoid boundaries
    i_start = n // 4
    i_end = 3 * n // 4

    # Compute correlation matrix using TeNPy
    sites2 = list(range(i_start, i_end))
    C = psi.correlation_function(op1, op2, sites1=[i_start], sites2=sites2)
    C = np.array(C).flatten()  # Ensure 1D

    # Connected correlator: subtract ⟨O⟩²
    exp_vals = psi.expectation_value(op1)
    exp_O = np.array([exp_vals[i] for i in sites2])
    exp_O_start = exp_O[0]
    C_conn = C - exp_O_start * exp_O

    # Extract r > 0 points
    r_vals = np.arange(1, len(C_conn))
    C_abs = np.abs(C_conn[1:]).flatten()

    # Filter out near-zero values
    mask = C_abs > 1e-14
    if np.sum(mask) < 3:
        return float('inf'), 0.0, []

    r_fit = r_vals[mask]
    C_fit = C_abs[mask]
    log_C = np.log(C_fit)

    # Linear fit: log|C| = -r/ξ + const
    try:
        coeffs = np.polyfit(r_fit, log_C, 1)
        xi = -1.0 / coeffs[0] if coeffs[0] < 0 else float('inf')
        # R² quality
        residuals = log_C - np.polyval(coeffs, r_fit)
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((log_C - np.mean(log_C))**2)
        r_squared = 1 - ss_res / ss_tot if ss_tot > 0 else 0
    except:
        xi = float('inf')
        r_squared = 0.0

    return xi, r_squared, list(zip(r_fit.tolist(), C_fit.tolist()))


def run_tfim(n, g, chi_max=60):
    """Run DMRG for TFIM and extract correlation length from correlator decay."""
    t0 = time.time()
    model = TFIChain({'L': n, 'J': 1.0, 'g': g, 'bc_MPS': 'finite'})
    psi = MPS.from_product_state(model.lat.mps_sites(), [0]*n, bc='finite')
    eng = dmrg.TwoSiteDMRGEngine(psi, model, {
        'mixer': True, 'max_E_err': 1e-10,
        'trunc_params': {'chi_max': chi_max, 'svd_min': 1e-12},
        'max_sweeps': 30,
    })
    E0, _ = eng.run()
    t_dmrg = time.time() - t0

    t1 = time.time()
    xi, r2, corr_data = extract_xi_from_correlations(psi, n)
    t_corr = time.time() - t1

    # Half-chain entropy
    S_half = float(psi.entanglement_entropy()[n//2 - 1])

    return float(E0), xi, r2, S_half, t_dmrg, t_corr


# === Timing test ===
print("=== Timing test: TFIM n=40, g=1.5 (disordered), chi=60 ===")
t0 = time.time()
E0, xi, r2, S, dt_dmrg, dt_corr = run_tfim(40, 1.5, chi_max=60)
dt = time.time() - t0
print(f"  E0={E0:.8f}, xi={xi:.4f}, R²={r2:.4f}, S_half={S:.4f}")
print(f"  time: DMRG={dt_dmrg:.1f}s, corr={dt_corr:.1f}s, total={dt:.1f}s")

# Exact TFIM: ξ_exact = 1/ln(g) for g>1 in thermodynamic limit
xi_exact = 1.0 / np.log(1.5)
print(f"  ξ_exact (thermo limit) = {xi_exact:.4f}, measured = {xi:.4f}")

# === Main sweep ===
# Dense near g_c=1.0, sparser away
g_values = [0.50, 0.70, 0.80, 0.90, 0.95, 0.98, 1.00,
            1.02, 1.05, 1.08, 1.10, 1.15, 1.20, 1.30, 1.50, 2.00]

results = {}
for n in [40, 80]:
    print(f"\n=== TFIM n={n}, chi=60 ===")
    results[n] = []
    t_start = time.time()

    for g in g_values:
        t0 = time.time()
        try:
            E0, xi, r2, S, dt_dmrg, dt_corr = run_tfim(n, g, chi_max=60)
            dt = time.time() - t0
            xi_str = f"{xi:.3f}" if xi < 1e6 else f"{xi:.2e}"
            print(f"  g={g:.2f}: xi={xi_str}, R²={r2:.3f}, S={S:.4f}, time={dt:.1f}s")
            results[n].append({
                'g': g, 'E0': E0, 'xi': xi if xi < 1e15 else 1e15,
                'R2': r2, 'S_half': S, 'time': dt, 'status': 'ok'
            })
        except Exception as e:
            dt = time.time() - t0
            print(f"  g={g:.2f}: FAILED ({e}), time={dt:.1f}s")
            results[n].append({'g': g, 'status': 'failed', 'error': str(e)[:200]})

        # Save incrementally
        with open('results/sprint_049a_tfim_xi.json', 'w') as f:
            json.dump({'sprint': '049a', 'model': 'TFIM', 'chi_max': 60,
                       'method': 'correlation_function_decay',
                       'sizes': {str(k): v for k, v in results.items()},
                       'g_c_exact': 1.0, 'nu_exact': 1.0}, f, indent=2)

# === Analysis: fit ν ===
print("\n=== Analysis: ν extraction from ξ(g) ===")
for n, data in results.items():
    # Disordered side: g > g_c, with good R² and ξ < n (finite-size saturated points excluded)
    ok = [d for d in data if d['status'] == 'ok' and d['g'] > 1.01
          and d['xi'] < n/2 and d['R2'] > 0.5 and d['xi'] > 0.1]
    if len(ok) < 3:
        print(f"  n={n}: insufficient data for fit ({len(ok)} good points)")
        continue

    g_arr = np.array([d['g'] for d in ok])
    xi_arr = np.array([d['xi'] for d in ok])

    # ξ ~ (g - g_c)^(-ν) => log(ξ) = -ν * log(g - g_c) + const
    log_dg = np.log(g_arr - 1.0)
    log_xi = np.log(xi_arr)
    coeffs = np.polyfit(log_dg, log_xi, 1)
    nu_fit = -coeffs[0]

    # Also compare to exact: ξ_exact = 1/ln(g)
    xi_exact = 1.0 / np.log(g_arr)
    ratio = xi_arr / xi_exact

    print(f"  n={n}: ν = {nu_fit:.3f} (exact=1.000), from {len(ok)} points")
    print(f"         g range [{g_arr[0]:.2f}, {g_arr[-1]:.2f}], ξ range [{xi_arr[0]:.2f}, {xi_arr[-1]:.2f}]")
    print(f"         ξ/ξ_exact ratios: {', '.join(f'{r:.3f}' for r in ratio)}")

# Exact comparison table
print("\n=== Exact comparison: ξ_measured vs ξ_exact = 1/ln(g) ===")
for n, data in results.items():
    print(f"  n={n}:")
    for d in data:
        if d['status'] == 'ok' and d['g'] > 1.0:
            xi_exact = 1.0 / np.log(d['g'])
            ratio = d['xi'] / xi_exact if xi_exact > 0 and d['xi'] < 1e10 else float('nan')
            print(f"    g={d['g']:.2f}: ξ={d['xi']:.3f}, ξ_exact={xi_exact:.3f}, ratio={ratio:.3f}, R²={d['R2']:.3f}")

print("\nResults saved to results/sprint_049a_tfim_xi.json")
