#!/usr/bin/env python3
"""Sprint 055a: iDMRG central charge extraction — TFIM validation.

At the TFIM critical point g=1.0, extract c from S vs ln(xi) at multiple chi.
CFT prediction: c = 1/2.
Uses TeNPy's built-in TFIChain with bc_MPS='infinite'.
"""
import numpy as np, json, time, warnings, sys
warnings.filterwarnings('ignore')

from tenpy.models.tf_ising import TFIChain
from tenpy.networks.mps import MPS
from tenpy.algorithms import dmrg

results = {'experiment': '055a', 'model': 'TFIM', 'q': 2, 'g_c': 1.0,
           'method': 'iDMRG S vs ln(xi)', 'c_exact': 0.5}

chi_values = [10, 20, 40, 80, 160]
data = []

for chi_max in chi_values:
    t0 = time.time()
    model = TFIChain({'L': 2, 'J': 1.0, 'g': 1.0, 'bc_MPS': 'infinite'})
    psi = MPS.from_lat_product_state(model.lat, [['up']])
    eng = dmrg.TwoSiteDMRGEngine(psi, model, {
        'mixer': True, 'max_E_err': 1e-10,
        'trunc_params': {'chi_max': chi_max, 'svd_min': 1e-14},
        'max_sweeps': 100,
        'norm_tol': 1e-3,       # relax norm tolerance for infinite MPS
        'norm_tol_final': 1e-5, # relax final norm tolerance
    })
    E, _ = eng.run()
    S = float(np.mean(psi.entanglement_entropy()))
    xi = float(psi.correlation_length2())
    chi_actual = max(psi.chi)
    dt = time.time() - t0

    print(f"chi={chi_max:4d}: E/site={E:.8f}, S={S:.6f}, xi={xi:.2f}, chi_actual={chi_actual}, time={dt:.1f}s", flush=True)
    data.append({
        'chi_max': chi_max, 'chi_actual': chi_actual,
        'E_per_site': float(E), 'S': S, 'xi': xi, 'time': dt
    })
    sys.stdout.flush()

results['data'] = data

# Extract c from S = (c/6)*ln(xi) + const
Ss = np.array([d['S'] for d in data])
xis = np.array([d['xi'] for d in data])
ln_xi = np.log(xis)

# Full fit
A = np.vstack([ln_xi, np.ones(len(ln_xi))]).T
slope, intercept = np.linalg.lstsq(A, Ss, rcond=None)[0]
c_full = 6 * slope

# Pairwise c
c_pairs = []
for i in range(len(data) - 1):
    c_p = 6 * (Ss[i+1] - Ss[i]) / (ln_xi[i+1] - ln_xi[i])
    c_pairs.append(float(c_p))

print(f"\nFull fit: c = {c_full:.4f} (exact: 0.500)")
print(f"Pairwise c: {' -> '.join(f'{c:.4f}' for c in c_pairs)}")

results['c_full_fit'] = float(c_full)
results['c_pairwise'] = c_pairs

with open('results/exp_055a.json', 'w') as f:
    json.dump(results, f, indent=2)
print("\nSaved results/exp_055a.json")
