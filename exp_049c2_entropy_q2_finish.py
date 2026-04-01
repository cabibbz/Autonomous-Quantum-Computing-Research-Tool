#!/usr/bin/env python3
"""Sprint 049c2: Complete q=2 TFIM entropy FSS (n=48 remainder + n=64)."""
import numpy as np, json, time, warnings
warnings.filterwarnings('ignore')
from tenpy.models.tf_ising import TFIChain
from tenpy.networks.mps import MPS
from tenpy.algorithms import dmrg

def run_entropy(n, g, chi_max=40):
    t0 = time.time()
    model = TFIChain({'L': n, 'J': 1.0, 'g': g, 'bc_MPS': 'finite'})
    psi = MPS.from_product_state(model.lat.mps_sites(), [i % 2 for i in range(n)], bc='finite')
    eng = dmrg.TwoSiteDMRGEngine(psi, model, {
        'mixer': True, 'max_E_err': 1e-10,
        'trunc_params': {'chi_max': chi_max, 'svd_min': 1e-12},
        'max_sweeps': 30,
    })
    E0, _ = eng.run()
    S = float(psi.entanglement_entropy()[n//2 - 1])
    return float(E0), S, time.time() - t0

with open('results/sprint_049b_entropy_fss.json') as f:
    results = json.load(f)

def save():
    with open('results/sprint_049b_entropy_fss.json', 'w') as f:
        json.dump(results, f, indent=2)

g_vals = [0.70, 0.80, 0.85, 0.90, 0.93, 0.95, 0.97, 0.99,
          1.00, 1.01, 1.03, 1.05, 1.07, 1.10, 1.15, 1.20, 1.30, 1.50]

# n=48: continue from where we stopped (already have g<=1.05)
existing_48 = results.get('q2', {}).get('48', [])
done_gs = {d['g'] for d in existing_48}
remaining_48 = [g for g in g_vals if g not in done_gs]

if remaining_48:
    print(f"=== q=2 n=48: {len(remaining_48)} points remaining ===", flush=True)
    for g in remaining_48:
        E0, S, dt = run_entropy(48, g)
        print(f"  g={g:.2f}: S={S:.4f}, t={dt:.1f}s", flush=True)
        existing_48.append({'g': g, 'E0': E0, 'S_half': S, 'time': dt})
    results['q2']['48'] = existing_48
    save()

# n=64
if '64' not in results.get('q2', {}):
    print(f"\n=== q=2 n=64 ===", flush=True)
    n64_data = []
    for g in g_vals:
        E0, S, dt = run_entropy(64, g, chi_max=60)
        print(f"  g={g:.2f}: S={S:.4f}, t={dt:.1f}s", flush=True)
        n64_data.append({'g': g, 'E0': E0, 'S_half': S, 'time': dt})
    results['q2']['64'] = n64_data
    save()

print("\nq=2 entropy FSS complete!", flush=True)
