#!/usr/bin/env python3
"""Sprint 055b: Central charge from entanglement entropy PROFILE — TFIM validation.

For OBC, S(l) = (c/6)*ln[(2n/pi)*sin(pi*l/n)] + const + oscillations.
Fit c from full S(l) profile at a single large n. Avoids FSS extrapolation.
"""
import numpy as np, json, time, warnings
warnings.filterwarnings('ignore')

from tenpy.models.tf_ising import TFIChain
from tenpy.networks.mps import MPS
from tenpy.algorithms import dmrg


def run_and_extract(n, g, chi_max=80):
    t0 = time.time()
    model = TFIChain({'L': n, 'J': 1.0, 'g': g, 'bc_MPS': 'finite'})
    psi = MPS.from_lat_product_state(model.lat, [['up']])
    eng = dmrg.TwoSiteDMRGEngine(psi, model, {
        'mixer': True, 'max_E_err': 1e-12,
        'trunc_params': {'chi_max': chi_max, 'svd_min': 1e-14},
        'max_sweeps': 50,
    })
    E0, _ = eng.run()
    S_profile = [float(s) for s in psi.entanglement_entropy()]
    chi_actual = max(psi.chi)
    dt = time.time() - t0
    return float(E0), S_profile, chi_actual, dt


def extract_c(S_profile, n):
    ls = np.arange(1, n)
    l_chord = (2 * n / np.pi) * np.sin(np.pi * ls / n)
    ln_chord = np.log(l_chord)
    S = np.array(S_profile)

    # Central half only (avoid edges)
    quarter = n // 4
    mask = (ls >= quarter) & (ls <= 3 * quarter)
    A = np.vstack([ln_chord[mask], np.ones(mask.sum())]).T
    slope, _ = np.linalg.lstsq(A, S[mask], rcond=None)[0]
    c_central = 6 * slope

    # Even bonds in central half (suppress oscillations)
    even_mask = mask & (ls % 2 == 0)
    if even_mask.sum() >= 3:
        A2 = np.vstack([ln_chord[even_mask], np.ones(even_mask.sum())]).T
        slope2, _ = np.linalg.lstsq(A2, S[even_mask], rcond=None)[0]
        c_even = 6 * slope2
    else:
        c_even = c_central

    # Odd bonds
    odd_mask = mask & (ls % 2 == 1)
    if odd_mask.sum() >= 3:
        A3 = np.vstack([ln_chord[odd_mask], np.ones(odd_mask.sum())]).T
        slope3, _ = np.linalg.lstsq(A3, S[odd_mask], rcond=None)[0]
        c_odd = 6 * slope3
    else:
        c_odd = c_central

    return float(c_central), float(c_even), float(c_odd)


# First: timing test at n=32
print("Timing test: n=32, chi=80")
E0, S_prof, chi_act, dt = run_and_extract(32, 1.0, 80)
c_c, c_e, c_o = extract_c(S_prof, 32)
print(f"  n=32: c(central)={c_c:.4f}, c(even)={c_e:.4f}, c(odd)={c_o:.4f}, chi={chi_act}, time={dt:.1f}s")

results = {'experiment': '055b', 'model': 'TFIM', 'g_c': 1.0, 'c_exact': 0.500, 'data': []}
data = []

for n in [32, 48, 64]:
    E0, S_prof, chi_act, dt = run_and_extract(n, 1.0, 80)
    c_c, c_e, c_o = extract_c(S_prof, n)
    print(f"n={n:3d}: c(central)={c_c:.4f}, c(even)={c_e:.4f}, c(odd)={c_o:.4f}, E/site={E0/n:.8f}, chi={chi_act}, time={dt:.1f}s")
    data.append({
        'n': n, 'E0': E0, 'S_half': S_prof[n//2-1],
        'c_central': c_c, 'c_even': c_e, 'c_odd': c_o,
        'chi_actual': chi_act, 'time': dt,
        'S_profile': S_prof,
    })

results['data'] = data

# Also test off-critical as sanity check
print("\nOff-critical g=0.5 (disordered), n=32:")
E0_off, S_off, chi_off, dt_off = run_and_extract(32, 0.5, 40)
c_c_off, _, _ = extract_c(S_off, 32)
print(f"  c(central)={c_c_off:.4f} (should be ~0, not critical)")
results['off_critical'] = {'g': 0.5, 'n': 32, 'c_central': c_c_off}

with open('results/exp_055b.json', 'w') as f:
    json.dump(results, f, indent=2)
print("\nSaved results/exp_055b.json")
