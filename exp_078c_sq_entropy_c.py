#!/usr/bin/env python3
"""Sprint 078c: Central charge c for S_q Potts from DMRG entropy profile.

Use entropy profile S(x) at g_c = 1/q for q=5,3 to extract c.
For open BC chain length L: S(x) = (c/6)·ln[(2L/π)·sin(πx/L)] + const (Calabrese-Cardy).
DMRG entropy profile gives S at each bond cut.

Compare S_q Potts c to:
- Hybrid model c (from Sprint 056: c ≈ 0.40·ln(q-1) + 0.55)
- Complex CFT predictions (Gorbenko et al.): Re(c) ≈ 1.138 for q=5
- Ma & He (PRB 2019): c_eff from Hermitian S_q chain

Also verify q=3 exact result: S_q ≡ hybrid at q=3, c should be same.
"""
import numpy as np
import json, time

from tenpy.models.model import CouplingMPOModel, NearestNeighborModel
from tenpy.networks.site import Site
from tenpy.linalg import np_conserved as npc
from tenpy.networks.mps import MPS
from tenpy.algorithms import dmrg

results = {
    'experiment': '078c_sq_entropy_c',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
}

def save():
    with open("results/sprint_078c_sq_entropy_c.json", "w") as f:
        json.dump(results, f, indent=2, default=str)


class SqPottsSite(Site):
    def __init__(self, q):
        leg = npc.LegCharge.from_trivial(q)
        Site.__init__(self, leg, [str(a) for a in range(q)], sort_charge=False)
        for a in range(q):
            P = np.zeros((q, q), dtype=complex); P[a, a] = 1.0
            self.add_op(f'P{a}', P)
        Sq_field = np.ones((q, q), dtype=complex) - np.eye(q, dtype=complex)
        self.add_op('SqField', Sq_field, hc='SqField')


class SqPottsChain(CouplingMPOModel, NearestNeighborModel):
    def init_sites(self, model_params):
        return SqPottsSite(model_params.get('q', 5))
    def init_terms(self, model_params):
        J = model_params.get('J', 1.0)
        g = model_params.get('g', 0.2)
        q_val = model_params.get('q', 5)
        for a in range(q_val):
            self.add_coupling(-J, 0, f'P{a}', 0, f'P{a}', 1)
        self.add_onsite(-g, 0, 'SqField')


def dmrg_entropy_profile(n, q, g, chi_max=80):
    """Get ground state entropy profile via DMRG."""
    model = SqPottsChain({'L': n, 'q': q, 'J': 1.0, 'g': g, 'bc_MPS': 'finite'})
    np.random.seed(42 + n + q * 100)
    init = [np.random.randint(q) for _ in range(n)]
    psi = MPS.from_product_state(model.lat.mps_sites(), init, bc='finite')
    eng = dmrg.TwoSiteDMRGEngine(psi, model, {
        'mixer': True, 'max_E_err': 1e-12,
        'trunc_params': {'chi_max': chi_max, 'svd_min': 1e-14},
        'max_sweeps': 40,
    })
    E0, _ = eng.run()
    S = [float(psi.entanglement_entropy()[i]) for i in range(n - 1)]
    return float(E0), S


def fit_cc(L, S_profile):
    """Fit Calabrese-Cardy formula to extract c.
    S(x) = (c/6)·ln[(2L/π)·sin(πx/L)] + s0
    where x = bond position (1, 2, ..., L-1).
    """
    x_vals = np.arange(1, L)
    chord = np.log((2 * L / np.pi) * np.sin(np.pi * x_vals / L))
    S_arr = np.array(S_profile)

    # Linear fit: S = (c/6)·chord + s0
    A = np.vstack([chord, np.ones_like(chord)]).T
    (slope, s0), residuals, _, _ = np.linalg.lstsq(A, S_arr, rcond=None)
    c = 6 * slope

    # R² goodness of fit
    S_pred = slope * chord + s0
    ss_res = np.sum((S_arr - S_pred) ** 2)
    ss_tot = np.sum((S_arr - np.mean(S_arr)) ** 2)
    R2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0

    return float(c), float(s0), float(R2)


# ---- MAIN ----
print("=" * 60, flush=True)
print("Sprint 078c: Central charge c for S_q Potts", flush=True)
print("=" * 60, flush=True)

# q=3 first (quick, verification: should match hybrid c ≈ 0.8)
print("\n--- q=3, g_c = 1/3 (S_q ≡ hybrid at q=3) ---", flush=True)
q3_data = {}
for n in [12, 16, 20, 24]:
    gc = 1.0 / 3
    chi = min(60 + n, 100)
    t0 = time.time()
    E0, S = dmrg_entropy_profile(n, 3, gc, chi_max=chi)
    dt = time.time() - t0
    c, s0, R2 = fit_cc(n, S)
    print(f"  n={n}, chi={chi}: c={c:.4f}, R²={R2:.6f}, S_mid={S[n//2-1]:.4f} [{dt:.1f}s]", flush=True)
    q3_data[f'n{n}'] = {'n': n, 'c': c, 's0': s0, 'R2': R2, 'S_mid': S[n//2-1],
                         'S_profile': S, 'E0': E0, 'time_s': dt}
    results['q3'] = q3_data
    save()
    if dt > 200:
        print(f"  Time limit, stopping q=3", flush=True)
        break

# q=5: the interesting case
print("\n--- q=5, g_c = 1/5 (complex CFT prediction: Re(c) ≈ 1.138) ---", flush=True)
q5_data = {}
for n in [12, 16, 20, 24]:
    gc = 0.2
    chi = min(60 + n, 100)
    t0 = time.time()
    E0, S = dmrg_entropy_profile(n, 5, gc, chi_max=chi)
    dt = time.time() - t0
    c, s0, R2 = fit_cc(n, S)
    print(f"  n={n}, chi={chi}: c={c:.4f}, R²={R2:.6f}, S_mid={S[n//2-1]:.4f} [{dt:.1f}s]", flush=True)
    q5_data[f'n{n}'] = {'n': n, 'c': c, 's0': s0, 'R2': R2, 'S_mid': S[n//2-1],
                         'S_profile': S, 'E0': E0, 'time_s': dt}
    results['q5'] = q5_data
    save()
    if dt > 200:
        print(f"  Time limit, stopping q=5", flush=True)
        break

# Summary
print(f"\n{'='*60}", flush=True)
print("SUMMARY: Central charge c", flush=True)
print(f"{'='*60}", flush=True)

print("\nq=3 (exact: c=4/5=0.800 for 3-state Potts):", flush=True)
for k, v in q3_data.items():
    print(f"  {k}: c = {v['c']:.4f} (R² = {v['R2']:.6f})", flush=True)

print("\nq=5 (complex CFT: Re(c) ≈ 1.138):", flush=True)
for k, v in q5_data.items():
    print(f"  {k}: c = {v['c']:.4f} (R² = {v['R2']:.6f})", flush=True)

# Hybrid model c for comparison (from Sprint 056: c ≈ 0.40·ln(q-1)+0.55)
c_hybrid_q5 = 0.40 * np.log(4) + 0.55
print(f"\nComparison:", flush=True)
print(f"  Hybrid c(q=5) ≈ {c_hybrid_q5:.3f} (from Sprint 056 formula)", flush=True)
print(f"  Complex CFT Re(c) for q=5 S_q: ≈ 1.138 (Gorbenko et al.)", flush=True)

save()
print("\nSaved to results/sprint_078c_sq_entropy_c.json", flush=True)

from db_utils import record
for q_val, data in [(3, q3_data), (5, q5_data)]:
    for k, v in data.items():
        record(sprint=78, model='sq_potts', q=q_val, n=v['n'],
               quantity='c_eff', value=v['c'],
               method=f'dmrg_CC_fit', notes=f'gc=1/{q_val}, R2={v["R2"]:.4f}')
print("Recorded to DB.", flush=True)
