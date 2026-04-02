#!/usr/bin/env python3
"""Sprint 071b: q=5 Ly=2 cylinder DMRG — transition scan.

Lean scan with chi=20. Find g_c via entropy peak, compare to q=2.
Key question: does entropy behavior differ qualitatively from q=2?
"""
import numpy as np
import json
import time

import tenpy
from tenpy.models.model import CouplingMPOModel
from tenpy.networks.site import Site
from tenpy.algorithms import dmrg
from tenpy.networks.mps import MPS
import tenpy.linalg.np_conserved as npc
from tenpy.models.lattice import Square


class ZqSite(Site):
    def __init__(self, q):
        ops = {}
        X = np.zeros((q, q), dtype=complex)
        for s in range(q):
            X[(s + 1) % q, s] = 1.0
        ops['XpXd'] = (X + X.conj().T).real.astype(float)
        for a in range(q):
            Pa = np.zeros((q, q))
            Pa[a, a] = 1.0
            ops[f'P{a}'] = Pa
        leg = npc.LegCharge.from_trivial(q)
        super().__init__(leg, ['s' + str(i) for i in range(q)], **ops)


class HybridPottsCylinder(CouplingMPOModel):
    def init_lattice(self, model_params):
        q = model_params.get('q', 5)
        Lx = model_params.get('Lx', 10)
        Ly = model_params.get('Ly', 2)
        site = ZqSite(q)
        lat = Square(Lx, Ly, site, bc=['open', 'periodic'], bc_MPS='finite')
        return lat

    def init_terms(self, model_params):
        q = model_params.get('q', 5)
        g = model_params.get('g', 1.0)
        for a in range(q):
            self.add_coupling(-1.0, 0, f'P{a}', 0, f'P{a}', [1, 0])
            self.add_coupling(-1.0, 0, f'P{a}', 0, f'P{a}', [0, 1])
        self.add_onsite(-g, 0, 'XpXd')


def run_dmrg(q, Lx, g, chi=20):
    M = HybridPottsCylinder({'q': q, 'g': float(g), 'Lx': Lx, 'Ly': 2})
    psi = MPS.from_lat_product_state(M.lat, [[['s0']]])
    eng = dmrg.TwoSiteDMRGEngine(psi, M, {
        'mixer': True, 'max_E_err': 1e-8,
        'trunc_params': {'chi_max': chi, 'svd_min': 1e-10, 'trunc_cut': None},
        'max_sweeps': 20,
        'max_trunc_err': 1.0,  # allow large truncation
    })
    E0, _ = eng.run()
    S = psi.entanglement_entropy()
    n = 2 * Lx
    S_mid = float(S[n // 2 - 1])
    return float(E0), S_mid, [float(s) for s in S]


# ---- MAIN ----
print("=" * 60, flush=True)
print("Sprint 071b: q=5 Ly=2 cylinder DMRG transition scan", flush=True)
print("=" * 60, flush=True)

results = {
    'experiment': '071b_cylinder_q5',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'Ly': 2, 'q': 5, 'chi': 20,
}

# g_c(1D)=0.441, g_c(2D)=1.588. Expect cylinder ~0.7-1.2
# Coarse scan
g_scan = np.linspace(0.3, 1.8, 16)
print(f"\n--- Coarse scan: 16 g-values ---", flush=True)

scan_data = {}
for Lx in [6, 8, 10]:
    print(f"\n  Lx={Lx} (n={2*Lx}):", flush=True)
    data = []
    t_total = time.time()
    for g in g_scan:
        t0 = time.time()
        E0, S_mid, _ = run_dmrg(5, Lx, g)
        dt = time.time() - t0
        data.append({'g': round(float(g), 4), 'E0': E0, 'S_mid': S_mid,
                     'E_per_site': E0 / (2 * Lx), 'time_s': round(dt, 1)})
    dt_total = time.time() - t_total
    peak = max(data, key=lambda d: d['S_mid'])
    print(f"  S_mid peak: {peak['S_mid']:.4f} at g={peak['g']:.3f} (total {dt_total:.0f}s)", flush=True)
    scan_data[Lx] = data

results['coarse_scan'] = {str(k): v for k, v in scan_data.items()}

# Fine scan
gc_peaks = [max(scan_data[Lx], key=lambda d: d['S_mid'])['g'] for Lx in scan_data]
rough_gc = np.mean(gc_peaks)
print(f"\n--- Fine scan around g_c ~ {rough_gc:.3f} ---", flush=True)

g_fine = np.linspace(rough_gc - 0.15, rough_gc + 0.15, 11)
fine_data = {}
for Lx in [6, 8, 10]:
    print(f"\n  Lx={Lx}:", flush=True)
    data = []
    for g in g_fine:
        t0 = time.time()
        E0, S_mid, _ = run_dmrg(5, Lx, g)
        dt = time.time() - t0
        data.append({'g': round(float(g), 4), 'E0': E0, 'S_mid': S_mid, 'time_s': round(dt, 1)})
    peak = max(data, key=lambda d: d['S_mid'])
    print(f"  S_mid peak: {peak['S_mid']:.4f} at g={peak['g']:.3f}", flush=True)
    fine_data[Lx] = data

results['fine_scan'] = {str(k): v for k, v in fine_data.items()}

gc_fine = [max(fine_data[Lx], key=lambda d: d['S_mid'])['g'] for Lx in fine_data]
gc_est = np.mean(gc_fine)
results['gc_estimate'] = gc_est
print(f"\n  g_c(Ly=2 cylinder, q=5) = {gc_est:.4f}", flush=True)

# Entropy scaling at g_c
print(f"\n--- Entropy scaling at g_c={gc_est:.4f} ---", flush=True)
Lx_list = [4, 6, 8, 10, 12]
S_scaling = []
for Lx in Lx_list:
    t0 = time.time()
    E0, S_mid, S_prof = run_dmrg(5, Lx, gc_est)
    dt = time.time() - t0
    print(f"  Lx={Lx}: S_mid={S_mid:.5f} ({dt:.1f}s)", flush=True)
    S_scaling.append({'Lx': Lx, 'S_mid': S_mid, 'E0': E0})

results['entropy_scaling'] = S_scaling

# Fit c
mask = [d['Lx'] >= 8 for d in S_scaling]
lnL_fit = np.log([d['Lx'] for d, m in zip(S_scaling, mask) if m])
S_fit = np.array([d['S_mid'] for d, m in zip(S_scaling, mask) if m])
if len(lnL_fit) >= 2:
    coeffs = np.polyfit(lnL_fit, S_fit, 1)
    c_eff = 6 * coeffs[0]
    results['c_eff'] = c_eff
    print(f"\n  c_eff = {c_eff:.4f} (q=5 1D: c≈1.10)", flush=True)

# Energy per site across transition (latent heat test)
print(f"\n--- Energy per site near g_c ---", flush=True)
Lx = 8
dg = 0.02
g_energy = np.arange(max(0.1, gc_est - 0.25), gc_est + 0.25, dg)
E_per_site = []
for g in g_energy:
    E0, _, _ = run_dmrg(5, Lx, g)
    E_per_site.append(E0 / (2 * Lx))

# Compute dE/dg
dEdg = np.gradient(E_per_site, dg)
d2Edg2 = np.gradient(dEdg, dg)

# Check for discontinuity in dE/dg
gc_idx = np.argmin(np.abs(g_energy - gc_est))
print(f"  dE/dg at g_c: {dEdg[gc_idx]:.6f}", flush=True)
print(f"  max |d²E/dg²|: {np.max(np.abs(d2Edg2)):.4f} at g={g_energy[np.argmax(np.abs(d2Edg2))]:.3f}", flush=True)

# Check dE/dg smoothness: max jump between consecutive dE/dg values
max_jump = np.max(np.abs(np.diff(dEdg)))
print(f"  max dE/dg jump: {max_jump:.6f} (first-order → large jump)", flush=True)

results['energy_deriv'] = {
    'Lx': Lx, 'g_values': g_energy.tolist(),
    'E_per_site': E_per_site,
    'dEdg': dEdg.tolist(),
    'max_dEdg_jump': float(max_jump),
}

# Save
with open("results/sprint_071b_cylinder_q5.json", "w") as f:
    json.dump(results, f, indent=2, default=str)
print("\nSaved to results/sprint_071b_cylinder_q5.json", flush=True)

from db_utils import record
record(sprint=71, model='hybrid_cyl', q=5, n=2,
       quantity='gc', value=gc_est, method='entropy_peak_Ly2',
       notes='Ly=2 cylinder chi=20')
if 'c_eff' in results:
    record(sprint=71, model='hybrid_cyl', q=5, n=2,
           quantity='c', value=results['c_eff'], method='S_vs_lnLx_Ly2',
           notes='Ly=2 cylinder chi=20')

print(f"\n--- Summary ---", flush=True)
print(f"g_c(1D)=0.441, g_c(2D)=1.588", flush=True)
print(f"g_c(Ly=2 cyl, q=5) = {gc_est:.4f}", flush=True)
print(f"Ratio g_c(cyl)/g_c(1D) = {gc_est/0.441:.3f}", flush=True)
print(f"q=2 ratio was {0.4533/0.250:.3f}", flush=True)
