#!/usr/bin/env python3
"""Sprint 071a: Validate Ly=2 cylinder DMRG on q=2 (Ising).

Lean version: chi=20, fewer g points. Find g_c via entropy peak, extract c.
"""
import numpy as np
import json
import time
import sys

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
        q = model_params.get('q', 2)
        Lx = model_params.get('Lx', 10)
        Ly = model_params.get('Ly', 2)
        site = ZqSite(q)
        lat = Square(Lx, Ly, site, bc=['open', 'periodic'], bc_MPS='finite')
        return lat

    def init_terms(self, model_params):
        q = model_params.get('q', 2)
        g = model_params.get('g', 1.0)
        for a in range(q):
            self.add_coupling(-1.0, 0, f'P{a}', 0, f'P{a}', [1, 0])
            self.add_coupling(-1.0, 0, f'P{a}', 0, f'P{a}', [0, 1])
        self.add_onsite(-g, 0, 'XpXd')


def run_dmrg(q, Lx, g, chi=20):
    """Run DMRG on Ly=2 cylinder."""
    M = HybridPottsCylinder({'q': q, 'g': float(g), 'Lx': Lx, 'Ly': 2})
    psi = MPS.from_lat_product_state(M.lat, [[['s0']]])
    eng = dmrg.TwoSiteDMRGEngine(psi, M, {
        'mixer': True, 'max_E_err': 1e-8,
        'trunc_params': {'chi_max': chi, 'svd_min': 1e-10},
        'max_sweeps': 20,
    })
    E0, _ = eng.run()
    S = psi.entanglement_entropy()
    n = 2 * Lx
    S_mid = float(S[n // 2 - 1])
    return float(E0), S_mid, [float(s) for s in S]


# ---- MAIN ----
print("=" * 60, flush=True)
print("Sprint 071a: Ly=2 cylinder DMRG — q=2 validation", flush=True)
print("=" * 60, flush=True)

results = {
    'experiment': '071a_cylinder_q2',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'Ly': 2, 'q': 2, 'chi': 20,
}

# --- Part 1: Coarse entropy scan to find g_c ---
# g_c should be between 1D (0.250) and 2D (0.771)
g_scan = np.linspace(0.2, 0.9, 15)
print(f"\n--- Coarse scan: 15 g-values, Lx=8,12,16 ---", flush=True)

scan_data = {}
for Lx in [8, 12, 16]:
    print(f"\n  Lx={Lx} (n={2*Lx}):", flush=True)
    data = []
    for g in g_scan:
        t0 = time.time()
        E0, S_mid, S_prof = run_dmrg(2, Lx, g)
        dt = time.time() - t0
        data.append({'g': round(float(g), 4), 'E0': E0, 'S_mid': S_mid, 'time_s': round(dt, 1)})
    # Find peak
    peak = max(data, key=lambda d: d['S_mid'])
    print(f"  S_mid peak: {peak['S_mid']:.4f} at g={peak['g']:.3f}", flush=True)
    scan_data[Lx] = data

results['coarse_scan'] = {str(k): v for k, v in scan_data.items()}

# --- Part 2: Fine scan around peak ---
gc_peaks = []
for Lx in scan_data:
    p = max(scan_data[Lx], key=lambda d: d['S_mid'])
    gc_peaks.append(p['g'])
rough_gc = np.mean(gc_peaks)
print(f"\n--- Fine scan around g_c ~ {rough_gc:.3f} ---", flush=True)

g_fine = np.linspace(rough_gc - 0.08, rough_gc + 0.08, 9)
fine_data = {}
for Lx in [12, 16, 20]:
    print(f"\n  Lx={Lx}:", flush=True)
    data = []
    for g in g_fine:
        t0 = time.time()
        E0, S_mid, _ = run_dmrg(2, Lx, g)
        dt = time.time() - t0
        data.append({'g': round(float(g), 4), 'E0': E0, 'S_mid': S_mid, 'time_s': round(dt, 1)})
    peak = max(data, key=lambda d: d['S_mid'])
    print(f"  S_mid peak: {peak['S_mid']:.4f} at g={peak['g']:.3f}", flush=True)
    fine_data[Lx] = data

results['fine_scan'] = {str(k): v for k, v in fine_data.items()}

# Best g_c estimate
gc_fine = []
for Lx in fine_data:
    p = max(fine_data[Lx], key=lambda d: d['S_mid'])
    gc_fine.append(p['g'])
gc_est = np.mean(gc_fine)
results['gc_estimate'] = gc_est
print(f"\n  g_c(Ly=2 cylinder, q=2) = {gc_est:.4f}", flush=True)

# --- Part 3: Entropy scaling at g_c ---
print(f"\n--- Entropy scaling at g_c={gc_est:.4f} ---", flush=True)
Lx_list = [6, 8, 10, 12, 14, 16, 18, 20]
S_scaling = []
for Lx in Lx_list:
    t0 = time.time()
    E0, S_mid, S_prof = run_dmrg(2, Lx, gc_est)
    dt = time.time() - t0
    print(f"  Lx={Lx}: S_mid={S_mid:.5f} ({dt:.1f}s)", flush=True)
    S_scaling.append({'Lx': Lx, 'S_mid': S_mid, 'E0': E0, 'S_profile': S_prof})

results['entropy_scaling'] = S_scaling

# Fit c from S_mid = (c/6)*ln(Lx) + const (cylinder scaling)
lnL = np.log([d['Lx'] for d in S_scaling])
S_vals = [d['S_mid'] for d in S_scaling]
# Use only Lx >= 10 to avoid small-size effects
mask = [d['Lx'] >= 10 for d in S_scaling]
lnL_fit = np.array([l for l, m in zip(lnL, mask) if m])
S_fit = np.array([s for s, m in zip(S_vals, mask) if m])
coeffs = np.polyfit(lnL_fit, S_fit, 1)
c_eff = 6 * coeffs[0]
print(f"\n  c_eff = {c_eff:.4f} (expect 0.5 for Ising)", flush=True)
results['c_eff'] = c_eff

# Save
with open("results/sprint_071a_cylinder_q2.json", "w") as f:
    json.dump(results, f, indent=2, default=str)
print("\nSaved to results/sprint_071a_cylinder_q2.json", flush=True)

from db_utils import record
record(sprint=71, model='hybrid_cyl', q=2, n=2,
       quantity='gc', value=gc_est, method='entropy_peak_Ly2',
       notes='Ly=2 cylinder chi=20')
record(sprint=71, model='hybrid_cyl', q=2, n=2,
       quantity='c', value=c_eff, method='S_vs_lnLx_Ly2',
       notes='Ly=2 cylinder chi=20')

print(f"\n--- Summary ---", flush=True)
print(f"g_c(1D) = 0.250, g_c(2D) = 0.771", flush=True)
print(f"g_c(Ly=2 cyl) = {gc_est:.4f}", flush=True)
print(f"c_eff = {c_eff:.4f} (Ising = 0.5)", flush=True)
