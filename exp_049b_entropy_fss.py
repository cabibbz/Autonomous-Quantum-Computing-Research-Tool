#!/usr/bin/env python3
"""Sprint 049b: Half-chain entropy FSS for ν extraction.

Fixed: use random initial MPS to avoid DMRG getting stuck in ordered state.
Focus on q=2 (validation) + q=3,4 (where ν matters most).
q=5,7,10 are too slow without charge conservation for a full sweep.
"""
import numpy as np
import json, time, sys
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


def run_entropy(q, n, g, chi_max=40):
    """Run DMRG and return half-chain entropy + energy.

    Uses alternating initial state to break symmetry and help convergence.
    """
    t0 = time.time()
    if q == 2:
        model = TFIChain({'L': n, 'J': 1.0, 'g': g, 'bc_MPS': 'finite'})
        # Alternating initial state: up-down-up-down
        init = [i % 2 for i in range(n)]
    else:
        model = PottsChain({'L': n, 'q': q, 'J': 1.0, 'g': g, 'bc_MPS': 'finite'})
        # Alternating initial state in q-state space
        init = [i % q for i in range(n)]

    psi = MPS.from_product_state(model.lat.mps_sites(), init, bc='finite')
    eng = dmrg.TwoSiteDMRGEngine(psi, model, {
        'mixer': True, 'max_E_err': 1e-10,
        'trunc_params': {'chi_max': chi_max, 'svd_min': 1e-12},
        'max_sweeps': 30,
    })
    E0, _ = eng.run()
    dt = time.time() - t0

    S_half = float(psi.entanglement_entropy()[n//2 - 1])
    return float(E0), S_half, dt


def entropy_fss_collapse(nu, g_c, g_arr, S_arr, n_arr):
    """Collapse quality for S(g, n)."""
    sizes = sorted(set(n_arr))
    if len(sizes) < 2:
        return 1e10

    curves = {}
    for s in sizes:
        mask = (n_arr == s)
        x = (g_arr[mask] - g_c) * s**(1.0/nu)
        y = S_arr[mask]
        order = np.argsort(x)
        curves[s] = (x[order], y[order])

    total_err = 0.0
    count = 0
    for i, s1 in enumerate(sizes):
        for s2 in sizes[i+1:]:
            x1, y1 = curves[s1]
            x2, y2 = curves[s2]
            x_min = max(x1.min(), x2.min())
            x_max = min(x1.max(), x2.max())
            if x_min >= x_max:
                continue
            x_grid = np.linspace(x_min, x_max, 50)
            y1_i = np.interp(x_grid, x1, y1)
            y2_i = np.interp(x_grid, x2, y2)
            diff = (y1_i - y2_i)**2
            norm = (y1_i**2 + y2_i**2) / 2
            total_err += np.sum(diff / (norm + 1e-15))
            count += len(x_grid)
    return total_err / max(count, 1)


# === Quick test: q=3 convergence with alternating init ===
print("=== Convergence test: q=3, n=16, g=1.0, chi=40 ===", flush=True)
E0, S, dt = run_entropy(3, 16, 1.0, chi_max=40)
print(f"  E0={E0:.6f}, S_half={S:.4f}, time={dt:.1f}s", flush=True)
if S < 0.1:
    print("  WARNING: S still very low — DMRG may not be converging", flush=True)

# ================================================================
# q=2 TFIM: thorough FSS (fast, well-converged)
# ================================================================
print("\n" + "="*60, flush=True)
print("q=2 TFIM: Entropy FSS", flush=True)
print("="*60, flush=True)

g_values_q2 = [0.70, 0.80, 0.85, 0.90, 0.93, 0.95, 0.97, 0.99,
               1.00, 1.01, 1.03, 1.05, 1.07, 1.10, 1.15, 1.20, 1.30, 1.50]
sizes_q2 = [16, 24, 32, 48, 64]

results_q2 = {}
for n in sizes_q2:
    print(f"\n  n={n}:", flush=True)
    results_q2[n] = []
    for g in g_values_q2:
        t0 = time.time()
        E0, S, dt = run_entropy(2, n, g, chi_max=40)
        print(f"    g={g:.2f}: S={S:.4f}, E={E0:.6f}, t={dt:.1f}s", flush=True)
        results_q2[n].append({'g': g, 'E0': E0, 'S_half': S, 'time': dt})

    # Save incrementally
    with open('results/sprint_049b_entropy_fss.json', 'w') as f:
        json.dump({'sprint': '049b', 'q2': {str(k): v for k, v in results_q2.items()}}, f, indent=2)

# q=2 Analysis
print("\n=== q=2 Analysis ===", flush=True)
g_all, S_all, n_all = [], [], []
for n, data in results_q2.items():
    for d in data:
        g_all.append(d['g']); S_all.append(d['S_half']); n_all.append(n)
g_all, S_all, n_all = np.array(g_all), np.array(S_all), np.array(n_all)

# Peak of dS/dg
for n in sizes_q2:
    mask = (n_all == n)
    gs = g_all[mask]; Ss = S_all[mask]
    order = np.argsort(gs); gs, Ss = gs[order], Ss[order]
    dSdg = np.gradient(Ss, gs)
    peak = np.argmax(np.abs(dSdg))
    print(f"  n={n}: S peak slope at g={gs[peak]:.3f}, dS/dg={dSdg[peak]:.3f}, S(g_c)={np.interp(1.0,gs,Ss):.4f}", flush=True)

# Central charge
S_at_gc = []
for n in sizes_q2:
    mask = (n_all == n)
    gs = g_all[mask]; Ss = S_all[mask]
    order = np.argsort(gs)
    S_gc = np.interp(1.0, gs[order], Ss[order])
    S_at_gc.append((n, S_gc))
ns = np.array([x[0] for x in S_at_gc])
Ss_gc = np.array([x[1] for x in S_at_gc])
c_fit = 6 * np.polyfit(np.log(ns), Ss_gc, 1)[0]
print(f"  Central charge c = {c_fit:.3f} (exact=0.500 for Ising)", flush=True)

# FSS collapse
nu_scan = np.linspace(0.3, 3.0, 55)
qualities = [entropy_fss_collapse(nu, 1.0, g_all, S_all, n_all) for nu in nu_scan]
best_idx = np.argmin(qualities)
print(f"  Best ν = {nu_scan[best_idx]:.2f} (quality={qualities[best_idx]:.6f})", flush=True)
print(f"  ν=1.0: quality={entropy_fss_collapse(1.0, 1.0, g_all, S_all, n_all):.6f}", flush=True)
print(f"  ν=0.5: quality={entropy_fss_collapse(0.5, 1.0, g_all, S_all, n_all):.6f}", flush=True)

# ================================================================
# q=3 Potts: moderate FSS
# ================================================================
print("\n" + "="*60, flush=True)
print("q=3 Potts: Entropy FSS", flush=True)
print("="*60, flush=True)

g_values_q3 = [0.70, 0.80, 0.85, 0.90, 0.93, 0.95, 0.97, 0.99,
               1.00, 1.01, 1.03, 1.05, 1.10, 1.20, 1.40]
sizes_q3 = [12, 16, 24]

results_q3 = {}
for n in sizes_q3:
    print(f"\n  n={n}:", flush=True)
    results_q3[n] = []
    t_start = time.time()
    for g in g_values_q3:
        if time.time() - t_start > 300:
            print(f"    Time limit reached", flush=True)
            break
        t0 = time.time()
        E0, S, dt = run_entropy(3, n, g, chi_max=40)
        print(f"    g={g:.2f}: S={S:.4f}, E={E0:.6f}, t={dt:.1f}s", flush=True)
        results_q3[n].append({'g': g, 'E0': E0, 'S_half': S, 'time': dt})

    with open('results/sprint_049b_entropy_fss.json', 'w') as f:
        json.dump({'sprint': '049b',
                   'q2': {str(k): v for k, v in results_q2.items()},
                   'q3': {str(k): v for k, v in results_q3.items()}}, f, indent=2)

# q=3 Analysis
print("\n=== q=3 Analysis ===", flush=True)
g_all3, S_all3, n_all3 = [], [], []
for n, data in results_q3.items():
    for d in data:
        g_all3.append(d['g']); S_all3.append(d['S_half']); n_all3.append(n)
g_all3, S_all3, n_all3 = np.array(g_all3), np.array(S_all3), np.array(n_all3)

for n in sizes_q3:
    mask = (n_all3 == n)
    gs = g_all3[mask]; Ss = S_all3[mask]
    order = np.argsort(gs); gs, Ss = gs[order], Ss[order]
    if len(gs) >= 3:
        dSdg = np.gradient(Ss, gs)
        peak = np.argmax(np.abs(dSdg))
        S_gc = np.interp(1.0, gs, Ss)
        print(f"  n={n}: peak at g={gs[peak]:.3f}, dS/dg={dSdg[peak]:.3f}, S(g_c)={S_gc:.4f}", flush=True)

if len(g_all3) >= 6:
    nu_scan = np.linspace(0.3, 3.0, 55)
    qualities = [entropy_fss_collapse(nu, 1.0, g_all3, S_all3, n_all3) for nu in nu_scan]
    best_idx = np.argmin(qualities)
    print(f"  Best ν = {nu_scan[best_idx]:.2f} (expected 5/6≈0.833)", flush=True)
    print(f"  ν=5/6: quality={entropy_fss_collapse(5/6, 1.0, g_all3, S_all3, n_all3):.6f}", flush=True)
    print(f"  ν=1.0: quality={entropy_fss_collapse(1.0, 1.0, g_all3, S_all3, n_all3):.6f}", flush=True)

# ================================================================
# q=4 Potts: the BKT candidate
# ================================================================
print("\n" + "="*60, flush=True)
print("q=4 Potts: Entropy FSS (BKT candidate)", flush=True)
print("="*60, flush=True)

g_values_q4 = [0.50, 0.60, 0.70, 0.75, 0.80, 0.83, 0.86, 0.89,
               0.92, 0.95, 1.00, 1.10, 1.20]
sizes_q4 = [8, 12, 16]

results_q4 = {}
for n in sizes_q4:
    print(f"\n  n={n}:", flush=True)
    results_q4[n] = []
    t_start = time.time()
    for g in g_values_q4:
        if time.time() - t_start > 300:
            print(f"    Time limit reached", flush=True)
            break
        t0 = time.time()
        E0, S, dt = run_entropy(4, n, g, chi_max=40)
        print(f"    g={g:.2f}: S={S:.4f}, E={E0:.6f}, t={dt:.1f}s", flush=True)
        results_q4[n].append({'g': g, 'E0': E0, 'S_half': S, 'time': dt})

    with open('results/sprint_049b_entropy_fss.json', 'w') as f:
        json.dump({'sprint': '049b',
                   'q2': {str(k): v for k, v in results_q2.items()},
                   'q3': {str(k): v for k, v in results_q3.items()},
                   'q4': {str(k): v for k, v in results_q4.items()}}, f, indent=2)

# q=4 Analysis
print("\n=== q=4 Analysis ===", flush=True)
g_all4, S_all4, n_all4 = [], [], []
for n, data in results_q4.items():
    for d in data:
        g_all4.append(d['g']); S_all4.append(d['S_half']); n_all4.append(n)
g_all4, S_all4, n_all4 = np.array(g_all4), np.array(S_all4), np.array(n_all4)

for n in sizes_q4:
    mask = (n_all4 == n)
    gs = g_all4[mask]; Ss = S_all4[mask]
    order = np.argsort(gs); gs, Ss = gs[order], Ss[order]
    if len(gs) >= 3:
        dSdg = np.gradient(Ss, gs)
        peak = np.argmax(np.abs(dSdg))
        print(f"  n={n}: peak at g={gs[peak]:.3f}, dS/dg={dSdg[peak]:.3f}", flush=True)

if len(g_all4) >= 6:
    # Test power law vs BKT
    for nu_test, label in [(1.0, 'Ising'), (5/6, 'Potts q=3'), (2.0, 'large-ν'), (3.0, 'BKT-proxy')]:
        q_val = entropy_fss_collapse(nu_test, 0.89, g_all4, S_all4, n_all4)
        print(f"  ν={nu_test:.2f} ({label}): quality={q_val:.6f}", flush=True)

    nu_scan = np.linspace(0.3, 5.0, 95)
    qualities = [entropy_fss_collapse(nu, 0.89, g_all4, S_all4, n_all4) for nu in nu_scan]
    best_idx = np.argmin(qualities)
    print(f"  Best ν = {nu_scan[best_idx]:.2f} (quality={qualities[best_idx]:.6f})", flush=True)

print("\nDone! Results saved to results/sprint_049b_entropy_fss.json", flush=True)
