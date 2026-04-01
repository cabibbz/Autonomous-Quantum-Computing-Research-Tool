#!/usr/bin/env python3
"""Sprint 065a v2: Save existing spectrum data + clock DMRG at n=16.

Exp A produced spectrum data before getting stuck on DMRG n=24.
This script:
1. Records the already-computed spectrum results
2. Runs ONLY clock DMRG at n=16 (skip n=24)
3. Produces comparison table
"""
import numpy as np
import json, time

# ============================================================
# RECORDED DATA FROM EXP 065a (before timeout)
# ============================================================

results = {
    "q": 5,
    "hybrid": {
        "gc": 0.441,
        "E0_per_site": {"4": -1.25116432, "6": -1.23560454, "8": -1.23033159},
        "delta1_n": {"4": 0.472613, "6": 0.477179, "8": 0.483454},
        "cx1_ratios": {"(4,6)": 11.2693, "(4,8)": 11.0314, "(6,8)": 10.7697},
        "dmrg": {"16": {"c": 1.2746, "chi": 30, "time": 262.1}},
    },
    "clock": {
        "gc": 0.52,
        "E0_per_site": {"4": -1.41623960, "6": -1.39961599, "8": -1.39391728},
        "delta1_n": {"4": 0.595146, "6": 0.594940, "8": 0.596475},
        "cx1_ratios": {"(4,6)": 9.6566, "(4,8)": 9.5805, "(6,8)": 9.4339},
    },
}

# ============================================================
# RUN CLOCK DMRG at n=16
# ============================================================

print("=" * 60)
print("CLOCK q=5 DMRG ENTROPY PROFILE at n=16")
print("=" * 60)

import warnings
warnings.filterwarnings('ignore')
from tenpy.models.model import CouplingMPOModel, NearestNeighborModel
from tenpy.networks.site import Site
from tenpy.linalg import np_conserved as npc
from tenpy.networks.mps import MPS
from tenpy.algorithms import dmrg

q = 5

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

class ClockChain(CouplingMPOModel, NearestNeighborModel):
    def init_sites(self, model_params):
        return PottsSite(model_params.get('q', 5))
    def init_terms(self, model_params):
        J = model_params.get('J', 1.0)
        g = model_params.get('g', 1.0)
        q = model_params.get('q', 5)
        for a in range(q):
            for b in range(q):
                coeff = np.cos(2*np.pi*(a-b)/q)
                if abs(coeff) > 1e-15:
                    self.add_coupling(-J * coeff, 0, f'P{a}', 0, f'P{b}', 1)
        self.add_onsite(-g, 0, 'Xphc')

def extract_c_profile(S_profile, n):
    ls = np.arange(1, n)
    l_chord = (2 * n / np.pi) * np.sin(np.pi * ls / n)
    ln_chord = np.log(l_chord)
    S = np.array(S_profile)
    quarter = n // 4
    mask = (ls >= quarter) & (ls <= 3 * quarter)
    A = np.vstack([ln_chord[mask], np.ones(mask.sum())]).T
    slope, _ = np.linalg.lstsq(A, S[mask], rcond=None)[0]
    return float(6 * slope)

gc_clock = 0.52
t0 = time.time()
model = ClockChain({'L': 16, 'q': q, 'J': 1.0, 'g': gc_clock, 'bc_MPS': 'finite'})
np.random.seed(42 + 16 + q)
init = [np.random.randint(q) for _ in range(16)]
psi = MPS.from_product_state(model.lat.mps_sites(), init, bc='finite')
eng = dmrg.TwoSiteDMRGEngine(psi, model, {
    'mixer': True, 'max_E_err': 1e-10,
    'trunc_params': {'chi_max': 30, 'svd_min': 1e-12},
    'max_sweeps': 30,
})
E0, _ = eng.run()
S_prof = [float(s) for s in psi.entanglement_entropy()]
chi_act = max(psi.chi)
dt = time.time() - t0
c_clock = extract_c_profile(S_prof, 16)
print(f"  n=16: c={c_clock:.4f}, chi={chi_act}, time={dt:.1f}s")
results["clock"]["dmrg"] = {"16": {"c": c_clock, "chi": chi_act, "time": round(dt, 1)}}

# ============================================================
# COMPARISON TABLE
# ============================================================

print(f"\n{'='*60}")
print("FINAL COMPARISON: HYBRID vs CLOCK at q=5")
print(f"{'='*60}")

h_cx1 = results["hybrid"]["cx1_ratios"]["(6,8)"]
c_cx1 = results["clock"]["cx1_ratios"]["(6,8)"]
h_c = results["hybrid"]["dmrg"]["16"]["c"]
c_c = results["clock"]["dmrg"]["16"]["c"]
h_x1 = h_c / h_cx1
c_x1 = c_c / c_cx1

print(f"\n{'Quantity':<15} {'Hybrid':>12} {'Clock':>12} {'Diff %':>10}")
print("-" * 52)
print(f"{'g_c':<15} {'0.441':>12} {'0.52':>12} {'+17.9%':>10}")
print(f"{'c/x₁ (6,8)':<15} {h_cx1:>12.4f} {c_cx1:>12.4f} {100*(c_cx1-h_cx1)/h_cx1:>+10.1f}%")
print(f"{'c (n=16)':<15} {h_c:>12.4f} {c_c:>12.4f} {100*(c_c-h_c)/h_c:>+10.1f}%")
print(f"{'x₁':<15} {h_x1:>12.4f} {c_x1:>12.4f} {100*(c_x1-h_x1)/h_x1:>+10.1f}%")
print(f"{'c·x₁':<15} {h_c*h_x1:>12.5f} {c_c*c_x1:>12.5f} {100*(c_c*c_x1-h_c*h_x1)/(h_c*h_x1):>+10.1f}%")

# Convergence trend
print(f"\n--- c/x₁ convergence: does difference shrink? ---")
for pair in ["(4,6)", "(4,8)", "(6,8)"]:
    h = results["hybrid"]["cx1_ratios"][pair]
    c = results["clock"]["cx1_ratios"][pair]
    pct = 100 * abs(c - h) / max(h, c)
    print(f"  {pair}: hybrid={h:.4f}, clock={c:.4f}, |diff|={pct:.1f}%")

results["comparison"] = {
    "hybrid_cx1_68": float(h_cx1), "clock_cx1_68": float(c_cx1),
    "hybrid_c_n16": float(h_c), "clock_c_n16": float(c_c),
    "hybrid_x1": float(h_x1), "clock_x1": float(c_x1),
    "hybrid_cx1_product": float(h_c * h_x1),
    "clock_cx1_product": float(c_c * c_x1),
}

with open("results/sprint_065a_hybrid_vs_clock_q5.json", "w") as f:
    json.dump(results, f, indent=2, default=str)
print("\nResults saved to results/sprint_065a_hybrid_vs_clock_q5.json")
