#!/usr/bin/env python3
"""Sprint 065b: Clock model at q=7 — first characterization.
REWRITTEN for speed: coarse scan n=4 only, n=6 only near crossing.
"""
import numpy as np
from scipy.sparse import csr_matrix, kron as sp_kron, eye as sp_eye, diags
from scipy.sparse.linalg import eigsh
import json, time

q = 7

def clock_hamiltonian_periodic(n, q, g):
    dim = q**n
    clock_2site = np.zeros(q**2)
    for a in range(q):
        for b in range(q):
            clock_2site[a*q + b] = np.cos(2*np.pi*(a-b)/q)
    clock_op = diags(clock_2site, 0, shape=(q**2, q**2), format='csr')
    X = np.zeros((q, q))
    for s in range(q):
        X[(s+1) % q, s] = 1.0
    XpXd = csr_matrix(X + X.T)
    H = csr_matrix((dim, dim))
    for i in range(n - 1):
        left = q**i; right = q**(n - i - 2)
        op = clock_op
        if left > 1: op = sp_kron(sp_eye(left), op, format='csr')
        if right > 1: op = sp_kron(op, sp_eye(right), format='csr')
        H = H - op
    diag_boundary = np.zeros(dim)
    for idx in range(dim):
        s0 = idx % q; sn = (idx // (q**(n-1))) % q
        diag_boundary[idx] = np.cos(2*np.pi*(sn - s0)/q)
    H = H - diags(diag_boundary, 0, shape=(dim, dim), format='csr')
    for i in range(n):
        left = q**i; right = q**(n - i - 1)
        op = XpXd
        if left > 1: op = sp_kron(sp_eye(left), op, format='csr')
        if right > 1: op = sp_kron(op, sp_eye(right), format='csr')
        H = H - g * op
    return H

def get_gap_n(n, q, g):
    H = clock_hamiltonian_periodic(n, q, g)
    k_eig = min(4, q**n - 2)
    evals, _ = eigsh(H, k=k_eig, which='SA')
    evals = np.sort(evals)
    return (evals[1] - evals[0]) * n, evals

# ============================================================
# STEP 1: Fast g_c scan with n=4 only
# ============================================================
print("=" * 60)
print(f"CLOCK q={q}: COARSE g_c SCAN (n=4 only)")
print("=" * 60)

t0_total = time.time()
gap4 = {}
for g in np.arange(0.30, 0.80, 0.02):
    g_f = round(float(g), 3)
    t0 = time.time()
    gap_val, _ = get_gap_n(4, q, g_f)
    gap4[g_f] = gap_val
    print(f"  g={g_f:.3f}: Δ·4={gap_val:.6f} ({time.time()-t0:.1f}s)")

# Find gap minimum (approximate g_c for n=4)
g_min = min(gap4, key=gap4.get)
print(f"\n  n=4 gap minimum at g ≈ {g_min:.3f}")

# Now compute n=6 at 5 points near the minimum
print(f"\n  Computing n=6 near g_min={g_min:.3f}...")
gap6 = {}
for g in np.arange(g_min - 0.04, g_min + 0.06, 0.02):
    g_f = round(float(g), 3)
    t0 = time.time()
    gap_val, _ = get_gap_n(6, q, g_f)
    gap6[g_f] = gap_val
    dt = time.time() - t0
    g4_val = gap4.get(g_f, None)
    if g4_val is None:
        g4_val, _ = get_gap_n(4, q, g_f)
        gap4[g_f] = g4_val
    print(f"  g={g_f:.3f}: Δ·4={g4_val:.6f}, Δ·6={gap_val:.6f} ({dt:.1f}s)")

# Find crossing
gc_clock = None
common_g = sorted(set(gap4.keys()) & set(gap6.keys()))
for i in range(len(common_g) - 1):
    g1, g2 = common_g[i], common_g[i+1]
    diff1 = gap4[g1] - gap6[g1]
    diff2 = gap4[g2] - gap6[g2]
    if diff1 * diff2 < 0:
        gc_clock = g1 + (-diff1) * (g2 - g1) / (diff2 - diff1)
        print(f"\n  Crossing at g ≈ {gc_clock:.4f}")
        break

if gc_clock is None:
    gc_clock = g_min
    print(f"\n  No crossing found. Using n=4 minimum: g ≈ {gc_clock:.3f}")

gc_corrected = gc_clock * 1.048
print(f"  Raw: {gc_clock:.4f}, FSS-corrected: {gc_corrected:.4f}")
print(f"  Scan time: {time.time()-t0_total:.1f}s")

results = {"q": q, "model": "clock", "gc_raw": float(gc_clock), "gc_corrected": float(gc_corrected)}

# ============================================================
# STEP 2: c/x₁ from spectrum at g_c
# ============================================================
print(f"\n{'='*60}")
print(f"c/x₁ AT g_c = {gc_corrected:.4f}")
print(f"{'='*60}")

E0_per_site = {}
delta1_n = {}

for n in [4, 6]:
    t0 = time.time()
    gap_val, evals = get_gap_n(n, q, gc_corrected)
    dt = time.time() - t0
    E0_per_site[n] = evals[0] / n
    delta1_n[n] = gap_val
    ratios = [(evals[i]-evals[0])/(evals[1]-evals[0]) for i in range(1, min(4, len(evals)))]
    print(f"  n={n}: E₀/N={E0_per_site[n]:.8f}, Δ·N={delta1_n[n]:.6f} ({dt:.1f}s)")
    print(f"    Ratios: {[f'{r:.3f}' for r in ratios]}")

n1, n2 = 4, 6
e_diff = E0_per_site[n1] - E0_per_site[n2]
inv_diff = 1.0/n1**2 - 1.0/n2**2
vc = 6 * e_diff / (np.pi * inv_diff)
vx1 = delta1_n[n2] / (2 * np.pi)
cx1_ratio = abs(vc / vx1)
print(f"\n  c/x₁ ({n1},{n2}): {cx1_ratio:.4f}")

results["cx1_ratio"] = float(cx1_ratio)
results["E0_per_site"] = {str(n): float(v) for n, v in E0_per_site.items()}
results["delta1_n"] = {str(n): float(v) for n, v in delta1_n.items()}

# ============================================================
# STEP 3: DMRG c at n=8
# ============================================================
print(f"\n{'='*60}")
print("DMRG c AT n=8")
print(f"{'='*60}")

import warnings
warnings.filterwarnings('ignore')
from tenpy.models.model import CouplingMPOModel, NearestNeighborModel
from tenpy.networks.site import Site
from tenpy.linalg import np_conserved as npc
from tenpy.networks.mps import MPS
from tenpy.algorithms import dmrg

class ClockSite(Site):
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

class ClockChainModel(CouplingMPOModel, NearestNeighborModel):
    def init_sites(self, model_params):
        return ClockSite(model_params.get('q', 7))
    def init_terms(self, model_params):
        J = model_params.get('J', 1.0)
        g = model_params.get('g', 1.0)
        q = model_params.get('q', 7)
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

t0 = time.time()
model = ClockChainModel({'L': 8, 'q': q, 'J': 1.0, 'g': gc_corrected, 'bc_MPS': 'finite'})
np.random.seed(42 + 8 + q)
init = [np.random.randint(q) for _ in range(8)]
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
c_val = extract_c_profile(S_prof, 8)
print(f"  n=8: c={c_val:.4f}, chi={chi_act}, time={dt:.1f}s")
results["dmrg_n8"] = {"c": c_val, "chi": chi_act, "time": round(dt, 1)}

# ============================================================
# COMPARISON
# ============================================================
print(f"\n{'='*60}")
print("CLOCK vs HYBRID at q=7")
print(f"{'='*60}")

hybrid_gc = 0.535
hybrid_cx1 = 15.109
hybrid_c = 1.462
hybrid_x1 = hybrid_c / hybrid_cx1

clock_c = c_val
clock_x1 = clock_c / cx1_ratio

print(f"\n{'Quantity':<15} {'Hybrid':>12} {'Clock':>12} {'Diff %':>10}")
print("-" * 52)
print(f"{'g_c':<15} {hybrid_gc:>12.4f} {gc_corrected:>12.4f} {100*(gc_corrected-hybrid_gc)/hybrid_gc:>+10.1f}%")
print(f"{'c/x₁':<15} {hybrid_cx1:>12.4f} {cx1_ratio:>12.4f} {100*(cx1_ratio-hybrid_cx1)/hybrid_cx1:>+10.1f}%")
print(f"{'c (DMRG n=8)':<15} {hybrid_c:>12.4f} {clock_c:>12.4f} {100*(clock_c-hybrid_c)/hybrid_c:>+10.1f}%")
print(f"{'x₁':<15} {hybrid_x1:>12.4f} {clock_x1:>12.4f} {100*(clock_x1-hybrid_x1)/hybrid_x1:>+10.1f}%")
print(f"{'c·x₁':<15} {hybrid_c*hybrid_x1:>12.5f} {clock_c*clock_x1:>12.5f}")

results["comparison"] = {
    "hybrid_gc": hybrid_gc, "hybrid_cx1": hybrid_cx1,
    "hybrid_c": hybrid_c, "hybrid_x1": float(hybrid_x1),
    "clock_gc": float(gc_corrected), "clock_cx1": float(cx1_ratio),
    "clock_c": float(clock_c), "clock_x1": float(clock_x1),
}

with open("results/sprint_065b_clock_q7.json", "w") as f:
    json.dump(results, f, indent=2, default=str)
print("\nResults saved to results/sprint_065b_clock_q7.json")
