#!/usr/bin/env python3
"""Sprint 065a: Head-to-head hybrid vs clock at q=5.

Compare c/x₁ ratio at SAME sizes (n=4,6,8) for both models.
If c/x₁ converges to same value → same universality class.
If different → distinct universality classes.

Also extract c from DMRG entropy profile at n=16,24.
"""
import numpy as np
from scipy.sparse import csr_matrix, kron as sp_kron, eye as sp_eye, diags
from scipy.sparse.linalg import eigsh
import json, time

q = 5

# ============================================================
# HAMILTONIAN BUILDERS
# ============================================================

def hybrid_hamiltonian_periodic(n, q, g):
    """Potts-clock hybrid: H = -J*delta(s_i,s_j) - g*(X+X†)"""
    dim = q**n
    eye_q = sp_eye(q, format='csr')

    # Potts delta coupling (q×q): M[a,b] = delta(a,b)
    potts_2site = np.zeros(q**2)
    for a in range(q):
        potts_2site[a*q + a] = 1.0
    potts_op = diags(potts_2site, 0, shape=(q**2, q**2), format='csr')

    # X + X† matrix
    X = np.zeros((q, q))
    for s in range(q):
        X[(s+1) % q, s] = 1.0
    XpXd = csr_matrix(X + X.T)

    H = csr_matrix((dim, dim))

    # Bulk bonds
    for i in range(n - 1):
        left = q**i
        right = q**(n - i - 2)
        op = potts_op
        if left > 1:
            op = sp_kron(sp_eye(left), op, format='csr')
        if right > 1:
            op = sp_kron(op, sp_eye(right), format='csr')
        H = H - op

    # Boundary bond (n-1, 0)
    diag_boundary = np.zeros(dim)
    for idx in range(dim):
        s0 = idx % q
        sn = (idx // (q**(n-1))) % q
        diag_boundary[idx] = 1.0 if s0 == sn else 0.0
    H = H - diags(diag_boundary, 0, shape=(dim, dim), format='csr')

    # Transverse field
    for i in range(n):
        left = q**i
        right = q**(n - i - 1)
        op = XpXd
        if left > 1:
            op = sp_kron(sp_eye(left), op, format='csr')
        if right > 1:
            op = sp_kron(op, sp_eye(right), format='csr')
        H = H - g * op

    return H


def clock_hamiltonian_periodic(n, q, g):
    """Clock model: H = -J*cos(2π(s_i-s_j)/q) - g*(X+X†)"""
    dim = q**n

    # Clock coupling (q×q)
    clock_2site = np.zeros(q**2)
    for a in range(q):
        for b in range(q):
            clock_2site[a*q + b] = np.cos(2*np.pi*(a-b)/q)
    clock_op = diags(clock_2site, 0, shape=(q**2, q**2), format='csr')

    # X + X†
    X = np.zeros((q, q))
    for s in range(q):
        X[(s+1) % q, s] = 1.0
    XpXd = csr_matrix(X + X.T)

    H = csr_matrix((dim, dim))

    for i in range(n - 1):
        left = q**i
        right = q**(n - i - 2)
        op = clock_op
        if left > 1:
            op = sp_kron(sp_eye(left), op, format='csr')
        if right > 1:
            op = sp_kron(op, sp_eye(right), format='csr')
        H = H - op

    # Boundary
    diag_boundary = np.zeros(dim)
    for idx in range(dim):
        s0 = idx % q
        sn = (idx // (q**(n-1))) % q
        diag_boundary[idx] = np.cos(2*np.pi*(sn - s0)/q)
    H = H - diags(diag_boundary, 0, shape=(dim, dim), format='csr')

    for i in range(n):
        left = q**i
        right = q**(n - i - 1)
        op = XpXd
        if left > 1:
            op = sp_kron(sp_eye(left), op, format='csr')
        if right > 1:
            op = sp_kron(op, sp_eye(right), format='csr')
        H = H - g * op

    return H

# ============================================================
# SPECTRUM EXTRACTION
# ============================================================

def extract_cft_data(ham_fn, q, gc, sizes, label):
    """Extract E0/N, Δ₁·N, c/x₁ for a model at given sizes."""
    print(f"\n{'='*60}")
    print(f"{label}: q={q}, g_c={gc:.4f}")
    print(f"{'='*60}")

    E0_per_site = {}
    delta1_n = {}

    for n in sizes:
        dim = q**n
        if dim > 500_000:
            print(f"  n={n}: dim={dim} too large, skipping")
            continue
        t0 = time.time()
        H = ham_fn(n, q, gc)
        k_eig = min(6, dim - 2)
        evals, _ = eigsh(H, k=k_eig, which='SA')
        evals = np.sort(evals)
        dt = time.time() - t0

        E0_per_site[n] = evals[0] / n
        delta1_n[n] = (evals[1] - evals[0]) * n
        print(f"  n={n}: E₀/N={E0_per_site[n]:.8f}, Δ₁·N={delta1_n[n]:.6f} ({dt:.1f}s)")

    # Pairwise c/x₁
    cx1_ratios = {}
    ns = sorted(E0_per_site.keys())
    for i in range(len(ns)):
        for j in range(i+1, len(ns)):
            n1, n2 = ns[i], ns[j]
            e_diff = E0_per_site[n1] - E0_per_site[n2]
            inv_diff = 1.0/n1**2 - 1.0/n2**2
            vc = 6 * e_diff / (np.pi * inv_diff)
            vx1 = delta1_n[n2] / (2 * np.pi)
            ratio = abs(vc / vx1)
            cx1_ratios[(n1, n2)] = ratio
            print(f"  c/x₁ pair ({n1},{n2}): {ratio:.4f}")

    return {
        "E0_per_site": {str(n): float(v) for n, v in E0_per_site.items()},
        "delta1_n": {str(n): float(v) for n, v in delta1_n.items()},
        "cx1_ratios": {f"({n1},{n2})": float(v) for (n1, n2), v in cx1_ratios.items()},
    }

# ============================================================
# RUN BOTH MODELS
# ============================================================

results = {"q": q}

# Known g_c values
gc_hybrid = 0.441  # from Sprint 051
gc_clock = 0.52    # from Sprint 063

sizes = [4, 6, 8]

hybrid_data = extract_cft_data(hybrid_hamiltonian_periodic, q, gc_hybrid, sizes, "HYBRID (Potts δ + clock field)")
clock_data = extract_cft_data(clock_hamiltonian_periodic, q, gc_clock, sizes, "CLOCK (cos coupling + clock field)")

results["hybrid"] = hybrid_data
results["hybrid"]["gc"] = gc_hybrid
results["clock"] = clock_data
results["clock"]["gc"] = gc_clock

# ============================================================
# HEAD-TO-HEAD COMPARISON
# ============================================================

print(f"\n{'='*60}")
print("HEAD-TO-HEAD COMPARISON: q=5")
print(f"{'='*60}")
print(f"{'Pair':<12} {'Hybrid c/x₁':>14} {'Clock c/x₁':>14} {'Diff %':>10}")
print("-" * 52)

for pair_key in sorted(hybrid_data["cx1_ratios"].keys()):
    h_val = hybrid_data["cx1_ratios"].get(pair_key, None)
    c_val = clock_data["cx1_ratios"].get(pair_key, None)
    if h_val and c_val:
        diff_pct = 100 * (c_val - h_val) / h_val
        print(f"{pair_key:<12} {h_val:>14.4f} {c_val:>14.4f} {diff_pct:>+10.1f}%")

# ============================================================
# DMRG ENTROPY PROFILE FOR c
# ============================================================

print(f"\n{'='*60}")
print("DMRG ENTROPY PROFILE FOR c")
print(f"{'='*60}")

import warnings
warnings.filterwarnings('ignore')
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

class HybridChain(CouplingMPOModel, NearestNeighborModel):
    def init_sites(self, model_params):
        return PottsSite(model_params.get('q', 5))
    def init_terms(self, model_params):
        J = model_params.get('J', 1.0)
        g = model_params.get('g', 1.0)
        q = model_params.get('q', 5)
        for a in range(q):
            self.add_coupling(-J, 0, f'P{a}', 0, f'P{a}', 1)
        self.add_onsite(-g, 0, 'Xphc')

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

def run_dmrg_c(model_class, q, gc, n_vals, label):
    """Extract c from DMRG entropy profile for multiple sizes."""
    print(f"\n--- {label}: q={q}, g_c={gc:.4f} ---")
    c_results = {}
    for n_dmrg in n_vals:
        t0 = time.time()
        model = model_class({'L': n_dmrg, 'q': q, 'J': 1.0, 'g': gc, 'bc_MPS': 'finite'})
        np.random.seed(42 + n_dmrg + q)
        init = [np.random.randint(q) for _ in range(n_dmrg)]
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

        c_val = extract_c_profile(S_prof, n_dmrg)
        print(f"  n={n_dmrg}: c={c_val:.4f}, chi={chi_act}, time={dt:.1f}s", flush=True)
        c_results[n_dmrg] = {"c": c_val, "chi": chi_act, "time": round(dt, 1)}

        if dt > 120:
            print("  Too slow, stopping DMRG.")
            break
    return c_results

# Run DMRG for both models
hybrid_dmrg = run_dmrg_c(HybridChain, q, gc_hybrid, [16, 24], "HYBRID")
clock_dmrg = run_dmrg_c(ClockChain, q, gc_clock, [16, 24], "CLOCK")

results["hybrid"]["dmrg"] = {str(n): v for n, v in hybrid_dmrg.items()}
results["clock"]["dmrg"] = {str(n): v for n, v in clock_dmrg.items()}

# ============================================================
# FINAL COMPARISON TABLE
# ============================================================

print(f"\n{'='*60}")
print("FINAL COMPARISON TABLE: q=5")
print(f"{'='*60}")

# Best c/x₁ (largest pair)
h_best_pair = max(hybrid_data["cx1_ratios"].keys())
c_best_pair = max(clock_data["cx1_ratios"].keys())
h_cx1 = hybrid_data["cx1_ratios"][h_best_pair]
c_cx1 = clock_data["cx1_ratios"][c_best_pair]

# Best c from DMRG
h_c_best = max(hybrid_dmrg.keys())
c_c_best = max(clock_dmrg.keys())
h_c = hybrid_dmrg[h_c_best]["c"]
c_c = clock_dmrg[c_c_best]["c"]

# x₁ = c / (c/x₁)
h_x1 = h_c / h_cx1
c_x1 = c_c / c_cx1

print(f"\n{'Quantity':<15} {'Hybrid':>12} {'Clock':>12} {'Diff %':>10}")
print("-" * 52)
print(f"{'g_c':<15} {gc_hybrid:>12.4f} {gc_clock:>12.4f} {100*(gc_clock-gc_hybrid)/gc_hybrid:>+10.1f}%")
print(f"{'c/x₁ best':<15} {h_cx1:>12.4f} {c_cx1:>12.4f} {100*(c_cx1-h_cx1)/h_cx1:>+10.1f}%")
print(f"{'c (DMRG)':<15} {h_c:>12.4f} {c_c:>12.4f} {100*(c_c-h_c)/h_c:>+10.1f}%")
print(f"{'x₁':<15} {h_x1:>12.4f} {c_x1:>12.4f} {100*(c_x1-h_x1)/h_x1:>+10.1f}%")
print(f"{'c·x₁':<15} {h_c*h_x1:>12.5f} {c_c*c_x1:>12.5f} {100*(c_c*c_x1-h_c*h_x1)/(h_c*h_x1):>+10.1f}%")

results["comparison"] = {
    "hybrid_cx1_ratio": float(h_cx1),
    "clock_cx1_ratio": float(c_cx1),
    "hybrid_c": float(h_c),
    "clock_c": float(c_c),
    "hybrid_x1": float(h_x1),
    "clock_x1": float(c_x1),
}

# Convergence analysis: does the difference shrink with n?
print(f"\n--- Convergence analysis: c/x₁ difference vs size ---")
for pair_key in sorted(hybrid_data["cx1_ratios"].keys()):
    h_val = hybrid_data["cx1_ratios"].get(pair_key, None)
    c_val = clock_data["cx1_ratios"].get(pair_key, None)
    if h_val and c_val:
        diff_pct = 100 * abs(c_val - h_val) / max(h_val, c_val)
        print(f"  {pair_key}: |diff| = {diff_pct:.1f}%")

with open("results/sprint_065a_hybrid_vs_clock_q5.json", "w") as f:
    json.dump(results, f, indent=2, default=str)
print("\nResults saved to results/sprint_065a_hybrid_vs_clock_q5.json")
