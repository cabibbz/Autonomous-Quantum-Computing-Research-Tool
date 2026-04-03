#!/usr/bin/env python3
"""Sprint 063b: Test c·x₁ for the Z_5 clock model (REWRITTEN for speed).

Clock model: H = -J Σ cos(2π(s_i-s_{i+1})/q) - g Σ (X + X†)
Build Hamiltonian using diagonal coupling matrix instead of q² projectors.

Plan:
1. Find g_c via energy gap Δ·N crossing (n=4,6)
2. Extract c/x₁ from Casimir energy + gap
3. Get c from DMRG entropy profile
4. Compute c·x₁
"""
import numpy as np
from scipy.sparse import csr_matrix, kron as sp_kron, eye as sp_eye, diags
from scipy.sparse.linalg import eigsh
import json, time

def clock_hamiltonian_periodic(n, q, g):
    """Clock model with periodic BC. Efficient construction."""
    dim = q**n
    eye_q = sp_eye(q, format='csr')

    # Clock coupling matrix (q×q): M[a,b] = cos(2π(a-b)/q)
    # This is diagonal in a×b basis: the 2-site operator is just the diagonal
    clock_2site = np.zeros(q**2)
    for a in range(q):
        for b in range(q):
            clock_2site[a*q + b] = np.cos(2*np.pi*(a-b)/q)
    clock_op = diags(clock_2site, 0, shape=(q**2, q**2), format='csr')

    # X + X† matrix
    X = np.zeros((q, q))
    for s in range(q):
        X[(s+1) % q, s] = 1.0
    XpXd = csr_matrix(X + X.T)

    H = csr_matrix((dim, dim))

    # Bulk bonds (i, i+1) for i=0..n-2
    for i in range(n - 1):
        left = q**i
        right = q**(n - i - 2)
        op = clock_op
        if left > 1:
            op = sp_kron(sp_eye(left), op, format='csr')
        if right > 1:
            op = sp_kron(op, sp_eye(right), format='csr')
        H = H - op

    # Boundary bond (n-1, 0) — build via permutation
    # Need to contract sites 0 and n-1 which are far apart in tensor product
    # Use explicit diagonal: for each basis state, add cos(2π(s_{n-1}-s_0)/q)
    diag_boundary = np.zeros(dim)
    for idx in range(dim):
        s0 = idx % q
        sn = (idx // (q**(n-1))) % q
        diag_boundary[idx] = np.cos(2*np.pi*(sn - s0)/q)
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


q = 5
print("=" * 70)
print(f"CLOCK q={q}: ENERGY GAP g_c SCAN")
print("=" * 70)

# Coarse scan
g_vals = np.arange(0.50, 1.00, 0.05)
gap_data = {4: {}, 6: {}}

t0_total = time.time()
for g in g_vals:
    for n in [4, 6]:
        dim = q**n
        t0 = time.time()
        H = clock_hamiltonian_periodic(n, q, g)
        evals, _ = eigsh(H, k=4, which='SA')
        evals = np.sort(evals)
        gap_data[n][g] = (evals[1] - evals[0]) * n
        dt = time.time() - t0

    print(f"  g={g:.3f}: Δ·4={gap_data[4][g]:.6f}, Δ·6={gap_data[6][g]:.6f}")

# Find crossing
gc_clock = None
g_list = sorted(gap_data[4].keys())
for i in range(len(g_list)-1):
    g1, g2 = g_list[i], g_list[i+1]
    diff1 = gap_data[4][g1] - gap_data[6][g1]
    diff2 = gap_data[4][g2] - gap_data[6][g2]
    if diff1 * diff2 < 0:
        gc_clock = g1 + (-diff1) * (g2 - g1) / (diff2 - diff1)
        print(f"\n  Coarse crossing at g ≈ {gc_clock:.4f}")

# Refine
if gc_clock:
    g_fine = np.arange(gc_clock - 0.05, gc_clock + 0.05, 0.005)
    for g in g_fine:
        for n in [4, 6]:
            H = clock_hamiltonian_periodic(n, q, g)
            evals, _ = eigsh(H, k=4, which='SA')
            evals = np.sort(evals)
            gap_data[n][g] = (evals[1] - evals[0]) * n

    g_list = sorted([g for g in gap_data[4].keys() if gc_clock - 0.06 < g < gc_clock + 0.06])
    for i in range(len(g_list)-1):
        g1, g2 = g_list[i], g_list[i+1]
        diff1 = gap_data[4][g1] - gap_data[6][g1]
        diff2 = gap_data[4][g2] - gap_data[6][g2]
        if diff1 * diff2 < 0:
            gc_clock = g1 + (-diff1) * (g2 - g1) / (diff2 - diff1)
            print(f"  Refined crossing at g_c ≈ {gc_clock:.4f}")
            break

gc_corrected = gc_clock * 1.048  # FSS correction
print(f"\n  Raw crossing: {gc_clock:.4f}")
print(f"  FSS-corrected (+4.8%): {gc_corrected:.4f}")

print(f"\n  Total scan time: {time.time()-t0_total:.1f}s")

# Extract c/x₁ at both g_c values
print(f"\n{'='*70}")
print("c/x₁ EXTRACTION FROM SPECTRUM")
print(f"{'='*70}")

results = {"q": q, "model": "clock", "gc_raw": float(gc_clock), "gc_corrected": float(gc_corrected)}

for gc_label, gc_use in [("raw", gc_clock), ("corrected", gc_corrected)]:
    print(f"\n--- g_c = {gc_use:.4f} ({gc_label}) ---")

    E0_per_site = {}
    delta1_n = {}

    for n in [4, 6, 8]:
        dim = q**n
        if dim > 500_000:
            continue
        t0 = time.time()
        H = clock_hamiltonian_periodic(n, q, gc_use)
        k_eig = min(6, dim-2)
        evals, _ = eigsh(H, k=k_eig, which='SA')
        evals = np.sort(evals)
        dt = time.time() - t0

        E0_per_site[n] = evals[0] / n
        delta1_n[n] = (evals[1] - evals[0]) * n
        print(f"  n={n}: E₀/N={E0_per_site[n]:.8f}, Δ₁·N={delta1_n[n]:.6f} ({dt:.1f}s)")

    sizes = sorted(E0_per_site.keys())
    for i in range(len(sizes)):
        for j in range(i+1, len(sizes)):
            n1, n2 = sizes[i], sizes[j]
            e_diff = E0_per_site[n1] - E0_per_site[n2]
            inv_diff = 1.0/n1**2 - 1.0/n2**2
            vc = 6 * e_diff / (np.pi * inv_diff)
            vx1 = delta1_n[n2] / (2 * np.pi)
            cx1_ratio = vc / vx1
            print(f"  Pair ({n1},{n2}): c/x₁ = {cx1_ratio:.3f}")

    results[f"E0_per_site_{gc_label}"] = {str(n): float(v) for n, v in E0_per_site.items()}
    results[f"delta1_n_{gc_label}"] = {str(n): float(v) for n, v in delta1_n.items()}

# DMRG entropy profile for c
print(f"\n{'='*70}")
print("c FROM DMRG ENTROPY PROFILE")
print(f"{'='*70}")

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

class ClockChain(CouplingMPOModel, NearestNeighborModel):
    def init_sites(self, model_params):
        return ClockSite(model_params.get('q', 5))
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

for gc_label, gc_use in [("raw", gc_clock)]:
    print(f"\n--- g_c = {gc_use:.4f} ({gc_label}) ---")
    for n_dmrg in [12, 16, 24]:
        t0 = time.time()
        model = ClockChain({'L': n_dmrg, 'q': q, 'J': 1.0, 'g': gc_use, 'bc_MPS': 'finite'})
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
        results[f"dmrg_n{n_dmrg}_{gc_label}"] = {"c": c_val, "chi": chi_act, "time": round(dt, 1)}

        if dt > 120:
            print("  Too slow, stopping DMRG.")
            break

# Summary
print(f"\n\n{'='*70}")
print("SUMMARY: CLOCK q=5 c·x₁")
print(f"{'='*70}")

# Get best c/x₁ from spectrum
for gc_label in ["raw"]:
    key_e = f"E0_per_site_{gc_label}"
    key_d = f"delta1_n_{gc_label}"
    sizes = sorted(results[key_e].keys(), key=int)
    if len(sizes) >= 2:
        n1, n2 = int(sizes[-2]), int(sizes[-1])
        e1 = results[key_e][str(n1)]
        e2 = results[key_e][str(n2)]
        vc = 6 * (e1 - e2) / (np.pi * (1/n1**2 - 1/n2**2))
        vx1 = results[key_d][str(n2)] / (2 * np.pi)
        cx1_ratio = vc / vx1
        results[f"cx1_ratio_{gc_label}"] = float(cx1_ratio)

        # Get best c from DMRG
        for n_key in [24, 16, 12]:
            dmrg_key = f"dmrg_n{n_key}_{gc_label}"
            if dmrg_key in results:
                c_clock = results[dmrg_key]["c"]
                x1_clock = c_clock / cx1_ratio
                cx1_clock = c_clock * x1_clock
                print(f"\nClock q=5 (g_c={results[f'gc_{gc_label}']:.4f}):")
                print(f"  c = {c_clock:.4f} (DMRG n={n_key})")
                print(f"  c/x₁ = {cx1_ratio:.3f} (spectrum ({n1},{n2}))")
                print(f"  x₁ = {x1_clock:.4f}")
                print(f"  c·x₁ = {cx1_clock:.5f}")
                results[f"c_clock"] = c_clock
                results[f"x1_clock"] = float(x1_clock)
                results[f"cx1_clock"] = float(cx1_clock)
                break

# Compare
print(f"\n--- Comparison ---")
print(f"  Potts q=5: c=1.10, x₁=0.1015, c·x₁=0.1117")
print(f"  Target: c·x₁ ≈ 0.112 ± 0.005")

with open("results/sprint_063b_clock_cx1.json", "w") as f:
    json.dump(results, f, indent=2, default=str)
print("\nResults saved to results/sprint_063b_clock_cx1.json")
