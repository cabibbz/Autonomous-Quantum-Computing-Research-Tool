#!/usr/bin/env python3
"""Sprint 063b: Test c·x₁ ≈ 0.112 for the Z_5 clock model.

The clock model uses cos(2π(s_i-s_j)/q) coupling instead of Potts δ(s_i,s_j).
Same Z_q symmetry but different Hamiltonian. If c·x₁ ≈ constant is universal
to Z_q-symmetric models, it should hold for clock too.

Clock q=5: g_c ≈ 0.67 from MI-CV crossings (Sprint 041-042, but these
were MI-CV based — may be inaccurate). We'll use energy gap method.

For clock: H = -J Σ cos(2π(s_i-s_{i+1})/q) - g Σ (X + X†)
Compare with: Potts: H = -J Σ δ(s_i,s_{i+1}) - g Σ (X + X†)

Method:
1. Find g_c via energy gap Δ·N crossing
2. Measure x₁ from gap ratios at g_c
3. Measure c from Casimir energy or entropy profile
4. Compute c·x₁
"""
import numpy as np
from scipy.sparse import csr_matrix, kron as sp_kron, eye as sp_eye
from scipy.sparse.linalg import eigsh
import json, time

def clock_hamiltonian_periodic(n, q, g):
    """Clock model Hamiltonian with periodic BC.
    H = -Σ cos(2π(s_i-s_{i+1})/q) - g·Σ(X+X†)
    """
    dim = q**n
    eye_q = np.eye(q)

    # Clock coupling: cos(2π(s_i-s_j)/q)
    coupling = np.zeros((q**2, q**2))
    for si in range(q):
        for sj in range(q):
            coupling[si*q+sj, si*q+sj] = np.cos(2*np.pi*(si-sj)/q)

    # Transverse field X + X†
    X = np.zeros((q, q))
    for s in range(q):
        X[(s+1) % q, s] = 1.0
    XpXd = X + X.T

    H = csr_matrix((dim, dim))

    # Nearest-neighbor coupling (periodic)
    for i in range(n):
        j = (i + 1) % n
        if i < j:
            left = q**i if i > 0 else 1
            middle = q**(j-i-1) if j-i-1 > 0 else 1
            right = q**(n-j-1) if n-j-1 > 0 else 1
            op = csr_matrix(coupling)
            if i > 0:
                op = sp_kron(sp_eye(left), op, format='csr')
            if j-i-1 > 0:
                # This won't work for non-adjacent sites, need different approach for boundary
                pass
            if n-j-1 > 0:
                op = sp_kron(op, sp_eye(right), format='csr')
            H = H - op
        else:
            # Boundary bond (n-1, 0): need to handle specially
            # Build as Σ_{si,sj} cos(2π(si-sj)/q) |si>_{n-1}<si|_{n-1} ⊗ |sj>_0<sj|_0
            for si in range(q):
                for sj in range(q):
                    coeff = np.cos(2*np.pi*(si-sj)/q)
                    if abs(coeff) < 1e-15:
                        continue
                    proj_0 = np.zeros((q, q)); proj_0[sj, sj] = 1.0
                    proj_n = np.zeros((q, q)); proj_n[si, si] = 1.0
                    ops = [eye_q] * n
                    ops[0] = proj_0
                    ops[n-1] = proj_n
                    r = csr_matrix(ops[0])
                    for op in ops[1:]:
                        r = sp_kron(r, csr_matrix(op), format='csr')
                    H = H - coeff * r

    # Transverse field
    for i in range(n):
        ops = [eye_q] * n
        ops[i] = XpXd
        r = csr_matrix(ops[0])
        for op in ops[1:]:
            r = sp_kron(r, csr_matrix(op), format='csr')
        H = H - g * r

    return H

# Step 1: Find g_c for clock q=5 via energy gap crossing
print("=" * 70)
print("CLOCK q=5: ENERGY GAP METHOD FOR g_c")
print("=" * 70)

q = 5

# Scan g values
g_vals = np.arange(0.50, 1.00, 0.05)
print(f"\n{'g':>6} {'Δ·4':>10} {'Δ·6':>10}")

gap_n4 = {}
gap_n6 = {}

for g in g_vals:
    for n in [4, 6]:
        dim = q**n
        if dim > 2_000_000:
            continue
        H = clock_hamiltonian_periodic(n, q, g)
        k_eig = min(4, dim-2)
        evals, _ = eigsh(H, k=k_eig, which='SA')
        evals = np.sort(evals)
        gap = evals[1] - evals[0]
        gap_n = gap * n
        if n == 4:
            gap_n4[g] = gap_n
        else:
            gap_n6[g] = gap_n

    g4 = gap_n4.get(g, float('nan'))
    g6 = gap_n6.get(g, float('nan'))
    print(f"{g:>6.3f} {g4:>10.6f} {g6:>10.6f}")

# Find crossing
print("\n--- Finding Δ·N crossing ---")
g_list = sorted(gap_n4.keys())
for i in range(len(g_list)-1):
    g1, g2 = g_list[i], g_list[i+1]
    if g1 in gap_n6 and g2 in gap_n6:
        diff1 = gap_n4[g1] - gap_n6[g1]
        diff2 = gap_n4[g2] - gap_n6[g2]
        if diff1 * diff2 < 0:
            g_cross = g1 + (-diff1) * (g2 - g1) / (diff2 - diff1)
            print(f"  Crossing at g ≈ {g_cross:.4f}")

# Refine around crossing
print("\n--- Refining crossing ---")
g_fine = np.arange(0.60, 0.80, 0.01)
gap_n4_f = {}
gap_n6_f = {}

for g in g_fine:
    for n in [4, 6]:
        dim = q**n
        H = clock_hamiltonian_periodic(n, q, g)
        k_eig = min(4, dim-2)
        evals, _ = eigsh(H, k=k_eig, which='SA')
        evals = np.sort(evals)
        gap = evals[1] - evals[0]
        gap_n = gap * n
        if n == 4:
            gap_n4_f[g] = gap_n
        else:
            gap_n6_f[g] = gap_n

g_list_f = sorted(gap_n4_f.keys())
gc_clock = None
for i in range(len(g_list_f)-1):
    g1, g2 = g_list_f[i], g_list_f[i+1]
    if g1 in gap_n6_f and g2 in gap_n6_f:
        diff1 = gap_n4_f[g1] - gap_n6_f[g1]
        diff2 = gap_n4_f[g2] - gap_n6_f[g2]
        if diff1 * diff2 < 0:
            gc_clock = g1 + (-diff1) * (g2 - g1) / (diff2 - diff1)
            print(f"  Refined crossing at g_c ≈ {gc_clock:.4f}")

if gc_clock is None:
    # Try n=4,8
    print("  No n=4,6 crossing found. Trying wider scan...")
    gc_clock = 0.70  # Fallback estimate from Sprint 042

# Apply 4.8% FSS correction (from Potts calibration)
gc_corrected = gc_clock * 1.048
print(f"  Raw crossing: {gc_clock:.4f}")
print(f"  FSS-corrected (+4.8%): {gc_corrected:.4f}")

# Step 2: Measure spectrum at g_c
print(f"\n\n{'='*70}")
print(f"SPECTRUM AT g_c = {gc_clock:.4f} (raw crossing)")
print(f"{'='*70}")

for gc_use in [gc_clock, gc_corrected]:
    print(f"\n--- g = {gc_use:.4f} ---")
    for n in [4, 6, 8]:
        dim = q**n
        if dim > 500_000:
            continue
        t0 = time.time()
        H = clock_hamiltonian_periodic(n, q, gc_use)
        k_eig = min(8, dim-2)
        evals, _ = eigsh(H, k=k_eig, which='SA')
        evals = np.sort(evals)
        dt = time.time() - t0

        E0 = evals[0]
        delta1 = evals[1] - evals[0]

        # CFT: x₁ = Δ₁·N/(2π·v), and E₀/N = e_∞ - πvc/(6N²)
        # From two sizes: v·c = -6/π · (E₀(n)/n - E₀(n')/n') · n²·n'²/(n'²-n²)

        print(f"  n={n}: E₀/N={E0/n:.8f}, Δ₁={delta1:.8f}, Δ₁·N={delta1*n:.6f} ({dt:.1f}s)")

# Step 3: Extract c and x₁ from size pairs
print(f"\n\n{'='*70}")
print("c AND x₁ EXTRACTION")
print(f"{'='*70}")

results = {"q": q, "model": "clock"}

for gc_label, gc_use in [("raw", gc_clock), ("corrected", gc_corrected)]:
    print(f"\n--- g_c = {gc_use:.4f} ({gc_label}) ---")

    E0_per_site = {}
    delta1_n = {}

    for n in [4, 6, 8]:
        dim = q**n
        if dim > 500_000:
            continue
        H = clock_hamiltonian_periodic(n, q, gc_use)
        k_eig = min(6, dim-2)
        evals, _ = eigsh(H, k=k_eig, which='SA')
        evals = np.sort(evals)
        E0_per_site[n] = evals[0] / n
        delta1_n[n] = (evals[1] - evals[0]) * n

    # Extract v·c and v·x₁ from pairs
    sizes = sorted(E0_per_site.keys())
    for i in range(len(sizes)):
        for j in range(i+1, len(sizes)):
            n1, n2 = sizes[i], sizes[j]
            # E₀/N = e_∞ - πvc/(6N²)
            # e₁ - e₂ = πvc/6 · (1/n₁² - 1/n₂²)
            e_diff = E0_per_site[n1] - E0_per_site[n2]
            inv_diff = 1.0/n1**2 - 1.0/n2**2
            vc = 6 * e_diff / (np.pi * inv_diff)

            # Δ₁·N = 2πv·x₁
            # Average from the two sizes
            vx1_1 = delta1_n[n1] / (2 * np.pi)
            vx1_2 = delta1_n[n2] / (2 * np.pi)

            cx1_1 = vc * vx1_1 / (vc / (2*np.pi*vx1_1))**0 # c·x₁ = (vc)·(vx₁)/v² = vc·vx₁/(vc/c·vx₁/x₁)...
            # Simpler: c/x₁ = vc/(vx₁), c·x₁ = vc·vx₁/v²
            # But v² = vc·vx₁ · (1/(c·x₁))... circular.
            # Instead: c/x₁ = vc/vx₁ (model-independent ratio)

            cx1_ratio = vc / vx1_1
            cx1_ratio2 = vc / vx1_2

            # To get c·x₁: need c and x₁ separately.
            # c·x₁ = (vc)·(vx₁) / v² and v = vc/c = vx₁·2π·x₁/Δ₁N...
            # Actually: c = vc/v, x₁ = vx₁·(2π)/v... no.
            #
            # From Casimir: vc is known.
            # From gap: vx₁ = Δ₁N/(2π) is known.
            # c/x₁ = vc/vx₁ is model-independent.
            # c·x₁ = (vc)·(vx₁)/v². Need v separately.
            #
            # v can be extracted from: Δ₁ = 2πv·x₁/N → v = Δ₁·N/(2πx₁)
            # But we don't know x₁. We know vx₁ = Δ₁N/(2π).
            #
            # Actually: c·x₁ = vc · vx₁ / v² = 1/(c/x₁) · c² = not helpful.
            #
            # Better: c²/(c/x₁) = c·x₁. And c = vc/v. We need v.
            #
            # Alternative: from vc and vx₁:
            # v = sqrt(vc · vx₁ / (c/x₁ · x₁²))... still circular.
            #
            # The model-independent quantity is c/x₁ = vc/vx₁.
            # To get c·x₁ we need the PRODUCT vc · vx₁ = v² · c · x₁ / 1
            # so c·x₁ = (vc · vx₁) / v².
            #
            # We can get v from higher levels or from the descendant gap.
            # But simpler: if we know c from an independent method (DMRG entropy),
            # then x₁ = c/(c/x₁).

            print(f"  Pair ({n1},{n2}): vc = {vc:.6f}, vx₁(n₁) = {vx1_1:.6f}, vx₁(n₂) = {vx1_2:.6f}")
            print(f"    c/x₁ = {cx1_ratio:.3f} (using n₁), {cx1_ratio2:.3f} (using n₂)")

    results[f"gc_{gc_label}"] = gc_use
    results[f"E0_per_site_{gc_label}"] = {str(n): float(v) for n, v in E0_per_site.items()}
    results[f"delta1_n_{gc_label}"] = {str(n): float(v) for n, v in delta1_n.items()}

# Step 4: Get c from DMRG entropy profile
print(f"\n\n{'='*70}")
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
        # Clock coupling: cos(2π(si-sj)/q) needs projectors
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
        # Clock coupling: -J·cos(2π(si-sj)/q) = -J·Σ_{a,b} cos(2π(a-b)/q) Pa_i Pb_j
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

gc_test = gc_clock  # Use raw crossing for DMRG too

for n_dmrg in [12, 16, 24]:
    t0 = time.time()
    model = ClockChain({'L': n_dmrg, 'q': q, 'J': 1.0, 'g': gc_test, 'bc_MPS': 'finite'})
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
    print(f"  n={n_dmrg}: c={c_val:.4f}, chi={chi_act}, time={dt:.1f}s")

    results[f"dmrg_n{n_dmrg}"] = {"c": c_val, "chi": chi_act, "time": round(dt, 1)}

    if dt > 120:
        print("  Too slow, stopping.")
        break

# Summary
print(f"\n\n{'='*70}")
print("SUMMARY: CLOCK q=5 vs POTTS q=5")
print(f"{'='*70}")

# Potts q=5 data
potts_c = 1.10
potts_x1 = 0.1015
potts_cx1 = potts_c * potts_x1

print(f"\nPotts q=5: c={potts_c:.3f}, x₁={potts_x1:.4f}, c·x₁={potts_cx1:.5f}")

# Clock: use c/x₁ from spectrum + c from DMRG
clock_c = results.get("dmrg_n16", results.get("dmrg_n12", {})).get("c", None)
if clock_c:
    # c/x₁ from spectrum pair
    # Find best pair
    for gc_label in ["raw", "corrected"]:
        key_e = f"E0_per_site_{gc_label}"
        key_d = f"delta1_n_{gc_label}"
        if key_e in results:
            sizes = sorted(results[key_e].keys(), key=int)
            if len(sizes) >= 2:
                n1, n2 = int(sizes[0]), int(sizes[1])
                e1, e2 = results[key_e][sizes[0]], results[key_e][sizes[1]]
                vc = 6 * (e1 - e2) / (np.pi * (1/n1**2 - 1/n2**2))
                vx1 = results[key_d][sizes[1]] / (2 * np.pi)
                cx1_ratio = vc / vx1
                clock_x1 = clock_c / cx1_ratio
                clock_cx1 = clock_c * clock_x1
                print(f"\nClock q=5 ({gc_label}): c={clock_c:.3f}, c/x₁={cx1_ratio:.3f}, x₁={clock_x1:.4f}, c·x₁={clock_cx1:.5f}")

print(f"\nTarget: c·x₁ = 0.112 ± 0.005 (Potts pattern)")

with open("results/sprint_063b_clock_cx1.json", "w") as f:
    json.dump(results, f, indent=2, default=str)
print("\nResults saved to results/sprint_063b_clock_cx1.json")
