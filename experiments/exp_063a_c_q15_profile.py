#!/usr/bin/env python3
"""Sprint 063a: Measure c(q=15) via entropy profile to test c·x₁ ≈ 1/9.

Known: x₁(q=15) ≈ 0.071 (Sprint 061, n=4 descendant gap).
If c·x₁ = 1/9, then c(15) = 1/(9·0.071) ≈ 1.57.
Log formula predicts c(15) ≈ 0.40·ln(14) + 0.55 ≈ 1.61.

Use DMRG entropy profile at open BC. Start with small n to check timing.
"""
import numpy as np, json, time, warnings
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


class PottsChain(CouplingMPOModel, NearestNeighborModel):
    def init_sites(self, model_params):
        return PottsSite(model_params.get('q', 3))
    def init_terms(self, model_params):
        J = model_params.get('J', 1.0); g = model_params.get('g', 1.0); q = model_params.get('q', 3)
        for a in range(q):
            self.add_coupling(-J, 0, f'P{a}', 0, f'P{a}', 1)
        self.add_onsite(-g, 0, 'Xphc')


def run_dmrg(q, n, g, chi_max=30):
    t0 = time.time()
    model = PottsChain({'L': n, 'q': q, 'J': 1.0, 'g': g, 'bc_MPS': 'finite'})
    np.random.seed(42 + n + q)
    init = [np.random.randint(q) for _ in range(n)]
    psi = MPS.from_product_state(model.lat.mps_sites(), init, bc='finite')
    eng = dmrg.TwoSiteDMRGEngine(psi, model, {
        'mixer': True, 'max_E_err': 1e-10,
        'trunc_params': {'chi_max': chi_max, 'svd_min': 1e-12},
        'max_sweeps': 40,
    })
    E0, _ = eng.run()
    S_profile = [float(s) for s in psi.entanglement_entropy()]
    chi_actual = max(psi.chi)
    dt = time.time() - t0
    return float(E0), S_profile, chi_actual, dt


def extract_c(S_profile, n):
    """Extract central charge from entropy profile using chord distance."""
    ls = np.arange(1, n)
    l_chord = (2 * n / np.pi) * np.sin(np.pi * ls / n)
    ln_chord = np.log(l_chord)
    S = np.array(S_profile)
    # Use central half of chain
    quarter = n // 4
    mask = (ls >= quarter) & (ls <= 3 * quarter)
    A = np.vstack([ln_chord[mask], np.ones(mask.sum())]).T
    slope, _ = np.linalg.lstsq(A, S[mask], rcond=None)[0]
    return float(6 * slope)


# g_c formula
def g_c(q):
    return 0.2 * np.sqrt(q - 1) + 0.05

# Known data for validation
x1_known = {
    3: 0.1337, 5: 0.1015, 7: 0.0860, 10: 0.0831,
    15: 0.071, 20: 0.055, 25: 0.045, 30: 0.038
}
c_known = {2: 0.500, 3: 0.800, 4: 1.000, 5: 1.10, 7: 1.30, 10: 1.40}

print("=" * 70)
print("c(q=15) MEASUREMENT VIA ENTROPY PROFILE")
print("=" * 70)
print(f"g_c(15) = {g_c(15):.4f}")
print(f"x₁(15) = {x1_known[15]}")
print(f"If c·x₁ = 1/9: c = {1/(9*x1_known[15]):.4f}")
print(f"Log formula: c = {0.40*np.log(14) + 0.55:.4f}")

# First: timing test at small n
print("\n--- Timing test at n=8 ---")
E0, S_prof, chi_act, dt = run_dmrg(15, 8, g_c(15), chi_max=30)
c8 = extract_c(S_prof, 8)
print(f"  n=8: c={c8:.4f}, chi={chi_act}, time={dt:.1f}s", flush=True)

results = {"q": 15, "g_c": float(g_c(15)), "x1": x1_known[15]}
results["n8"] = {"c": c8, "chi": chi_act, "time": round(dt, 1), "E0": E0, "S_half": S_prof[3]}

# If timing ok, go to n=12
if dt < 60:
    print("\n--- n=12 ---")
    E0, S_prof, chi_act, dt = run_dmrg(15, 12, g_c(15), chi_max=30)
    c12 = extract_c(S_prof, 12)
    print(f"  n=12: c={c12:.4f}, chi={chi_act}, time={dt:.1f}s", flush=True)
    results["n12"] = {"c": c12, "chi": chi_act, "time": round(dt, 1), "E0": E0, "S_half": S_prof[5]}

    if dt < 60:
        print("\n--- n=16 ---")
        E0, S_prof, chi_act, dt = run_dmrg(15, 16, g_c(15), chi_max=30)
        c16 = extract_c(S_prof, 16)
        print(f"  n=16: c={c16:.4f}, chi={chi_act}, time={dt:.1f}s", flush=True)
        results["n16"] = {"c": c16, "chi": chi_act, "time": round(dt, 1), "E0": E0, "S_half": S_prof[7]}

        if dt < 60:
            print("\n--- n=24 ---")
            E0, S_prof, chi_act, dt = run_dmrg(15, 24, g_c(15), chi_max=30)
            c24 = extract_c(S_prof, 24)
            print(f"  n=24: c={c24:.4f}, chi={chi_act}, time={dt:.1f}s", flush=True)
            results["n24"] = {"c": c24, "chi": chi_act, "time": round(dt, 1), "E0": E0, "S_half": S_prof[11]}

# Also do q=20 if fast enough
print("\n\n" + "=" * 70)
print("c(q=20) MEASUREMENT VIA ENTROPY PROFILE")
print("=" * 70)
print(f"g_c(20) = {g_c(20):.4f}")
print(f"x₁(20) = {x1_known[20]}")
print(f"If c·x₁ = 1/9: c = {1/(9*x1_known[20]):.4f}")

print("\n--- n=8 ---")
E0, S_prof, chi_act, dt = run_dmrg(20, 8, g_c(20), chi_max=30)
c8_20 = extract_c(S_prof, 8)
print(f"  n=8: c={c8_20:.4f}, chi={chi_act}, time={dt:.1f}s", flush=True)
results["q20_n8"] = {"c": c8_20, "chi": chi_act, "time": round(dt, 1), "E0": E0}

if dt < 60:
    print("\n--- n=12 ---")
    E0, S_prof, chi_act, dt = run_dmrg(20, 12, g_c(20), chi_max=30)
    c12_20 = extract_c(S_prof, 12)
    print(f"  n=12: c={c12_20:.4f}, chi={chi_act}, time={dt:.1f}s", flush=True)
    results["q20_n12"] = {"c": c12_20, "chi": chi_act, "time": round(dt, 1), "E0": E0}

    if dt < 60:
        print("\n--- n=16 ---")
        E0, S_prof, chi_act, dt = run_dmrg(20, 16, g_c(20), chi_max=30)
        c16_20 = extract_c(S_prof, 16)
        print(f"  n=16: c={c16_20:.4f}, chi={chi_act}, time={dt:.1f}s", flush=True)
        results["q20_n16"] = {"c": c16_20, "chi": chi_act, "time": round(dt, 1), "E0": E0}

# Summary
print("\n\n" + "=" * 70)
print("SUMMARY: c·x₁ test")
print("=" * 70)
print(f"{'q':>3} {'n':>3} {'c_raw':>8} {'x₁':>7} {'c·x₁':>7} {'1/9':>7} {'diff%':>6}")

# Include old data for comparison
for q in [3, 4, 5, 7, 10]:
    c = c_known[q]
    x1 = x1_known.get(q, None)
    if x1:
        cx1 = c * x1
        print(f"{q:>3} {'best':>3} {c:>8.4f} {x1:>7.4f} {cx1:>7.5f} {1/9:>7.5f} {(cx1-1/9)/(1/9)*100:>6.1f}")

# Add new q=15 results
for key in sorted(results.keys()):
    if key.startswith("n") and "c" in results[key]:
        n_val = int(key[1:])
        c_val = results[key]["c"]
        cx1 = c_val * x1_known[15]
        print(f" 15 {n_val:>3} {c_val:>8.4f} {x1_known[15]:>7.4f} {cx1:>7.5f} {1/9:>7.5f} {(cx1-1/9)/(1/9)*100:>6.1f}")

for key in sorted(results.keys()):
    if key.startswith("q20") and "c" in results[key]:
        n_val = int(key.split("n")[1])
        c_val = results[key]["c"]
        cx1 = c_val * x1_known[20]
        print(f" 20 {n_val:>3} {c_val:>8.4f} {x1_known[20]:>7.4f} {cx1:>7.5f} {1/9:>7.5f} {(cx1-1/9)/(1/9)*100:>6.1f}")

with open("results/sprint_063a_c_q15_q20.json", "w") as f:
    json.dump(results, f, indent=2)
print("\nResults saved to results/sprint_063a_c_q15_q20.json")
