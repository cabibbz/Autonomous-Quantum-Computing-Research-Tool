#!/usr/bin/env python3
"""Sprint 050a: Find true q=3 Potts critical point via entropy FSS.

Sprint 049e gave n=8,12 data. Now add n=16,24,32 at dense g near critical region.
At g_c, S(n) = (c/6)ln(n) + const with c=4/5 for q=3 Potts.
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
        J = model_params.get('J', 1.0)
        g = model_params.get('g', 1.0)
        q = model_params.get('q', 3)
        for a in range(q):
            self.add_coupling(-J, 0, f'P{a}', 0, f'P{a}', 1)
        self.add_onsite(-g, 0, 'Xphc')


def run_dmrg(q, n, g, chi_max=60):
    t0 = time.time()
    model = PottsChain({'L': n, 'q': q, 'J': 1.0, 'g': g, 'bc_MPS': 'finite'})
    np.random.seed(42 + n + int(g * 1000))
    init = [np.random.randint(q) for _ in range(n)]
    psi = MPS.from_product_state(model.lat.mps_sites(), init, bc='finite')
    eng = dmrg.TwoSiteDMRGEngine(psi, model, {
        'mixer': True, 'max_E_err': 1e-12,
        'trunc_params': {'chi_max': chi_max, 'svd_min': 1e-14},
        'max_sweeps': 50,
    })
    E0, _ = eng.run()
    S = float(psi.entanglement_entropy()[n // 2 - 1])
    chi_used = max(psi.chi)
    return float(E0), S, chi_used, time.time() - t0


# Load existing n=8,12 data from Sprint 049e
try:
    with open('results/sprint_049e_potts_critical.json') as f:
        old = json.load(f)
    print("Loaded Sprint 049e data for n=8,12", flush=True)
    existing = {}
    for nk, data in old['data'].items():
        existing[int(nk)] = {d['g']: d['S_half'] for d in data}
except:
    existing = {}

# Dense g scan in the critical region
g_values = [0.18, 0.20, 0.22, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29,
            0.30, 0.32, 0.34, 0.36, 0.40, 0.50]

results = {}

# Time a single point first
print("=== Timing test: n=16, g=0.26, chi=60 ===", flush=True)
E, S, chi, dt = run_dmrg(3, 16, 0.26, chi_max=60)
print(f"  S={S:.6f}, chi={chi}, t={dt:.1f}s", flush=True)

# Decide chi based on timing
if dt > 60:
    chi_large = 40
    print(f"  Using chi=40 for n>=16 (too slow at chi=60)", flush=True)
else:
    chi_large = 60
    print(f"  Using chi=60 for all sizes", flush=True)

# Run n=16, 24, 32
for n in [16, 24, 32]:
    chi = chi_large if n <= 24 else min(chi_large, 40)
    print(f"\n=== q=3 Potts n={n}, chi={chi} ===", flush=True)
    results[n] = []
    t_start = time.time()
    for g in g_values:
        if time.time() - t_start > 220:
            print(f"  Time limit reached at n={n}", flush=True)
            break
        E, S, chi_used, dt = run_dmrg(3, n, g, chi_max=chi)
        print(f"  g={g:.2f}: S={S:.6f}, chi={chi_used}, t={dt:.1f}s", flush=True)
        results[n].append({'g': float(g), 'E': float(E), 'S_half': float(S),
                          'chi_used': chi_used, 'time': dt})

    # Save incrementally
    with open('results/sprint_050a_q3_gc.json', 'w') as f:
        json.dump({'sprint': '050a', 'q': 3, 'method': 'DMRG',
                   'data': {str(k): v for k, v in results.items()}}, f, indent=2)

# === Central charge analysis ===
print("\n=== Central charge c(g) from S(n) = (c/6)ln(n) + const ===", flush=True)

# Combine existing + new data
all_data = {}
for n, gdict in existing.items():
    for g, s in gdict.items():
        g_r = round(g, 2)
        if g_r not in all_data:
            all_data[g_r] = {}
        all_data[g_r][n] = s
for n, data_list in results.items():
    for d in data_list:
        g_r = round(d['g'], 2)
        if g_r not in all_data:
            all_data[g_r] = {}
        all_data[g_r][n] = d['S_half']

print(f"\n{'g':>6s} | {'c':>6s} | sizes | S values", flush=True)
print("-" * 70, flush=True)

c_values = []
for g in sorted(all_data.keys()):
    ns = sorted(all_data[g].keys())
    if len(ns) < 3:
        continue
    ns_a = np.array(ns, dtype=float)
    Ss_a = np.array([all_data[g][n] for n in ns])
    # Fit S = (c/6) ln(n) + const
    coeffs = np.polyfit(np.log(ns_a), Ss_a, 1)
    c = 6 * coeffs[0]
    s_str = ", ".join([f"S({n})={all_data[g][n]:.4f}" for n in ns])
    print(f"  {g:.2f} | {c:>6.3f} | {ns} | {s_str}", flush=True)
    c_values.append((g, c))

# Find where c is closest to 4/5 = 0.800
if c_values:
    best_g, best_c = min(c_values, key=lambda x: abs(x[1] - 0.800))
    print(f"\n*** Best match to c=4/5: g={best_g:.2f}, c={best_c:.3f} ***", flush=True)

    # Also check pairwise c estimates for convergence
    print("\n=== Pairwise c estimates at best g ===", flush=True)
    ns = sorted(all_data[best_g].keys())
    for i in range(len(ns) - 1):
        n1, n2 = ns[i], ns[i+1]
        S1, S2 = all_data[best_g][n1], all_data[best_g][n2]
        c_pair = 6 * (S2 - S1) / (np.log(n2) - np.log(n1))
        print(f"  c({n1},{n2}) = {c_pair:.4f}", flush=True)

print("\nDone!", flush=True)
