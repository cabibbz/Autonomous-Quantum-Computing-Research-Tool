#!/usr/bin/env python3
"""Sprint 050b: q=4 Potts entropy via DMRG at n=8,12.

q=4 exact diag n=6 (dim=4096) shows steepest drop at g≈0.35.
q=4 exact diag n=8 (dim=65536) is too slow.
Use DMRG for n=8,12 to get c(g) from entropy FSS and locate g_c.
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


# g values — dense near the transition (g≈0.30-0.40 based on n=6 data)
g_values = [0.15, 0.20, 0.25, 0.28, 0.30, 0.32, 0.34, 0.36, 0.38, 0.40, 0.45, 0.50, 0.70]

# Load q=4 n=6 exact diag data from exp_050a3
try:
    with open('results/sprint_050a3_gc_verify.json') as f:
        old = json.load(f)
    n6_data = {round(d['g'], 2): d['S'] for d in old['data'].get('q4_n6', [])}
    print(f"Loaded q=4 n=6 exact diag data: {len(n6_data)} points", flush=True)
except:
    n6_data = {}

# Timing test
print("=== Timing test: q=4 n=8, g=0.32, chi=60 ===", flush=True)
E, S, chi, dt = run_dmrg(4, 8, 0.32, chi_max=60)
print(f"  S={S:.6f}, chi={chi}, t={dt:.1f}s", flush=True)

results = {}

for n in [8, 12]:
    chi = 60 if n == 8 else 40
    print(f"\n=== q=4 Potts n={n}, chi={chi} (DMRG) ===", flush=True)
    results[n] = []
    t_start = time.time()
    for g in g_values:
        if time.time() - t_start > 230:
            print(f"  Time limit reached at n={n}", flush=True)
            break
        E, S, chi_used, dt = run_dmrg(4, n, g, chi_max=chi)
        print(f"  g={g:.2f}: S={S:.6f}, chi={chi_used}, t={dt:.1f}s", flush=True)
        results[n].append({'g': float(g), 'E': float(E), 'S_half': float(S),
                          'chi_used': chi_used, 'time': dt})

    # Save incrementally
    with open('results/sprint_050b_q4_dmrg.json', 'w') as f:
        json.dump({'sprint': '050b', 'q': 4, 'method': 'DMRG',
                   'data': {str(k): v for k, v in results.items()}}, f, indent=2)

# === Analysis: central charge c(g) ===
print("\n=== Central charge c(g) for q=4 Potts ===", flush=True)

all_data = {}
# n=6 exact diag
for g, s in n6_data.items():
    if g not in all_data:
        all_data[g] = {}
    all_data[g][6] = s
# DMRG
for n, data_list in results.items():
    for d in data_list:
        g_r = round(d['g'], 2)
        if g_r not in all_data:
            all_data[g_r] = {}
        all_data[g_r][n] = d['S_half']

print(f"\n{'g':>6s} | {'c(6,8)':>7s} | {'c(8,12)':>7s} | S(6)     S(8)     S(12)", flush=True)
print("-" * 70, flush=True)

c_results = []
for g in sorted(all_data.keys()):
    ns = sorted(all_data[g].keys())
    Ss = {n: all_data[g][n] for n in ns}
    c_68, c_812 = None, None
    if 6 in Ss and 8 in Ss:
        c_68 = 6 * (Ss[8] - Ss[6]) / (np.log(8) - np.log(6))
    if 8 in Ss and 12 in Ss:
        c_812 = 6 * (Ss[12] - Ss[8]) / (np.log(12) - np.log(8))
    s_str = "  ".join([f"S({n})={Ss[n]:.4f}" for n in ns])
    c_68_str = f"{c_68:.3f}" if c_68 is not None else "  --  "
    c_812_str = f"{c_812:.3f}" if c_812 is not None else "  --  "
    print(f"  {g:.2f} | {c_68_str:>7s} | {c_812_str:>7s} | {s_str}", flush=True)
    if c_68 is not None or c_812 is not None:
        c_results.append((g, c_68, c_812))

# dS/dg analysis for n=8
if results.get(8):
    gs = [d['g'] for d in results[8]]
    Ss = [d['S_half'] for d in results[8]]
    dSdg = np.gradient(Ss, gs)
    min_idx = np.argmin(dSdg)
    print(f"\n  n=8 steepest drop: g={gs[min_idx]:.2f}, dS/dg={dSdg[min_idx]:.3f}", flush=True)

print("\nDone!", flush=True)
