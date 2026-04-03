#!/usr/bin/env python3
"""Sprint 049e: Find q=3 Potts critical point precisely using DMRG.

Key finding: g_c ≈ 0.25-0.35, NOT 1.0. Scan the critical region at multiple n.
DMRG verified against exact diag at n=8.
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
    # Random initial state for better convergence
    np.random.seed(42 + n + int(g * 1000))
    init = [np.random.randint(q) for _ in range(n)]
    psi = MPS.from_product_state(model.lat.mps_sites(), init, bc='finite')
    eng = dmrg.TwoSiteDMRGEngine(psi, model, {
        'mixer': True, 'max_E_err': 1e-12,
        'trunc_params': {'chi_max': chi_max, 'svd_min': 1e-14},
        'max_sweeps': 40,
    })
    E0, _ = eng.run()
    S = float(psi.entanglement_entropy()[n // 2 - 1])
    return float(E0), S, time.time() - t0


# Dense scan near the true critical region
g_values = [0.10, 0.15, 0.20, 0.22, 0.24, 0.26, 0.28, 0.30,
            0.32, 0.34, 0.36, 0.38, 0.40, 0.45, 0.50, 0.60, 0.80, 1.00]

results = {}
for n in [8, 12, 16, 24]:
    print(f"\n=== q=3 Potts n={n}, chi=60 ===", flush=True)
    results[n] = []
    t_start = time.time()
    for g in g_values:
        if time.time() - t_start > 250:
            print(f"  Time limit", flush=True)
            break
        E, S, dt = run_dmrg(3, n, g, chi_max=60)
        print(f"  g={g:.2f}: S={S:.6f}, E={E:.6f}, t={dt:.1f}s", flush=True)
        results[n].append({'g': float(g), 'E': float(E), 'S_half': float(S), 'time': dt})

    # Save incrementally
    with open('results/sprint_049e_potts_critical.json', 'w') as f:
        json.dump({'sprint': '049e', 'q': 3, 'chi_max': 60,
                   'data': {str(k): v for k, v in results.items()}}, f, indent=2)

# === Analysis ===
print("\n=== Entropy FSS analysis ===", flush=True)
sizes = sorted(results.keys())
for g in g_values:
    line = f"  g={g:.2f}: "
    vals = []
    for n in sizes:
        matching = [d for d in results.get(n, []) if abs(d['g'] - g) < 0.001]
        if matching:
            vals.append((n, matching[0]['S_half']))
            line += f"S(n={n})={matching[0]['S_half']:.4f}  "
    if len(vals) >= 2:
        # Check if S increases or saturates with n
        Ss = [v[1] for v in vals]
        if Ss[-1] > Ss[0] * 1.02:
            line += " ↑ (grows)"
        elif Ss[-1] < Ss[0] * 0.98:
            line += " ↓ (decreases)"
        else:
            line += " = (saturated)"
    print(line, flush=True)

# Central charge estimation where S grows with n
print("\n=== Central charge from S(n) at g values where S grows ===", flush=True)
for g in g_values:
    ns = []
    Ss = []
    for n in sizes:
        matching = [d for d in results.get(n, []) if abs(d['g'] - g) < 0.001]
        if matching:
            ns.append(n)
            Ss.append(matching[0]['S_half'])
    if len(ns) >= 3:
        ns_a = np.array(ns, dtype=float)
        Ss_a = np.array(Ss)
        # Fit S = (c/6) ln(n) + const
        coeffs = np.polyfit(np.log(ns_a), Ss_a, 1)
        c = 6 * coeffs[0]
        if abs(c) > 0.01:
            print(f"  g={g:.2f}: c = {c:.3f}", flush=True)

print("\nDone!", flush=True)
