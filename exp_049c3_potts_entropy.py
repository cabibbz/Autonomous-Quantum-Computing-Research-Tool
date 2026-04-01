#!/usr/bin/env python3
"""Sprint 049c3: q=3 and q=4 Potts entropy FSS with improved convergence.

Key fix: use random initial MPS (not alternating = ordered state!) and
increase max_sweeps + mixer strength for Potts at criticality.
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


def run_potts_entropy(q, n, g, chi_max=40):
    """DMRG with random-ish initial state to avoid getting stuck in ordered state."""
    t0 = time.time()
    model = PottsChain({'L': n, 'q': q, 'J': 1.0, 'g': g, 'bc_MPS': 'finite'})
    # Use a non-trivial initial state: alternate between 0 and 1 only
    # This breaks the q-fold symmetry without being the fully ordered state
    np.random.seed(42 + n + int(g * 1000))
    init = [np.random.randint(q) for _ in range(n)]
    psi = MPS.from_product_state(model.lat.mps_sites(), init, bc='finite')
    eng = dmrg.TwoSiteDMRGEngine(psi, model, {
        'mixer': True, 'max_E_err': 1e-10,
        'trunc_params': {'chi_max': chi_max, 'svd_min': 1e-12},
        'max_sweeps': 40,
    })
    E0, _ = eng.run()
    S = float(psi.entanglement_entropy()[n//2 - 1])
    return float(E0), S, time.time() - t0


# === Test convergence fix ===
print("=== Convergence test: q=3 n=16 g=1.0 (critical) ===", flush=True)
E, S, dt = run_potts_entropy(3, 16, 1.0, 40)
print(f"  E={E:.6f}, S={S:.4f}, t={dt:.1f}s", flush=True)
print(f"  Expected S > 0.3 for critical q=3 Potts", flush=True)

print("=== Convergence test: q=3 n=12 g=1.0 ===", flush=True)
E, S, dt = run_potts_entropy(3, 12, 1.0, 40)
print(f"  E={E:.6f}, S={S:.4f}, t={dt:.1f}s", flush=True)

# Load existing
with open('results/sprint_049b_entropy_fss.json') as f:
    results = json.load(f)

def save():
    with open('results/sprint_049b_entropy_fss.json', 'w') as f:
        json.dump(results, f, indent=2)

# ================================================================
# q=3 Potts: n=12, 16 (skip n=24 — too slow at ~110s/point)
# ================================================================
g_values_q3 = [0.70, 0.80, 0.85, 0.90, 0.93, 0.95, 0.97, 0.99,
               1.00, 1.01, 1.03, 1.05, 1.10, 1.20, 1.40]

for n in [12, 16]:
    print(f"\n=== q=3 n={n} ===", flush=True)
    n_data = []
    t_start = time.time()
    for g in g_values_q3:
        if time.time() - t_start > 200:
            print(f"  Time limit", flush=True)
            break
        E0, S, dt = run_potts_entropy(3, n, g, chi_max=40)
        print(f"  g={g:.2f}: S={S:.4f}, E={E0:.6f}, t={dt:.1f}s", flush=True)
        n_data.append({'g': g, 'E0': E0, 'S_half': S, 'time': dt})
    if 'q3' not in results:
        results['q3'] = {}
    results['q3'][str(n)] = n_data
    save()

# ================================================================
# q=4 Potts: n=8, 12 (n=16 at 60s/point is borderline)
# ================================================================
g_values_q4 = [0.50, 0.60, 0.70, 0.75, 0.80, 0.83, 0.86, 0.89,
               0.92, 0.95, 1.00, 1.10, 1.20]

for n in [8, 12, 16]:
    print(f"\n=== q=4 n={n} ===", flush=True)
    n_data = []
    t_start = time.time()
    tlimit = 120 if n <= 12 else 200
    for g in g_values_q4:
        if time.time() - t_start > tlimit:
            print(f"  Time limit", flush=True)
            break
        E0, S, dt = run_potts_entropy(4, n, g, chi_max=40)
        print(f"  g={g:.2f}: S={S:.4f}, E={E0:.6f}, t={dt:.1f}s", flush=True)
        n_data.append({'g': g, 'E0': E0, 'S_half': S, 'time': dt})
    if 'q4' not in results:
        results['q4'] = {}
    results['q4'][str(n)] = n_data
    save()

print("\nPotts entropy FSS saved!", flush=True)
