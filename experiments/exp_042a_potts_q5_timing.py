"""
Sprint 042a: True q=5 Potts model — custom Kronecker-delta coupling.
Timing test at n=8, single g point. Validates model + measures runtime.

Potts Hamiltonian: H = -J Σ δ(s_i, s_j) - g Σ (X_i + X_i†)
where δ(s_i,s_j) = Σ_a |a><a|_i ⊗ |a><a|_j (Kronecker delta)
and X is the q-state shift operator.

Key difference from clock: Clock uses cos(2π(s_i-s_j)/q) (graduated),
Potts uses δ (binary same/different).
"""

import numpy as np
from scipy import linalg as la
import json, time, warnings, sys
warnings.filterwarnings('ignore')
import logging; logging.getLogger('tenpy').setLevel(logging.WARNING)

from tenpy.models.model import CouplingMPOModel, NearestNeighborModel
from tenpy.networks.site import Site
from tenpy.algorithms import dmrg
from tenpy.networks.mps import MPS
import tenpy.linalg.np_conserved as npc

t0 = time.time()
q = 5

# --- Custom Potts Site ---
class PottsSite(Site):
    def __init__(self, q, conserve=None):
        d = q
        leg = npc.LegCharge.from_trivial(d)
        names = [str(a) for a in range(d)]
        Site.__init__(self, leg, names, sort_charge=False)
        # Projectors P_a = |a><a| for Kronecker-delta coupling
        for a in range(d):
            P = np.zeros((d, d), dtype=complex)
            P[a, a] = 1.0
            self.add_op(f'P{a}', P)
        # Shift operator X and X+Xdagger (transverse field)
        X = np.zeros((d, d), dtype=complex)
        for a in range(d):
            X[(a + 1) % d, a] = 1.0
        self.add_op('X', X, hc='Xhc')
        self.add_op('Xhc', X.conj().T, hc='X')
        self.add_op('Xphc', X + X.conj().T, hc='Xphc')


# --- Custom Potts Model ---
class PottsModel(CouplingMPOModel):
    def init_sites(self, model_params):
        q = model_params.get('q', 5)
        return PottsSite(q)

    def init_terms(self, model_params):
        J = model_params.get('J', 1.0)
        g = model_params.get('g', 1.0)
        q = model_params.get('q', 5)
        # Bond: -J Σ_a P_a ⊗ P_a (Kronecker delta)
        for a in range(q):
            self.add_coupling(-J, 0, f'P{a}', 0, f'P{a}', 1)
        # Field: -g (X + X†)
        self.add_onsite(-g, 0, 'Xphc')


class PottsChain(PottsModel, NearestNeighborModel):
    pass


# --- Gell-Mann basis for d=5 ---
gm_mats = []
for j in range(q):
    for k in range(j+1, q):
        m = np.zeros((q, q), dtype=complex)
        m[j, k] = 1; m[k, j] = 1
        gm_mats.append((f'sym_{j}{k}', m.copy()))
for j in range(q):
    for k in range(j+1, q):
        m = np.zeros((q, q), dtype=complex)
        m[j, k] = -1j; m[k, j] = 1j
        gm_mats.append((f'asym_{j}{k}', m.copy()))
for l in range(1, q):
    m = np.zeros((q, q), dtype=complex)
    norm = np.sqrt(2.0 / (l * (l + 1)))
    for j in range(l):
        m[j, j] = norm
    m[l, l] = -l * norm
    gm_mats.append((f'diag_{l}', m.copy()))

op_names = [name for name, _ in gm_mats]
op_mats_all = [np.eye(q, dtype=complex)] + [mat for _, mat in gm_mats]


def entropy(rho):
    ev = la.eigvalsh(rho); ev = ev[ev > 1e-15]
    return float(-np.sum(ev * np.log2(ev)))


def mi_cv_potts(n, g_J, chi_max=15):
    model = PottsChain({'L': n, 'q': q, 'J': 1.0, 'g': g_J,
                         'bc_MPS': 'finite', 'conserve': None})
    site = model.lat.site(0)
    for name, mat in gm_mats:
        if name not in site.opnames:
            site.add_op(name, mat)

    psi = MPS.from_product_state(model.lat.mps_sites(), [0]*n, bc='finite')
    eng = dmrg.TwoSiteDMRGEngine(psi, model, {
        'mixer': True, 'max_E_err': 1e-6,
        'trunc_params': {'chi_max': chi_max}, 'max_sweeps': 6,
    })
    E0, _ = eng.run()

    exp_vals = {name: psi.expectation_value(name) for name in op_names}
    corr = {}
    for a in op_names:
        for b in op_names:
            corr[(a, b)] = psi.correlation_function(a, b)

    mi_vals = []
    for i in range(n):
        for j in range(i+1, n):
            rho_ij = np.eye(q*q, dtype=complex) / (q*q)
            for idx_a, a_name in enumerate(op_names):
                a_mat = op_mats_all[idx_a + 1]
                ev_i = complex(exp_vals[a_name][i])
                ev_j = complex(exp_vals[a_name][j])
                rho_ij += (ev_i / (2*q)) * np.kron(a_mat, np.eye(q))
                rho_ij += (ev_j / (2*q)) * np.kron(np.eye(q), a_mat)
                for idx_b, b_name in enumerate(op_names):
                    b_mat = op_mats_all[idx_b + 1]
                    c_val = complex(corr[(a_name, b_name)][i, j])
                    rho_ij += (c_val / 4.0) * np.kron(a_mat, b_mat)
            rho_ij = (rho_ij + rho_ij.conj().T) / 2
            ev = la.eigvalsh(rho_ij)
            if np.min(ev) < -1e-10:
                ev_pos = np.maximum(ev, 0)
                evec = la.eigh(rho_ij)[1]
                rho_ij = evec @ np.diag(ev_pos) @ evec.T
                rho_ij /= np.trace(rho_ij)
            rho_d = rho_ij.reshape(q, q, q, q)
            rho_i = np.trace(rho_d, axis1=1, axis2=3)
            rho_j = np.trace(rho_d, axis1=0, axis2=2)
            mi = entropy(rho_i) + entropy(rho_j) - entropy(rho_ij)
            mi_vals.append(max(float(np.real(mi)), 0))

    mi_pos = [m for m in mi_vals if m > 1e-10]
    if len(mi_pos) < 2: return 0.0, E0
    return float(np.std(mi_pos) / np.mean(mi_pos)), E0


print("=== Sprint 042a: True q=5 Potts Model — Timing Test ===\n")
sys.stdout.flush()

# Single test point: g=0.8 (should be in ordered phase)
t1 = time.time()
cv, E0 = mi_cv_potts(8, 0.8)
dt = time.time() - t1
print(f"  g/J=0.80: CV={cv:.4f}, E0={E0:.6f} ({dt:.1f}s)")

# Compare with clock result: q=5 clock at g=0.8 n=8 had CV=0.3616
print(f"  Clock q=5 at g=0.8 n=8: CV=0.3616 (from Sprint 041)")
print(f"  Potts q=5 at g=0.8 n=8: CV={cv:.4f}")
print(f"  Difference: {cv - 0.3616:.4f}")

results = {
    'experiment': '042a',
    'model': 'q=5 Potts (Kronecker delta)',
    'description': 'Timing test for custom Potts model with true delta coupling',
    'n': 8, 'q': 5, 'chi_max': 15,
    'data': {
        '0.8': {'cv': float(cv), 'E0': float(E0), 'time': float(dt)}
    },
    'comparison_clock': {'g_0.8_n8_cv': 0.3616},
    'total_runtime': time.time() - t0
}

with open('results/sprint_042a_potts_q5_timing.json', 'w') as f:
    json.dump(results, f, indent=2)

print(f"\nTotal runtime: {time.time()-t0:.1f}s")
print("Budget: ~16.5s/point for clock; Potts should be similar (same d=5)")
