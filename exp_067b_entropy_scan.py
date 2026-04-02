#!/usr/bin/env python3
"""Sprint 067b: Entanglement entropy scan across full phase diagram.

Use DMRG (TeNPy) at q=5 to measure S(n/2) at n=16,24 across g=[0.1, 1.5].
A floating phase would show: S ~ (c_eff/3)*ln(n) with c_eff=1 (free boson).
Compare effective c at each g from the two-size entropy difference.

Also exact diag at n=8 for q=5,7 for wider coverage.
"""
import numpy as np, json, time, sys

# ============================================================
# Part 1: Exact diag entropy at n=8 for q=5
# ============================================================
from scipy.sparse import csr_matrix, kron as sp_kron, eye as sp_eye, diags
from scipy.sparse.linalg import eigsh
from scipy.linalg import svd

def build_H_open(n, q, g):
    """Open BC Hamiltonian for entropy measurement."""
    dim = q**n
    potts_2site = np.zeros(q**2)
    for a in range(q):
        potts_2site[a*q + a] = 1.0
    potts_op = diags(potts_2site, 0, shape=(q**2, q**2), format='csr')
    X = np.zeros((q, q))
    for s in range(q):
        X[(s+1) % q, s] = 1.0
    XpXd = csr_matrix(X + X.T)

    H = csr_matrix((dim, dim))
    for i in range(n - 1):
        left = q**i; right = q**(n - i - 2)
        op = potts_op
        if left > 1: op = sp_kron(sp_eye(left), op, format='csr')
        if right > 1: op = sp_kron(op, sp_eye(right), format='csr')
        H = H - op
    for i in range(n):
        left = q**i; right = q**(n - i - 1)
        op = XpXd
        if left > 1: op = sp_kron(sp_eye(left), op, format='csr')
        if right > 1: op = sp_kron(op, sp_eye(right), format='csr')
        H = H - g * op
    return H

def half_chain_entropy(psi, n, q):
    """Entanglement entropy of left half."""
    nA = n // 2
    dimA = q**nA
    dimB = q**(n - nA)
    psi_mat = psi.reshape(dimA, dimB)
    s_vals = svd(psi_mat, compute_uv=False)
    s_vals = s_vals[s_vals > 1e-14]
    return -np.sum(s_vals**2 * np.log(s_vals**2))

print("=" * 60)
print("Part 1: Exact diag entropy scan at q=5, n=8")
print("=" * 60)

q = 5; n = 8
g_values_ed = np.concatenate([
    np.linspace(0.10, 0.35, 6),
    np.linspace(0.38, 0.55, 10),   # Dense near g_c=0.441
    np.linspace(0.60, 1.50, 10),
])
g_values_ed = np.sort(np.unique(g_values_ed))

ed_data = []
print(f"q={q}, n={n}, {len(g_values_ed)} g-points", flush=True)
t0 = time.time()
for g in g_values_ed:
    H = build_H_open(n, q, g)
    evals, evecs = eigsh(H, k=1, which='SA')
    psi = evecs[:, 0]
    S = half_chain_entropy(psi, n, q)
    ed_data.append({"g": float(g), "S_half": float(S), "E0": float(evals[0])})
    print(f"  g={g:.3f}: S={S:.4f}", flush=True)
dt_ed = time.time() - t0
print(f"Done in {dt_ed:.1f}s", flush=True)

# Also n=6 for effective c extraction
print(f"\nq={q}, n=6 for c extraction", flush=True)
ed_data_n6 = []
t0 = time.time()
for g in g_values_ed:
    H = build_H_open(6, q, g)
    evals, evecs = eigsh(H, k=1, which='SA')
    psi = evecs[:, 0]
    S = half_chain_entropy(psi, 6, q)
    ed_data_n6.append({"g": float(g), "S_half": float(S), "E0": float(evals[0])})
dt_ed6 = time.time() - t0
print(f"Done in {dt_ed6:.1f}s", flush=True)

# Effective c from two sizes: S(n/2) ~ (c/6)*ln(n) + const
# c_eff = 6*(S8 - S6) / (ln(8) - ln(6))
print(f"\nEffective central charge from n=6,8 exact diag:")
print(f"{'g':>8} {'S(n=6)':>10} {'S(n=8)':>10} {'c_eff':>10}")
print("-" * 42)
c_eff_ed = []
for i in range(len(g_values_ed)):
    S6 = ed_data_n6[i]["S_half"]
    S8 = ed_data[i]["S_half"]
    g = g_values_ed[i]
    # c = 6 * (S8 - S6) / (ln(8) - ln(6)) for open BC half-chain
    # More precisely: S(l, n) = (c/6)*ln[(2n/pi)*sin(pi*l/n)] + const
    # For half chain: l=n/2, so chord = (2n/pi)*sin(pi/2) = 2n/pi
    # S(n/2) = (c/6)*ln(2n/pi) + const = (c/6)*ln(n) + (c/6)*ln(2/pi) + const
    # c_eff = 6*(S8 - S6) / (ln(8) - ln(6))
    c_eff = 6 * (S8 - S6) / (np.log(8) - np.log(6))
    c_eff_ed.append(float(c_eff))
    print(f"{g:8.3f} {S6:10.4f} {S8:10.4f} {c_eff:10.4f}")

# ============================================================
# Part 2: DMRG entropy scan at q=5, n=16,24
# ============================================================
print("\n" + "=" * 60)
print("Part 2: DMRG entropy scan at q=5")
print("=" * 60)

import tenpy
from tenpy.models.model import CouplingMPOModel
from tenpy.networks.site import SpinSite, Site
from tenpy.models.lattice import Chain
from tenpy.algorithms import dmrg
from tenpy.networks.mps import MPS
import tenpy.linalg.np_conserved as npc

class ZqSite(Site):
    """Z_q site for clock/Potts models."""
    def __init__(self, q):
        d = q
        # Operators
        # Identity
        ops = {}
        # Z_q clock operator: Z|s> = omega^s |s>
        omega = np.exp(2j * np.pi / q)
        Z = np.diag([omega**s for s in range(q)])
        ops['Z'] = Z
        ops['Zd'] = Z.conj()
        # X shift operator: X|s> = |s+1 mod q>
        X = np.zeros((q, q), dtype=complex)
        for s in range(q):
            X[(s+1) % q, s] = 1.0
        ops['X'] = X
        ops['Xd'] = X.conj().T
        ops['XpXd'] = (X + X.conj().T).real.astype(float)  # Hermitian
        # Projectors for Potts coupling: sum_a |a><a| x |a><a|
        for a in range(q):
            Pa = np.zeros((q, q))
            Pa[a, a] = 1.0
            ops[f'P{a}'] = Pa
        leg = npc.LegCharge.from_trivial(d)
        super().__init__(leg, ['s' + str(i) for i in range(d)], **ops)

class HybridPottsModel(CouplingMPOModel):
    def init_sites(self, model_params):
        q = model_params.get('q', 5)
        return ZqSite(q)

    def init_terms(self, model_params):
        q = model_params.get('q', 5)
        g = model_params.get('g', 1.0)
        # Potts coupling: -sum_a P_a(i) P_a(i+1)
        for a in range(q):
            self.add_coupling(-1.0, 0, f'P{a}', 0, f'P{a}', 1)
        # Transverse field: -g*(X + X†)
        self.add_onsite(-g, 0, 'XpXd')

# Scan g values for DMRG
g_values_dmrg = np.concatenate([
    np.linspace(0.15, 0.35, 5),
    np.linspace(0.38, 0.55, 8),   # Dense near g_c
    np.linspace(0.60, 1.20, 7),
])
g_values_dmrg = np.sort(np.unique(g_values_dmrg))

dmrg_data = {}
for n_dmrg in [16, 24]:
    print(f"\nDMRG n={n_dmrg}, chi=30, {len(g_values_dmrg)} g-points", flush=True)
    dmrg_data[n_dmrg] = []
    for g in g_values_dmrg:
        t0 = time.time()
        model_params = {
            'q': 5, 'g': float(g),
            'L': n_dmrg,
            'bc_MPS': 'finite',
            'conserve': None,
        }
        try:
            M = HybridPottsModel(model_params)
            psi = MPS.from_lat_product_state(M.lat, [['s0']])
            dmrg_params = {
                'mixer': True,
                'max_E_err': 1e-8,
                'trunc_params': {'chi_max': 30, 'svd_min': 1e-10},
                'max_sweeps': 20,
            }
            eng = dmrg.TwoSiteDMRGEngine(psi, M, dmrg_params)
            E0, _ = eng.run()
            # Half-chain entropy
            S_half = psi.entanglement_entropy()[n_dmrg // 2 - 1]
            dt = time.time() - t0
            dmrg_data[n_dmrg].append({"g": float(g), "S_half": float(S_half),
                                       "E0": float(E0), "time_s": float(dt)})
            print(f"  g={g:.3f}: S={S_half:.4f} ({dt:.1f}s)", flush=True)
        except Exception as e:
            dt = time.time() - t0
            print(f"  g={g:.3f}: FAILED ({e}) ({dt:.1f}s)", flush=True)
            dmrg_data[n_dmrg].append({"g": float(g), "S_half": None, "E0": None,
                                       "time_s": float(dt), "error": str(e)})

# Effective c from DMRG n=16,24
print(f"\nEffective c from DMRG n=16,24:")
print(f"{'g':>8} {'S(n=16)':>10} {'S(n=24)':>10} {'c_eff':>10}")
print("-" * 42)
c_eff_dmrg = []
for i in range(len(g_values_dmrg)):
    d16 = dmrg_data[16][i]
    d24 = dmrg_data[24][i]
    g = g_values_dmrg[i]
    if d16["S_half"] is not None and d24["S_half"] is not None:
        c_eff = 6 * (d24["S_half"] - d16["S_half"]) / (np.log(24) - np.log(16))
        c_eff_dmrg.append({"g": float(g), "c_eff": float(c_eff)})
        print(f"{g:8.3f} {d16['S_half']:10.4f} {d24['S_half']:10.4f} {c_eff:10.4f}")
    else:
        c_eff_dmrg.append({"g": float(g), "c_eff": None})
        print(f"{g:8.3f} {'FAIL':>10} {'FAIL':>10} {'—':>10}")

# ============================================================
# Save results
# ============================================================
results = {
    "experiment": "067b_entropy_scan",
    "ed_q5_n8": ed_data,
    "ed_q5_n6": ed_data_n6,
    "c_eff_ed": [{"g": float(g_values_ed[i]), "c_eff": c_eff_ed[i]} for i in range(len(g_values_ed))],
    "dmrg_q5_n16": dmrg_data.get(16, []),
    "dmrg_q5_n24": dmrg_data.get(24, []),
    "c_eff_dmrg": c_eff_dmrg,
}

with open("results/sprint_067b_entropy_scan.json", "w") as f:
    json.dump(results, f, indent=2, default=str)
print("\nSaved to results/sprint_067b_entropy_scan.json", flush=True)

# Record to DB
from db_utils import record
# Record peak c_eff and its location
for dataset, label in [(c_eff_dmrg, "dmrg_16_24")]:
    valid = [(d["g"], d["c_eff"]) for d in dataset if d["c_eff"] is not None]
    if valid:
        g_peak, c_peak = max(valid, key=lambda x: x[1])
        record(sprint=67, model='hybrid', q=5, n=24, quantity='c_eff_peak',
               value=c_peak, method=label, notes=f'g_peak={g_peak:.3f}')
        record(sprint=67, model='hybrid', q=5, n=24, quantity='c_eff_peak_g',
               value=g_peak, method=label)

print("Recorded to DB.", flush=True)
