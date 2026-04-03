#!/usr/bin/env python3
"""Sprint 067b v2: Targeted entropy scan at q=5 using DMRG.

Only 6 key g values at n=16,24 to detect floating phase.
Floating phase would show c_eff=1 in a range of g ABOVE g_c.
Use chi=20 for speed (sufficient for c extraction at n≤24).
"""
import numpy as np, json, time, sys

import tenpy
from tenpy.models.model import CouplingMPOModel
from tenpy.networks.site import Site
from tenpy.algorithms import dmrg
from tenpy.networks.mps import MPS
import tenpy.linalg.np_conserved as npc

class ZqSite(Site):
    def __init__(self, q):
        ops = {}
        X = np.zeros((q, q), dtype=complex)
        for s in range(q):
            X[(s+1) % q, s] = 1.0
        ops['XpXd'] = (X + X.conj().T).real.astype(float)
        for a in range(q):
            Pa = np.zeros((q, q))
            Pa[a, a] = 1.0
            ops[f'P{a}'] = Pa
        leg = npc.LegCharge.from_trivial(q)
        super().__init__(leg, ['s' + str(i) for i in range(q)], **ops)

class HybridPottsModel(CouplingMPOModel):
    def init_sites(self, model_params):
        return ZqSite(model_params.get('q', 5))
    def init_terms(self, model_params):
        q = model_params.get('q', 5)
        g = model_params.get('g', 1.0)
        for a in range(q):
            self.add_coupling(-1.0, 0, f'P{a}', 0, f'P{a}', 1)
        self.add_onsite(-g, 0, 'XpXd')

def run_dmrg_entropy(q, n, g, chi=20):
    """Run DMRG and return half-chain entropy."""
    model_params = {'q': q, 'g': float(g), 'L': n, 'bc_MPS': 'finite'}
    M = HybridPottsModel(model_params)
    psi = MPS.from_lat_product_state(M.lat, [['s0']])
    dmrg_params = {
        'mixer': True,
        'max_E_err': 1e-8,
        'trunc_params': {'chi_max': chi, 'svd_min': 1e-10},
        'max_sweeps': 15,
    }
    eng = dmrg.TwoSiteDMRGEngine(psi, M, dmrg_params)
    E0, _ = eng.run()
    S_half = psi.entanglement_entropy()[n // 2 - 1]
    return float(E0), float(S_half)

# Key g values: ordered, sub-critical, critical, super-critical, disordered, far-disordered
gc = 0.441
g_targets = [0.20, 0.35, 0.44, 0.55, 0.80, 1.20]

results = {"experiment": "067b_v2_entropy_scan", "q": 5, "gc": gc, "g_values": g_targets}

for n in [16, 24]:
    print(f"\nDMRG n={n}, chi=20", flush=True)
    data = []
    for g in g_targets:
        t0 = time.time()
        try:
            E0, S = run_dmrg_entropy(5, n, g, chi=20)
            dt = time.time() - t0
            print(f"  g={g:.3f}: S_half={S:.4f}, E0={E0:.6f} ({dt:.1f}s)", flush=True)
            data.append({"g": float(g), "S_half": S, "E0": E0, "time_s": dt})
        except Exception as e:
            dt = time.time() - t0
            print(f"  g={g:.3f}: FAILED ({e}) ({dt:.1f}s)", flush=True)
            data.append({"g": float(g), "S_half": None, "E0": None, "time_s": dt, "error": str(e)})
    results[f"n{n}"] = data

# Compute effective c from n=16,24 pairs
print("\n" + "=" * 60)
print("Effective central charge from DMRG n=16 vs n=24")
print("=" * 60)
# For open BC half chain: S(n/2) ≈ (c/6)*ln[(2n/π)*sin(π/2)] + const = (c/6)*ln(2n/π) + const
# c_eff = 6*(S24 - S16) / [ln(2*24/π) - ln(2*16/π)] = 6*(S24-S16)/[ln(24)-ln(16)]
print(f"{'g':>8} {'S(16)':>10} {'S(24)':>10} {'c_eff':>10} {'interpretation':>20}")
print("-" * 62)

c_eff_data = []
for i, g in enumerate(g_targets):
    s16 = results["n16"][i]["S_half"]
    s24 = results["n24"][i]["S_half"]
    if s16 is not None and s24 is not None:
        c_eff = 6 * (s24 - s16) / (np.log(24) - np.log(16))
        interp = ""
        if c_eff < 0.1:
            interp = "gapped"
        elif 0.4 < c_eff < 0.6:
            interp = "c≈0.5 (Ising?)"
        elif 0.85 < c_eff < 1.15:
            interp = "c≈1 (FLOAT?)"
        elif c_eff > 1.0:
            interp = f"c>{c_eff:.1f} (critical)"
        else:
            interp = f"intermediate"
        c_eff_data.append({"g": float(g), "c_eff": float(c_eff)})
        print(f"{g:8.3f} {s16:10.4f} {s24:10.4f} {c_eff:10.4f} {interp:>20}")
    else:
        c_eff_data.append({"g": float(g), "c_eff": None})
        print(f"{g:8.3f} {'FAIL':>10} {'FAIL':>10} {'—':>10}")

results["c_eff"] = c_eff_data

# Also include the exact diag c_eff we already computed
# Re-compute from n=6,8 for same g values
from scipy.sparse import csr_matrix, kron as sp_kron, eye as sp_eye, diags
from scipy.sparse.linalg import eigsh
from scipy.linalg import svd

def build_H_open(n, q, g):
    dim = q**n
    potts_2site = np.zeros(q**2)
    for a in range(q): potts_2site[a*q + a] = 1.0
    potts_op = diags(potts_2site, 0, shape=(q**2, q**2), format='csr')
    X = np.zeros((q, q))
    for s in range(q): X[(s+1) % q, s] = 1.0
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

def half_entropy(psi, n, q):
    nA = n // 2; dimA = q**nA; dimB = q**(n - nA)
    s = svd(psi.reshape(dimA, dimB), compute_uv=False)
    s = s[s > 1e-14]
    return -np.sum(s**2 * np.log(s**2))

print("\n\nExact diag comparison (n=6,8):")
print(f"{'g':>8} {'c_eff(6,8)':>12} {'c_eff(16,24)':>14}")
for g in g_targets:
    H6 = build_H_open(6, 5, g); H8 = build_H_open(8, 5, g)
    _, v6 = eigsh(H6, k=1, which='SA'); _, v8 = eigsh(H8, k=1, which='SA')
    S6 = half_entropy(v6[:,0], 6, 5); S8 = half_entropy(v8[:,0], 8, 5)
    c68 = 6*(S8-S6)/(np.log(8)-np.log(6))
    # Find DMRG c_eff
    cd = next((d["c_eff"] for d in c_eff_data if abs(d["g"]-g)<0.001), None)
    cd_str = f"{cd:.4f}" if cd is not None else "—"
    print(f"{g:8.3f} {c68:12.4f} {cd_str:>14}")

with open("results/sprint_067b_entropy_scan.json", "w") as f:
    json.dump(results, f, indent=2, default=str)
print("\nSaved to results/sprint_067b_entropy_scan.json", flush=True)

from db_utils import record
for d in c_eff_data:
    if d["c_eff"] is not None:
        record(sprint=67, model='hybrid', q=5, n=24, quantity='c_eff',
               value=d["c_eff"], method='dmrg_16_24',
               notes=f'g={d["g"]:.3f}')
print("Recorded to DB.", flush=True)
