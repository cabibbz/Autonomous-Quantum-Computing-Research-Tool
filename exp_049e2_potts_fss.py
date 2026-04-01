#!/usr/bin/env python3
"""Sprint 049e2: q=3 Potts entropy FSS at the TRUE critical region g≈0.25-0.35."""
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
        for a in range(q): X[(a+1)%q, a] = 1.0
        self.add_op('X', X, hc='Xhc')
        self.add_op('Xhc', X.conj().T, hc='X')
        self.add_op('Xphc', X+X.conj().T, hc='Xphc')

class PottsChain(CouplingMPOModel, NearestNeighborModel):
    def init_sites(self, mp): return PottsSite(mp.get('q',3))
    def init_terms(self, mp):
        J,g,q = mp.get('J',1.0), mp.get('g',1.0), mp.get('q',3)
        for a in range(q): self.add_coupling(-J,0,f'P{a}',0,f'P{a}',1)
        self.add_onsite(-g,0,'Xphc')

def run(n, g, chi=40):
    t0 = time.time()
    m = PottsChain({'L':n,'q':3,'J':1.0,'g':g,'bc_MPS':'finite'})
    np.random.seed(42+n+int(g*1000))
    psi = MPS.from_product_state(m.lat.mps_sites(), [np.random.randint(3) for _ in range(n)], bc='finite')
    eng = dmrg.TwoSiteDMRGEngine(psi, m, {
        'mixer': True, 'max_E_err': 1e-10,
        'trunc_params': {'chi_max': chi, 'svd_min': 1e-12}, 'max_sweeps': 30})
    E, _ = eng.run()
    S = float(psi.entanglement_entropy()[n//2-1])
    return float(E), S, time.time()-t0

# Critical region only
g_vals = [0.15, 0.20, 0.24, 0.26, 0.28, 0.30, 0.32, 0.36, 0.40, 0.50]
results = {}

# n=8 already done, load from 049e
try:
    with open('results/sprint_049e_potts_critical.json') as f:
        prev = json.load(f)
    results[8] = prev['data'].get('8', [])
    print(f"Loaded n=8: {len(results[8])} points", flush=True)
except:
    pass

for n in [12, 16, 24]:
    print(f"\n=== q=3 n={n}, chi=40 ===", flush=True)
    results[n] = []
    t_start = time.time()
    for g in g_vals:
        if time.time() - t_start > 250:
            print(f"  Time limit", flush=True)
            break
        E, S, dt = run(n, g, chi=40)
        print(f"  g={g:.2f}: S={S:.6f}, E={E:.6f}, t={dt:.1f}s", flush=True)
        results[n].append({'g': float(g), 'E': float(E), 'S_half': float(S), 'time': dt})

    with open('results/sprint_049e_potts_critical.json', 'w') as f:
        json.dump({'sprint': '049e', 'q': 3, 'chi_max': 40,
                   'note': 'n=8 at chi=60, others chi=40',
                   'data': {str(k): v for k, v in results.items()}}, f, indent=2)

# Analysis
print("\n=== Entropy vs n at each g ===", flush=True)
for g in g_vals:
    line = f"g={g:.2f}: "
    for n in sorted(results.keys()):
        match = [d for d in results.get(n,[]) if abs(d['g']-g)<0.001]
        if match:
            line += f"S({n})={match[0]['S_half']:.4f} "
    print(f"  {line}", flush=True)

# Central charge
print("\n=== Central charge c from S vs ln(n) ===", flush=True)
for g in g_vals:
    ns, Ss = [], []
    for n in sorted(results.keys()):
        match = [d for d in results.get(n,[]) if abs(d['g']-g)<0.001]
        if match:
            ns.append(n); Ss.append(match[0]['S_half'])
    if len(ns) >= 3:
        c_fit = 6 * np.polyfit(np.log(ns), Ss, 1)[0]
        print(f"  g={g:.2f}: c = {c_fit:.3f} (from n={ns})", flush=True)

# dS/dg peak
print("\n=== dS/dg peak location ===", flush=True)
for n in sorted(results.keys()):
    data = sorted(results[n], key=lambda d: d['g'])
    gs = np.array([d['g'] for d in data])
    Ss = np.array([d['S_half'] for d in data])
    if len(gs) >= 4:
        dSdg = np.gradient(Ss, gs)
        peak = np.argmin(dSdg)  # Most negative = steepest drop
        print(f"  n={n}: steepest drop at g={gs[peak]:.3f}, dS/dg={dSdg[peak]:.3f}", flush=True)

print("\nDone!", flush=True)
