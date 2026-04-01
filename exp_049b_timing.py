#!/usr/bin/env python3
"""Quick timing test for Potts DMRG at various q and n."""
import sys, time, numpy as np, warnings
warnings.filterwarnings('ignore')
from tenpy.models.tf_ising import TFIChain
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
            P = np.zeros((q,q),dtype=complex); P[a,a]=1.0
            self.add_op(f'P{a}', P)
        X = np.zeros((q,q),dtype=complex)
        for a in range(q): X[(a+1)%q,a]=1.0
        self.add_op('X',X,hc='Xhc'); self.add_op('Xhc',X.conj().T,hc='X')
        self.add_op('Xphc',X+X.conj().T,hc='Xphc')

class PottsChain(CouplingMPOModel, NearestNeighborModel):
    def init_sites(self, mp): return PottsSite(mp.get('q',3))
    def init_terms(self, mp):
        J,g,q = mp.get('J',1.0),mp.get('g',1.0),mp.get('q',3)
        for a in range(q): self.add_coupling(-J,0,f'P{a}',0,f'P{a}',1)
        self.add_onsite(-g,0,'Xphc')

tests = [
    (2, 32, 1.0, 40), (2, 48, 1.0, 40),
    (3, 12, 1.0, 30), (3, 16, 1.0, 30), (3, 24, 1.0, 30),
    (4, 12, 0.89, 30), (4, 16, 0.89, 30),
    (5, 12, 0.45, 30), (5, 16, 0.45, 30),
    (7, 8, 0.26, 30), (7, 12, 0.26, 30),
    (10, 8, 0.25, 40), (10, 12, 0.25, 40),
]

for q, n, g, chi in tests:
    t0 = time.time()
    if q == 2:
        model = TFIChain({'L':n,'J':1.0,'g':g,'bc_MPS':'finite'})
    else:
        model = PottsChain({'L':n,'q':q,'J':1.0,'g':g,'bc_MPS':'finite'})
    psi = MPS.from_product_state(model.lat.mps_sites(),[0]*n,bc='finite')
    eng = dmrg.TwoSiteDMRGEngine(psi,model,{'mixer':True,'max_E_err':1e-8,
        'trunc_params':{'chi_max':chi,'svd_min':1e-10},'max_sweeps':15})
    E0,_ = eng.run()
    S = float(psi.entanglement_entropy()[n//2-1])
    dt = time.time()-t0
    print(f'q={q:2d} n={n:3d} g={g:.2f} chi={chi}: S={S:.4f} E={E0:.6f} t={dt:.1f}s', flush=True)
    if dt > 55:
        print(f'  TOO SLOW — skip larger n for q={q}', flush=True)
