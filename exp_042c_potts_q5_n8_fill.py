"""
Sprint 042c: True q=5 Potts n=8 — fill transition region g=0.35-0.45.
CV jumps from 0.16 (g=0.3) to 1.06 (g=0.5). Need finer resolution.
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

class PottsSite(Site):
    def __init__(self, q):
        leg = npc.LegCharge.from_trivial(q)
        Site.__init__(self, leg, [str(a) for a in range(q)], sort_charge=False)
        for a in range(q):
            P = np.zeros((q,q), dtype=complex); P[a,a]=1.0
            self.add_op(f'P{a}', P)
        X = np.zeros((q,q), dtype=complex)
        for a in range(q): X[(a+1)%q, a]=1.0
        self.add_op('X', X, hc='Xhc')
        self.add_op('Xhc', X.conj().T, hc='X')
        self.add_op('Xphc', X+X.conj().T, hc='Xphc')

class PottsModel(CouplingMPOModel):
    def init_sites(self, model_params):
        return PottsSite(model_params.get('q', 5))
    def init_terms(self, model_params):
        J, g, qq = model_params.get('J',1.), model_params.get('g',1.), model_params.get('q',5)
        for a in range(qq): self.add_coupling(-J, 0, f'P{a}', 0, f'P{a}', 1)
        self.add_onsite(-g, 0, 'Xphc')

class PottsChain(PottsModel, NearestNeighborModel):
    pass

gm_mats = []
for j in range(q):
    for k in range(j+1,q):
        m=np.zeros((q,q),dtype=complex); m[j,k]=1; m[k,j]=1
        gm_mats.append((f'sym_{j}{k}', m))
for j in range(q):
    for k in range(j+1,q):
        m=np.zeros((q,q),dtype=complex); m[j,k]=-1j; m[k,j]=1j
        gm_mats.append((f'asym_{j}{k}', m))
for l in range(1,q):
    m=np.zeros((q,q),dtype=complex); norm=np.sqrt(2./(l*(l+1)))
    for j in range(l): m[j,j]=norm
    m[l,l]=-l*norm
    gm_mats.append((f'diag_{l}', m))
op_names = [n for n,_ in gm_mats]
op_mats_all = [np.eye(q,dtype=complex)] + [m for _,m in gm_mats]

def entropy(rho):
    ev=la.eigvalsh(rho); ev=ev[ev>1e-15]; return float(-np.sum(ev*np.log2(ev)))

def mi_cv_potts(n, g_J, chi_max=15):
    model = PottsChain({'L':n,'q':q,'J':1.0,'g':g_J,'bc_MPS':'finite','conserve':None})
    site = model.lat.site(0)
    for name, mat in gm_mats:
        if name not in site.opnames: site.add_op(name, mat)
    psi = MPS.from_product_state(model.lat.mps_sites(), [0]*n, bc='finite')
    eng = dmrg.TwoSiteDMRGEngine(psi, model, {
        'mixer':True,'max_E_err':1e-6,'trunc_params':{'chi_max':chi_max},'max_sweeps':6})
    E0,_ = eng.run()
    exp_vals = {name: psi.expectation_value(name) for name in op_names}
    corr = {(a,b): psi.correlation_function(a,b) for a in op_names for b in op_names}
    mi_vals = []
    for i in range(n):
        for j in range(i+1,n):
            rho_ij=np.eye(q*q,dtype=complex)/(q*q)
            for ia,an in enumerate(op_names):
                am=op_mats_all[ia+1]; ei=complex(exp_vals[an][i]); ej=complex(exp_vals[an][j])
                rho_ij+=(ei/(2*q))*np.kron(am,np.eye(q))+(ej/(2*q))*np.kron(np.eye(q),am)
                for ib,bn in enumerate(op_names):
                    bm=op_mats_all[ib+1]; cv_=complex(corr[(an,bn)][i,j])
                    rho_ij+=(cv_/4.)*np.kron(am,bm)
            rho_ij=(rho_ij+rho_ij.conj().T)/2
            ev=la.eigvalsh(rho_ij)
            if np.min(ev)<-1e-10:
                ev_pos=np.maximum(ev,0); evec=la.eigh(rho_ij)[1]
                rho_ij=evec@np.diag(ev_pos)@evec.T; rho_ij/=np.trace(rho_ij)
            rd=rho_ij.reshape(q,q,q,q)
            ri=np.trace(rd,axis1=1,axis2=3); rj=np.trace(rd,axis1=0,axis2=2)
            mi=entropy(ri)+entropy(rj)-entropy(rho_ij)
            mi_vals.append(max(float(np.real(mi)),0))
    mi_pos=[m for m in mi_vals if m>1e-10]
    if len(mi_pos)<2: return 0.0
    return float(np.std(mi_pos)/np.mean(mi_pos))

print("=== Sprint 042c: q=5 Potts n=8 — Transition Fill ===\n")
sys.stdout.flush()

results = {'experiment':'042c','model':'q=5 Potts','n':8,'chi_max':15,'data':{}}
for g_J in [0.35, 0.40, 0.45]:
    elapsed=time.time()-t0
    if elapsed>50: print(f"  TIMEOUT {elapsed:.0f}s"); break
    t1=time.time()
    try:
        cv=mi_cv_potts(8, g_J)
        dt=time.time()-t1
        print(f"  g/J={g_J:.2f}: CV={cv:.4f} ({dt:.1f}s)")
        results['data'][str(g_J)]={'cv':float(cv),'time':float(dt)}
    except Exception as e:
        print(f"  g/J={g_J:.2f}: FAILED: {e}")
        results['data'][str(g_J)]={'cv':None,'error':str(e)}
    sys.stdout.flush()
    results['total_runtime']=time.time()-t0
    with open('results/sprint_042c_potts_q5_n8_fill.json','w') as f:
        json.dump(results, f, indent=2)

print(f"\nTotal: {time.time()-t0:.1f}s")
print("\nFull n=8 picture so far:")
print("  g=0.10: 0.0117  g=0.30: 0.1567  g=0.50: 1.0645  g=0.80: 1.4922")
