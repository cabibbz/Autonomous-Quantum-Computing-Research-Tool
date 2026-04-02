#!/usr/bin/env python3
"""Sprint 068b part 1: 2D gap*L for q=2,3 only (fast)."""
import numpy as np, json, time
from scipy.sparse import csr_matrix, kron as sp_kron, eye as sp_eye, diags
from scipy.sparse.linalg import eigsh

def build_H_2d(Lx, Ly, q, g):
    n = Lx * Ly; dim = q**n
    potts_2site = np.zeros(q**2)
    for a in range(q): potts_2site[a*q+a] = 1.0
    potts_op = diags(potts_2site, 0, shape=(q**2, q**2), format='csr')
    X = np.zeros((q,q))
    for s in range(q): X[(s+1)%q, s] = 1.0
    XpXd = csr_matrix(X + X.T)
    H = csr_matrix((dim, dim))
    bonds = []
    for y in range(Ly):
        for x in range(Lx):
            site = y*Lx+x
            nbr_h = y*Lx+(x+1)%Lx
            if nbr_h != site: bonds.append((min(site,nbr_h), max(site,nbr_h)))
            nbr_v = ((y+1)%Ly)*Lx+x
            if nbr_v != site: bonds.append((min(site,nbr_v), max(site,nbr_v)))
    bonds = list(set(bonds))
    for (i,j) in bonds:
        if j == i+1:
            left = q**i; right = q**(n-i-2); op = potts_op
            if left > 1: op = sp_kron(sp_eye(left), op, format='csr')
            if right > 1: op = sp_kron(op, sp_eye(right), format='csr')
        else:
            dv = np.zeros(dim)
            for idx in range(dim):
                si = (idx//(q**i))%q; sj = (idx//(q**j))%q
                dv[idx] = 1.0 if si == sj else 0.0
            op = diags(dv, 0, shape=(dim,dim), format='csr')
        H = H - op
    for i in range(n):
        left = q**i; right = q**(n-i-1); op = XpXd.copy()
        if left > 1: op = sp_kron(sp_eye(left), op, format='csr')
        if right > 1: op = sp_kron(op, sp_eye(right), format='csr')
        H = H - g * op
    return H

def scan(L, q, g_values):
    n = L*L; dim = q**n; data = []
    for g in g_values:
        H = build_H_2d(L,L,q,g)
        if dim < 500: evals = np.linalg.eigvalsh(H.toarray())[:4]
        else: evals = np.sort(eigsh(H, k=4, which='SA')[0])
        gap = evals[1]-evals[0]
        data.append({'g':round(g,5),'gap':round(float(gap),10),'gap_x_L':round(float(gap*L),8),
                      'E0_per_site':round(float(evals[0]/n),8)})
    return data

def find_crossing(d1, d2):
    from scipy.interpolate import interp1d
    g1=[d['g'] for d in d1]; y1=[d['gap_x_L'] for d in d1]
    g2=[d['g'] for d in d2]; y2=[d['gap_x_L'] for d in d2]
    gm=max(min(g1),min(g2)); gx=min(max(g1),max(g2))
    gg=np.linspace(gm,gx,500)
    f1=interp1d(g1,y1,kind='cubic',fill_value='extrapolate')
    f2=interp1d(g2,y2,kind='cubic',fill_value='extrapolate')
    diff=f1(gg)-f2(gg); cr=[]
    for i in range(len(diff)-1):
        if diff[i]*diff[i+1]<0:
            gc=gg[i]-diff[i]*(gg[i+1]-gg[i])/(diff[i+1]-diff[i])
            cr.append(round(float(gc),4))
    return cr

results = {'experiment':'068b_q2q3','timestamp':time.strftime('%Y-%m-%d %H:%M:%S')}

# q=2: L=2,3,4
print("q=2: L=2,3,4")
g_scan = np.linspace(0.2, 1.5, 40)
q2 = {}
for L in [2,3,4]:
    t0=time.time(); q2[L]=scan(L,2,g_scan)
    print(f"  L={L} (dim={2**(L*L)}): {time.time()-t0:.1f}s")

for L1,L2 in [(2,3),(3,4),(2,4)]:
    cr=find_crossing(q2[L1],q2[L2])
    print(f"  Crossing L={L1} vs {L2}: {cr}")
    results[f'q2_cross_{L1}_{L2}']=cr
results['q2']={str(L):v for L,v in q2.items()}

# q=3: L=2,3
print("\nq=3: L=2,3")
g_scan3 = np.linspace(0.3, 2.0, 40)
q3 = {}
for L in [2,3]:
    t0=time.time(); q3[L]=scan(L,3,g_scan3)
    print(f"  L={L} (dim={3**(L*L)}): {time.time()-t0:.1f}s")
cr=find_crossing(q3[2],q3[3])
print(f"  Crossing L=2 vs 3: {cr}")
results['q3_cross_2_3']=cr
results['q3']={str(L):v for L,v in q3.items()}

# Print gap*L at crossing for convergence check
print("\nGap*L at crossings:")
from scipy.interpolate import interp1d
for q_val,qdata,key in [(2,q2,'q2_cross_3_4'),(3,q3,'q3_cross_2_3')]:
    cr = results.get(key,[])
    if cr:
        gc = cr[0]
        for L,data in sorted(qdata.items()):
            g_arr=[d['g'] for d in data]; gL_arr=[d['gap_x_L'] for d in data]
            f=interp1d(g_arr,gL_arr,kind='cubic',fill_value='extrapolate')
            print(f"  q={q_val}, L={L}: gap*L(g_c={gc:.3f}) = {float(f(gc)):.4f}")

with open("results/sprint_068b_q2q3.json","w") as f:
    json.dump(results,f,indent=2,default=str)
print("\nSaved results/sprint_068b_q2q3.json")
