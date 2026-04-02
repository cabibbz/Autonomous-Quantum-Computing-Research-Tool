#!/usr/bin/env python3
"""Sprint 068b part 2: 2D gap*L for q=5 — THE KEY TEST.

q=5 is the first q where 2D classical Potts is first-order.
If our hybrid model remains continuous in 2D at q=5, that's remarkable.

L=2 (dim=625): full scan, 40 g points
L=3 (dim=1953125): coarse scan, 10 g points (~22s each = 220s total)
"""
import numpy as np, json, time
from scipy.sparse import csr_matrix, kron as sp_kron, eye as sp_eye, diags
from scipy.sparse.linalg import eigsh

def build_H_2d(Lx, Ly, q, g):
    n = Lx*Ly; dim = q**n
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
    for gi, g in enumerate(g_values):
        t0 = time.time()
        H = build_H_2d(L,L,q,g)
        if dim < 500: evals = np.linalg.eigvalsh(H.toarray())[:6]
        else: evals = np.sort(eigsh(H, k=min(6, dim-2), which='SA')[0])
        gap = evals[1]-evals[0]
        gap2 = evals[2]-evals[0] if len(evals)>2 else None
        dt = time.time()-t0
        data.append({'g':round(g,5),'gap':round(float(gap),10),'gap_x_L':round(float(gap*L),8),
                      'gap2':round(float(gap2),10) if gap2 else None,
                      'E0_per_site':round(float(evals[0]/n),8),'time_s':round(dt,1)})
        print(f"    g={g:.3f}: gap={gap:.6e}, gap*L={gap*L:.4f} ({dt:.1f}s)")
    return data

results = {'experiment':'068b_q5','timestamp':time.strftime('%Y-%m-%d %H:%M:%S')}

# L=2 (dim=625)
print("q=5, L=2 (dim=625)")
g_scan_small = np.linspace(0.3, 3.0, 40)
t0=time.time()
data_L2 = scan(2, 5, g_scan_small)
print(f"  Total: {time.time()-t0:.1f}s")
results['L2'] = data_L2

# Find gap minimum at L=2
gaps_L2 = [d['gap'] for d in data_L2]
min_idx = np.argmin(gaps_L2)
print(f"  Gap min: {gaps_L2[min_idx]:.6e} at g={data_L2[min_idx]['g']:.3f}")

# L=3 (dim=1953125) — coarse scan near expected g_c
# From L=2 data, identify transition region
# In 1D, g_c(q=5) = 0.441. 2D ratio ~2-3x based on q=2,3 results.
# Expect g_c(2D) ~ 1.0-1.5 for q=5
print(f"\nq=5, L=3 (dim=1953125)")
g_scan_large = np.linspace(0.6, 2.0, 10)  # Coarse, 10 points
t0=time.time()
data_L3 = scan(3, 5, g_scan_large)
dt = time.time()-t0
print(f"  Total: {dt:.1f}s ({dt/10:.1f}s per point)")
results['L3'] = data_L3

# Find crossing
from scipy.interpolate import interp1d
g2=[d['g'] for d in data_L2]; y2=[d['gap_x_L'] for d in data_L2]
g3=[d['g'] for d in data_L3]; y3=[d['gap_x_L'] for d in data_L3]
gm=max(min(g2),min(g3)); gx=min(max(g2),max(g3))
gg=np.linspace(gm,gx,500)
f2=interp1d(g2,y2,kind='cubic',fill_value='extrapolate')
f3=interp1d(g3,y3,kind='cubic',fill_value='extrapolate')
diff=f2(gg)-f3(gg)
crossings=[]
for i in range(len(diff)-1):
    if diff[i]*diff[i+1]<0:
        gc=gg[i]-diff[i]*(gg[i+1]-gg[i])/(diff[i+1]-diff[i])
        crossings.append(round(float(gc),4))

print(f"\ngap*L crossing (q=5, L=2 vs L=3): {crossings}")
results['crossing_2_3'] = crossings

# First-order test: gap minimum at L=3
gaps_L3 = [d['gap'] for d in data_L3]
min_idx3 = np.argmin(gaps_L3)
print(f"Gap min L=2: {min(gaps_L2):.6e}")
print(f"Gap min L=3: {min(gaps_L3):.6e}")
if min(gaps_L2) > 0 and min(gaps_L3) > 0:
    ratio = min(gaps_L3) / min(gaps_L2)
    # First-order: gap ~ exp(-α·L²), ratio ~ exp(-α·(9-4)) = exp(-5α)
    # Continuous: gap ~ L^{-z}, ratio ~ (3/2)^{-z}
    print(f"Ratio gap_min(L=3)/gap_min(L=2) = {ratio:.4f}")
    print(f"  If continuous (z=1): expect (2/3)^1 = 0.667")
    print(f"  If first-order: expect exp(-5α) << 0.667")

# Also: degeneracy pattern at transition
print("\nDegeneracy check near transition:")
for data, L in [(data_L2, 2), (data_L3, 3)]:
    for d in data:
        if d['gap2'] is not None and d['gap'] > 1e-10:
            r = d['gap2'] / d['gap']
            if abs(d['g'] - (crossings[0] if crossings else 1.0)) < 0.3:
                print(f"  L={L}, g={d['g']:.3f}: gap2/gap1 = {r:.3f}")

with open("results/sprint_068b_q5.json","w") as f:
    json.dump(results,f,indent=2,default=str)
print("\nSaved results/sprint_068b_q5.json")
