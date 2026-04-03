#!/usr/bin/env python3
"""Sprint 071c: Exact diag energy gap on Ly=2 cylinder for q=2,5.

Vectorized bond construction for speed.
q=2: Lx=4,5,6,7 (dim up to 16k)
q=5: Lx=3,4 (dim up to 390k). Lx=5 (10M) only if fast enough.
"""
import numpy as np
import json
import time
from scipy.sparse import csr_matrix, kron as sp_kron, eye as sp_eye, diags
from gpu_utils import eigsh


def build_H_cylinder(Lx, Ly, q, g):
    """Build hybrid Hamiltonian on Lx x Ly cylinder. Vectorized."""
    n = Lx * Ly
    dim = q ** n

    # Adjacent Potts coupling operator (for j = i+1)
    potts_2site = np.zeros(q ** 2)
    for a in range(q):
        potts_2site[a * q + a] = 1.0
    potts_op = diags(potts_2site, 0, shape=(q ** 2, q ** 2), format='csr')

    # Clock field: X + X†
    X = np.zeros((q, q))
    for s in range(q):
        X[(s + 1) % q, s] = 1.0
    XpXd = csr_matrix(X + X.T)

    H = csr_matrix((dim, dim))

    # Build bond list
    bonds = set()
    for y in range(Ly):
        for x in range(Lx):
            site = y * Lx + x
            if x + 1 < Lx:
                nbr = y * Lx + (x + 1)
                bonds.add((min(site, nbr), max(site, nbr)))
            nbr_v = ((y + 1) % Ly) * Lx + x
            if nbr_v != site:
                bonds.add((min(site, nbr_v), max(site, nbr_v)))

    for (i, j) in bonds:
        if j == i + 1:
            # Adjacent sites — use Kronecker product
            left = q ** i
            right = q ** (n - j - 1)
            op = potts_op
            if left > 1:
                op = sp_kron(sp_eye(left), op, format='csr')
            if right > 1:
                op = sp_kron(op, sp_eye(right), format='csr')
        else:
            # Non-adjacent: vectorized diagonal construction
            # For each basis state idx, extract s_i and s_j
            indices = np.arange(dim)
            si = (indices // (q ** i)) % q
            sj = (indices // (q ** j)) % q
            diag_vals = (si == sj).astype(float)
            op = diags(diag_vals, 0, shape=(dim, dim), format='csr')
        H = H - op

    for i in range(n):
        left = q ** i
        right = q ** (n - i - 1)
        op = XpXd.copy()
        if left > 1:
            op = sp_kron(sp_eye(left), op, format='csr')
        if right > 1:
            op = sp_kron(op, sp_eye(right), format='csr')
        H = H - g * op

    return H


def get_gap_and_E0(Lx, Ly, q, g):
    """Get energy gap and ground state energy."""
    H = build_H_cylinder(Lx, Ly, q, g)
    dim = H.shape[0]
    if dim < 500:
        evals = np.linalg.eigvalsh(H.toarray())
        return evals[1] - evals[0], evals[0]
    k = min(6, dim - 1)
    evals, _ = eigsh(H, k=k, which='SA')
    evals = np.sort(evals)
    return evals[1] - evals[0], evals[0]


def compute_order_param(Lx, Ly, q, g):
    """Compute <delta(s_i, s_j)> for middle rung bond."""
    H = build_H_cylinder(Lx, Ly, q, g)
    dim = H.shape[0]
    n = Lx * Ly
    if dim < 500:
        evals, evecs = np.linalg.eigh(H.toarray())
        psi = evecs[:, 0]
    else:
        _, evecs = eigsh(H, k=1, which='SA')
        psi = evecs[:, 0]

    mid_x = Lx // 2
    i = 0 * Lx + mid_x
    j = 1 * Lx + mid_x

    indices = np.arange(dim)
    si = (indices // (q ** i)) % q
    sj = (indices // (q ** j)) % q
    return float(np.sum(psi ** 2 * (si == sj)))


# ---- MAIN ----
print("=" * 60, flush=True)
print("Sprint 071c: Exact diag energy gap on Ly=2 cylinder", flush=True)
print("=" * 60, flush=True)

results = {
    'experiment': '071c_cylinder_gap',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'Ly': 2,
}

# --- Part 1: q=2 gap*Lx crossing ---
print("\n--- q=2: gap*Lx crossing ---", flush=True)
g_range_q2 = np.linspace(0.3, 0.7, 41)
gap_data_q2 = {}

for Lx in [4, 5, 6, 7]:
    n = Lx * 2
    dim = 2 ** n
    t0 = time.time()
    data = []
    for g in g_range_q2:
        gap, E0 = get_gap_and_E0(Lx, 2, 2, g)
        data.append({'g': round(float(g), 4), 'gap': float(gap),
                     'gap_Lx': float(gap * Lx), 'E0': float(E0)})
    dt = time.time() - t0
    gapLx = [d['gap_Lx'] for d in data]
    print(f"  Lx={Lx} (dim={dim}): done {dt:.1f}s, gap*Lx=[{min(gapLx):.3f}, {max(gapLx):.3f}]", flush=True)
    gap_data_q2[Lx] = data

results['q2_gap'] = {str(k): v for k, v in gap_data_q2.items()}

# Find crossings
print("\n  Crossings:", flush=True)
q2_crossings = []
Lx_list = sorted(gap_data_q2.keys())
for i in range(len(Lx_list) - 1):
    Lx1, Lx2 = Lx_list[i], Lx_list[i + 1]
    g1 = np.array([d['gap_Lx'] for d in gap_data_q2[Lx1]])
    g2 = np.array([d['gap_Lx'] for d in gap_data_q2[Lx2]])
    diff = g1 - g2
    for j in range(len(diff) - 1):
        if diff[j] * diff[j + 1] < 0:
            g_cross = g_range_q2[j] - diff[j] * (g_range_q2[j + 1] - g_range_q2[j]) / (diff[j + 1] - diff[j])
            print(f"  Lx={Lx1},{Lx2}: g_c = {g_cross:.4f}", flush=True)
            q2_crossings.append({'Lx_pair': [Lx1, Lx2], 'g_c': float(g_cross)})
            break

results['q2_crossings'] = q2_crossings
if q2_crossings:
    gc_q2 = np.mean([c['g_c'] for c in q2_crossings])
    results['q2_gc'] = gc_q2
    print(f"  Average g_c(q=2, Ly=2 cyl) = {gc_q2:.4f}", flush=True)

# --- Part 2: q=5 gap*Lx crossing ---
print("\n--- q=5: gap*Lx crossing ---", flush=True)
g_range_q5 = np.linspace(0.5, 1.2, 36)

gap_data_q5 = {}
for Lx in [3, 4]:
    n = Lx * 2
    dim = 5 ** n
    t0 = time.time()
    data = []
    for gi, g in enumerate(g_range_q5):
        gap, E0 = get_gap_and_E0(Lx, 2, 5, g)
        data.append({'g': round(float(g), 4), 'gap': float(gap),
                     'gap_Lx': float(gap * Lx), 'E0': float(E0)})
        if (gi + 1) % 12 == 0:
            print(f"    Lx={Lx}: {gi+1}/{len(g_range_q5)} done", flush=True)
    dt = time.time() - t0
    gapLx = [d['gap_Lx'] for d in data]
    print(f"  Lx={Lx} (dim={dim}): done {dt:.1f}s, gap*Lx=[{min(gapLx):.3f}, {max(gapLx):.3f}]", flush=True)
    gap_data_q5[Lx] = data

# Try Lx=5 (10M dim) with timing test first
print("\n  Testing Lx=5 (dim=9.8M)...", flush=True)
t0 = time.time()
try:
    gap_test, E0_test = get_gap_and_E0(5, 2, 5, 0.8)
    dt_test = time.time() - t0
    print(f"    Single point: {dt_test:.1f}s", flush=True)
    if dt_test < 120:  # If under 2 min/point, do full scan
        npts = max(10, int(120 * 36 / dt_test))  # scale down points
        npts = min(npts, 20)
        g_range_q5_big = np.linspace(0.55, 1.1, npts)
        print(f"    Running {npts} points for Lx=5...", flush=True)
        data = []
        for gi, g in enumerate(g_range_q5_big):
            gap, E0 = get_gap_and_E0(5, 2, 5, g)
            data.append({'g': round(float(g), 4), 'gap': float(gap),
                         'gap_Lx': float(gap * 5), 'E0': float(E0)})
            if (gi + 1) % 5 == 0:
                print(f"      {gi+1}/{npts} done", flush=True)
        gapLx = [d['gap_Lx'] for d in data]
        print(f"  Lx=5 (dim=9.8M): gap*Lx=[{min(gapLx):.3f}, {max(gapLx):.3f}]", flush=True)
        gap_data_q5[5] = data
    else:
        print(f"    Too slow ({dt_test:.0f}s/pt), skipping full scan", flush=True)
except Exception as e:
    dt_test = time.time() - t0
    print(f"    Failed ({e}) after {dt_test:.0f}s", flush=True)

results['q5_gap'] = {str(k): v for k, v in gap_data_q5.items()}

# Find crossings
print("\n  Crossings:", flush=True)
q5_crossings = []
Lx_list_5 = sorted(gap_data_q5.keys())
for i in range(len(Lx_list_5) - 1):
    Lx1, Lx2 = Lx_list_5[i], Lx_list_5[i + 1]
    # Interpolate to common g grid
    g1_data = gap_data_q5[Lx1]
    g2_data = gap_data_q5[Lx2]
    g_common = np.array([d['g'] for d in g1_data])
    gLx1 = np.array([d['gap_Lx'] for d in g1_data])
    gLx2_interp = np.interp(g_common, [d['g'] for d in g2_data], [d['gap_Lx'] for d in g2_data])
    diff = gLx1 - gLx2_interp
    found = False
    for j in range(len(diff) - 1):
        if diff[j] * diff[j + 1] < 0:
            g_cross = g_common[j] - diff[j] * (g_common[j + 1] - g_common[j]) / (diff[j + 1] - diff[j])
            print(f"  Lx={Lx1},{Lx2}: g_c = {g_cross:.4f}", flush=True)
            q5_crossings.append({'Lx_pair': [Lx1, Lx2], 'g_c': float(g_cross)})
            found = True
            break
    if not found:
        print(f"  Lx={Lx1},{Lx2}: NO crossing found", flush=True)

results['q5_crossings'] = q5_crossings
if q5_crossings:
    gc_q5 = np.mean([c['g_c'] for c in q5_crossings])
    results['q5_gc'] = gc_q5
    print(f"  Average g_c(q=5, Ly=2 cyl) = {gc_q5:.4f}", flush=True)

# --- Part 3: Order parameter ---
print("\n--- Order parameter <delta(s_i,s_j)> ---", flush=True)

# q=2
gc_q2 = results.get('q2_gc', 0.45)
g_op_q2 = np.linspace(0.2, 0.8, 13)
print(f"\n  q=2, Lx=6:", flush=True)
op_q2 = []
for g in g_op_q2:
    op = compute_order_param(6, 2, 2, g)
    op_q2.append({'g': round(float(g), 4), 'delta_ij': float(op)})
print(f"  Near g_c={gc_q2:.3f}:", flush=True)
for d in op_q2:
    marker = " <-- g_c" if abs(d['g'] - gc_q2) < 0.05 else ""
    print(f"    g={d['g']:.3f}: <delta>={d['delta_ij']:.5f}{marker}", flush=True)
results['q2_order_param'] = {'Lx': 6, 'data': op_q2}

# q=5
gc_q5_val = results.get('q5_gc', 0.8)
g_op_q5 = np.linspace(0.3, 1.3, 21)
print(f"\n  q=5, Lx=4:", flush=True)
op_q5 = []
for g in g_op_q5:
    t0 = time.time()
    op = compute_order_param(4, 2, 5, g)
    dt = time.time() - t0
    op_q5.append({'g': round(float(g), 4), 'delta_ij': float(op)})
print(f"  Near g_c={gc_q5_val:.3f}:", flush=True)
for d in op_q5:
    marker = " <-- g_c" if abs(d['g'] - gc_q5_val) < 0.08 else ""
    print(f"    g={d['g']:.3f}: <delta>={d['delta_ij']:.5f}{marker}", flush=True)
results['q5_order_param'] = {'Lx': 4, 'data': op_q5}

# Check for jump (first-order signal)
print(f"\n  Order parameter analysis:", flush=True)
for label, op_data, qv, gc_val in [('q2', op_q2, 2, gc_q2), ('q5', op_q5, 5, gc_q5_val)]:
    gs = np.array([d['g'] for d in op_data])
    ops = np.array([d['delta_ij'] for d in op_data])
    dop_dg = np.gradient(ops, gs)
    max_slope = np.max(np.abs(dop_dg))
    max_slope_g = gs[np.argmax(np.abs(dop_dg))]
    # Value at g_c
    op_at_gc = np.interp(gc_val, gs, ops)
    print(f"  q={qv}: <delta> at g_c = {op_at_gc:.4f} (random={1/qv:.3f}, ordered=1)", flush=True)
    print(f"       max |d<delta>/dg| = {max_slope:.3f} at g={max_slope_g:.3f}", flush=True)
    results[f'{label}_op_analysis'] = {
        'op_at_gc': float(op_at_gc), 'max_slope': float(max_slope),
        'max_slope_g': float(max_slope_g)
    }

# Save
with open("results/sprint_071c_cylinder_gap.json", "w") as f:
    json.dump(results, f, indent=2, default=str)
print("\nSaved to results/sprint_071c_cylinder_gap.json", flush=True)

from db_utils import record
if 'q2_gc' in results:
    record(sprint=71, model='hybrid_cyl', q=2, n=2,
           quantity='gc_gap', value=results['q2_gc'],
           method='gapLx_crossing_Ly2', notes='Ly=2 cylinder exact diag')
if 'q5_gc' in results:
    record(sprint=71, model='hybrid_cyl', q=5, n=2,
           quantity='gc_gap', value=results['q5_gc'],
           method='gapLx_crossing_Ly2', notes='Ly=2 cylinder exact diag')

print(f"\n{'='*60}", flush=True)
print(f"SUMMARY", flush=True)
print(f"{'='*60}", flush=True)
if 'q2_gc' in results:
    print(f"g_c(q=2, Ly=2 cyl) = {results['q2_gc']:.4f}  [1D: 0.250, 2D: 0.771]", flush=True)
    print(f"  Ratio to 1D: {results['q2_gc']/0.250:.3f}", flush=True)
if 'q5_gc' in results:
    print(f"g_c(q=5, Ly=2 cyl) = {results['q5_gc']:.4f}  [1D: 0.441, 2D: 1.588]", flush=True)
    print(f"  Ratio to 1D: {results['q5_gc']/0.441:.3f}", flush=True)
if 'q2_gc' in results and 'q5_gc' in results:
    r2 = results['q2_gc'] / 0.250
    r5 = results['q5_gc'] / 0.441
    print(f"  Cyl/1D ratios: q=2 {r2:.3f}, q=5 {r5:.3f}", flush=True)
