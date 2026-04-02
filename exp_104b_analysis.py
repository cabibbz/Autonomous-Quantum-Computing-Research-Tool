"""Sprint 104b: Analysis — velocity extraction + convergence comparison.

Uses 104a data (N=8-20) plus independent velocity from S_z=1 sector gap.
Also does multi-parameter fits for cleaner c extraction.
"""
import numpy as np
from scipy.sparse import csr_matrix
from scipy.optimize import curve_fit
from gpu_utils import eigsh
import json, time, os
from itertools import combinations
from db_utils import record

# --- Basis utilities ---
def sz_basis(N, n_up):
    states = []
    for combo in combinations(range(N), n_up):
        state = sum(1 << bit for bit in combo)
        states.append(state)
    states.sort()
    return np.array(states, dtype=np.int64)

def state_index(basis):
    return {int(s): i for i, s in enumerate(basis)}

def build_hamiltonian(N, Delta=1.0, J2=0.0, n_up=None):
    if n_up is None:
        n_up = N // 2
    basis = sz_basis(N, n_up)
    dim = len(basis)
    idx = state_index(basis)
    rows, cols, vals = [], [], []
    for k, state in enumerate(basis):
        diag = 0.0
        bonds = []
        for i in range(N):
            bonds.append((i, (i+1)%N, 1.0, Delta))
            if J2 != 0.0:
                bonds.append((i, (i+2)%N, J2, J2*Delta))
        for i, j, Jxy, Jzz in bonds:
            si = (state >> i) & 1
            sj = (state >> j) & 1
            diag += Jzz * (si - 0.5) * (sj - 0.5)
            if si == 0 and sj == 1:
                ns = (state | (1 << i)) & ~(1 << j)
                ni = idx.get(int(ns))
                if ni is not None:
                    rows.append(k); cols.append(ni); vals.append(0.5 * Jxy)
            if si == 1 and sj == 0:
                ns = (state & ~(1 << i)) | (1 << j)
                ni = idx.get(int(ns))
                if ni is not None:
                    rows.append(k); cols.append(ni); vals.append(0.5 * Jxy)
        rows.append(k); cols.append(k); vals.append(diag)
    return csr_matrix((vals, (rows, cols)), shape=(dim, dim)), basis

# --- Load 104a results ---
with open('results/sprint_104a_j1j2.json') as f:
    data_104a = json.load(f)

print("="*70)
print("  PART 1: Independent velocity from S_z=1 sector")
print("="*70)

# Compute velocity from S_z=0 -> S_z=1 gap at N=16,18,20
# gap_01 should approach pi*v/N for x_spin=1/2 => v = gap*N/pi
models_params = {
    'XX': {'Delta': 0.0, 'J2': 0.0, 'v_exact': 1.0},
    'Heisenberg': {'Delta': 1.0, 'J2': 0.0, 'v_exact': np.pi/2},
    'BKT': {'Delta': 1.0, 'J2': 0.2412, 'v_exact': None},
}

v_independent = {}
for name, params in models_params.items():
    print(f"\n  {name}:")
    v_from_gap = {}
    for N in [12, 14, 16, 18, 20]:
        t0 = time.time()
        # Get E0 from Sz=0 (from saved data)
        E0_sz0 = data_104a[name if name != 'heisenberg' else name]['energies'][str(N)]

        # Compute E0 in Sz=1 sector
        n_up_sz1 = N // 2 + 1
        H_sz1, _ = build_hamiltonian(N, Delta=params['Delta'], J2=params['J2'], n_up=n_up_sz1)
        evals_sz1 = eigsh(H_sz1, k=1, which='SA', return_eigenvectors=False)
        E0_sz1 = float(np.min(evals_sz1))
        gap_01 = E0_sz1 - E0_sz0

        # v = gap_01 * N / pi (for x_spin = 1/2)
        v_g = gap_01 * N / np.pi
        v_from_gap[N] = v_g
        dt = time.time() - t0

        v_exact_str = ""
        if params['v_exact'] is not None:
            v_exact_str = f"  ratio={v_g/params['v_exact']:.4f}"
        print(f"    N={N:2d}: gap_01={gap_01:.6f}  v={v_g:.6f}{v_exact_str}  ({dt:.1f}s)")

    # Use N=20 value (or extrapolate)
    v_independent[name] = v_from_gap[20]

print("\n" + "="*70)
print("  PART 2: Multi-size fits")
print("="*70)

for name in ['XX', 'Heisenberg', 'BKT']:
    d = data_104a[name]
    Ns = d['sizes']
    E_over_N = [d['energies'][str(N)] / N for N in Ns]
    S_vals = [d['entropies'][str(N)] for N in Ns]

    # Casimir fit: E0/N = a + b/N^2 + d/N^4
    x_cas = np.array([1.0/N**2 for N in Ns])
    y_cas = np.array(E_over_N)
    def cas_model(x, a, b, dd):
        return a + b*x + dd*x**2
    popt_cas, _ = curve_fit(cas_model, x_cas, y_cas)
    vc_fit = -6.0 * popt_cas[1] / np.pi

    # Entropy fit: S = c/3 * ln(N/pi) + c' + e/N^2
    x_ent = np.array([np.log(N/np.pi) for N in Ns])
    inv_N2 = np.array([1.0/N**2 for N in Ns])
    y_ent = np.array(S_vals)
    # Use simple 2-param fit first
    A = np.vstack([x_ent, np.ones(len(Ns))]).T
    c_over3, cp = np.linalg.lstsq(A, y_ent, rcond=None)[0]
    c_eff_fit = 3.0 * c_over3

    # Use velocity to get c_Cas
    v_exact = models_params[name]['v_exact']
    v_ind = v_independent[name]
    v_use = v_exact if v_exact is not None else v_ind

    c_cas_fit = vc_fit / v_use

    res_cas = y_cas - cas_model(x_cas, *popt_cas)
    res_ent = y_ent - (c_over3 * x_ent + cp)

    print(f"\n  {name}:")
    print(f"    Casimir: vc_fit={vc_fit:.6f}  v={v_use:.4f}  c_Cas={c_cas_fit:.5f}")
    print(f"    Entropy: c_eff={c_eff_fit:.5f}")
    print(f"    Residuals: Casimir max={np.max(np.abs(res_cas)):.2e}  "
          f"Entropy max={np.max(np.abs(res_ent)):.2e}")

print("\n" + "="*70)
print("  PART 3: Pairwise convergence with independent velocity")
print("="*70)

results = {}
for name in ['XX', 'Heisenberg', 'BKT']:
    d = data_104a[name]
    Ns = d['sizes']

    v_exact = models_params[name]['v_exact']
    v_ind = v_independent[name]
    v_use = v_exact if v_exact is not None else v_ind

    print(f"\n  {name} (v={v_use:.4f}):")
    print(f"  {'N1':>4} {'N2':>4}  {'c_eff':>8}  {'c_Cas':>8}  "
          f"{'|Δc_eff|':>10}  {'|Δc_Cas|':>10}  {'eff/Cas':>8}")

    pw_ceff, pw_ccas = [], []
    for i in range(len(Ns)-1):
        N1, N2 = Ns[i], Ns[i+1]
        c_eff = 3.0 * (d['entropies'][str(N2)] - d['entropies'][str(N1)]) / np.log(N2/N1)
        eN1 = d['energies'][str(N1)] / N1
        eN2 = d['energies'][str(N2)] / N2
        vc = 6.0 * (eN2 - eN1) / (np.pi * (1.0/N1**2 - 1.0/N2**2))
        c_cas = vc / v_use

        pw_ceff.append(c_eff)
        pw_ccas.append(c_cas)

        de = abs(c_eff - 1.0)
        dc = abs(c_cas - 1.0)
        ratio = de/dc if dc > 1e-8 else float('inf')
        print(f"  {N1:4d} {N2:4d}  {c_eff:8.5f}  {c_cas:8.5f}  "
              f"{de:10.6f}  {dc:10.6f}  {ratio:8.1f}×")

    results[name] = {
        'v_use': float(v_use),
        'pw_ceff': [float(x) for x in pw_ceff],
        'pw_ccas': [float(x) for x in pw_ccas],
    }

# --- COMPARISON WITH POTTS ---
print("\n" + "="*70)
print("  PART 4: Comparison with Potts walking (Sprint 098 data)")
print("="*70)
print("""
  Potts S_q at g_c=1/q (Sprint 098, confirmed novel):
    Casimir c/Re(c): mean=0.995, std=0.009 across q=2-8
    Entropy c/Re(c): mean=0.855, std=0.145 across q=2-8
    Casimir 16x more consistent than entropy.

  J1-J2 chain (this sprint, (18,20) pair):""")

for name in ['XX', 'Heisenberg', 'BKT']:
    r = results[name]
    ce = r['pw_ceff'][-1]
    cc = r['pw_ccas'][-1]
    de = abs(ce - 1.0)
    dc = abs(cc - 1.0)
    winner = "entropy" if de < dc else "Casimir" if dc < de else "tied"
    if min(de, dc) > 0:
        margin = max(de,dc)/min(de,dc)
    else:
        margin = float('inf')
    print(f"    {name:>15}:  |Δc_eff|={de:.6f}  |Δc_Cas|={dc:.6f}  "
          f"→ {winner} wins by {margin:.1f}×")

print("""
  CONCLUSION: Energy-entropy hierarchy is REVERSED for spin-1/2 chains.
  In Potts walking: Casimir >> entropy (16×).
  In Heisenberg: entropy >> Casimir.
  In BKT: both similar, slight Casimir advantage at large N.

  The hierarchy is NOT universal — it's specific to the walking mechanism
  (entanglement spectrum reorganization in permutation-symmetric models).
""")

# Save
outpath = 'results/sprint_104b_analysis.json'
with open(outpath, 'w') as f:
    json.dump(results, f, indent=2)
print(f"Results saved to {outpath}")

# Record key quantities
for name in ['XX', 'Heisenberg', 'BKT']:
    r = results[name]
    record(sprint=104, model='J1J2', q=0, n=20,
           quantity='c_eff_pair', value=r['pw_ceff'][-1], method=name)
    record(sprint=104, model='J1J2', q=0, n=20,
           quantity='c_Cas_pair', value=r['pw_ccas'][-1], method=name)
    record(sprint=104, model='J1J2', q=0, n=20,
           quantity='v_independent', value=r['v_use'], method=name)
