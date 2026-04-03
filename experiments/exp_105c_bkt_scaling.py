"""Sprint 105c: chi_F at exactly J2_c=0.2412 for all N.

Extract the BKT scaling: how does chi_F(J2_c, N) grow with N?
Literature predicts chi_F ~ N^2 / (ln N)^2 for BKT.
Compare with Potts walking alpha=2.09 and MG saturation.

Also: compute chi_F at J2=0.5 (MG exact) to verify saturation.
"""
import numpy as np
from scipy.sparse import csr_matrix
from gpu_utils import eigsh
import json, time, os
from itertools import combinations
from db_utils import record

results = {
    'experiment': '105c_bkt_scaling',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'bkt_exact': {},
    'mg_exact': {},
}

def save():
    outpath = os.path.join(os.path.dirname(__file__), 'results', 'sprint_105c_bkt_scaling.json')
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    with open(outpath, 'w') as f:
        json.dump(results, f, indent=2, default=str)

def sz0_basis(N):
    n_up = N // 2
    states = []
    for combo in combinations(range(N), n_up):
        state = sum(1 << bit for bit in combo)
        states.append(state)
    states.sort()
    return np.array(states, dtype=np.int64)

def state_index(basis):
    return {int(s): i for i, s in enumerate(basis)}

def build_j1j2_parts(N):
    basis = sz0_basis(N)
    dim = len(basis)
    idx = state_index(basis)
    def build_bond_matrix(bond_list):
        rows, cols, vals = [], [], []
        for k, state in enumerate(basis):
            diag = 0.0
            for i, j in bond_list:
                si = (state >> i) & 1
                sj = (state >> j) & 1
                diag += (si - 0.5) * (sj - 0.5)
                if si == 0 and sj == 1:
                    new_state = (state | (1 << i)) & ~(1 << j)
                    ni = idx.get(int(new_state))
                    if ni is not None:
                        rows.append(k); cols.append(ni); vals.append(0.5)
                if si == 1 and sj == 0:
                    new_state = (state & ~(1 << i)) | (1 << j)
                    ni = idx.get(int(new_state))
                    if ni is not None:
                        rows.append(k); cols.append(ni); vals.append(0.5)
            rows.append(k); cols.append(k); vals.append(diag)
        return csr_matrix((vals, (rows, cols)), shape=(dim, dim))

    nn_bonds = [(i, (i + 1) % N) for i in range(N)]
    nnn_bonds = [(i, (i + 2) % N) for i in range(N)]
    return build_bond_matrix(nn_bonds), build_bond_matrix(nnn_bonds), basis

def chi_F_at(H_NN, H_NNN, J2, dJ2, N):
    H1 = H_NN + J2 * H_NNN
    H2 = H_NN + (J2 + dJ2) * H_NNN
    _, v1 = eigsh(H1, k=1, which='SA')
    _, v2 = eigsh(H2, k=1, which='SA')
    overlap = abs(np.vdot(v1[:, 0], v2[:, 0]))
    overlap = min(overlap, 1.0 - 1e-15)
    return 2.0 * (1.0 - overlap) / (dJ2**2 * N)

dJ2 = 1e-4
J2_c = 0.2412
J2_mg = 0.5

sizes = [8, 10, 12, 14, 16, 18, 20]

print("=" * 65)
print("Sprint 105c: chi_F scaling at exact J2 values")
print("=" * 65)

# Measure at both J2_c and J2_mg for each size
chi_bkt = {}
chi_mg = {}

for N in sizes:
    t0 = time.time()
    H_NN, H_NNN, basis = build_j1j2_parts(N)
    dim = H_NN.shape[0]

    # At BKT point
    cf_bkt = chi_F_at(H_NN, H_NNN, J2_c, dJ2, N)
    # At MG point (slightly offset to avoid exact level crossing)
    cf_mg_below = chi_F_at(H_NN, H_NNN, 0.498, dJ2, N)
    cf_mg_above = chi_F_at(H_NN, H_NNN, 0.502, dJ2, N)

    dt = time.time() - t0

    chi_bkt[N] = float(cf_bkt)
    chi_mg[N] = float(cf_mg_below)  # below MG to avoid crossing

    results['bkt_exact'][str(N)] = {'chi_F': float(cf_bkt), 'dim': dim}
    results['mg_exact'][str(N)] = {
        'chi_F_below': float(cf_mg_below),
        'chi_F_above': float(cf_mg_above),
        'dim': dim,
    }

    record(sprint=105, model='J1J2', q=0, n=N,
           quantity='chi_F_at_J2c', value=cf_bkt, method='BKT_exact')

    print(f"N={N:2d}  dim={dim:>6d}  chi_F(J2_c)={cf_bkt:.6f}  "
          f"chi_F(0.498)={cf_mg_below:.4f}  chi_F(0.502)={cf_mg_above:.4f}  ({dt:.1f}s)")

    save()
    if dt > 30:
        print("  Skipping larger sizes")
        break

# Analysis
print(f"\n{'=' * 65}")
print("chi_F(J2_c=0.2412) scaling")
print("=" * 65)

Ns = sorted(chi_bkt.keys())
log_N = np.log([float(n) for n in Ns])
log_chi_bkt = np.log([chi_bkt[n] for n in Ns])
log_chi_mg = np.log([chi_mg[n] for n in Ns])

# BKT fit: chi_F ~ N^alpha
coeffs_bkt = np.polyfit(log_N, log_chi_bkt, 1)
alpha_bkt = coeffs_bkt[0]

# BKT fit with log correction: chi_F ~ N^2 / (ln N)^beta
# log(chi) = 2*log(N) - beta*log(log(N)) + const
# Test: is alpha consistent with 2 if we add log correction?
loglog_N = np.log(np.log([float(n) for n in Ns]))
# Fit: log(chi) = alpha*log(N) + gamma*log(log(N)) + c
A = np.column_stack([log_N, loglog_N, np.ones(len(Ns))])
fit_3param, _, _, _ = np.linalg.lstsq(A, log_chi_bkt, rcond=None)
alpha_corr, gamma, const = fit_3param

print(f"\nSimple power law: chi_F ~ N^{alpha_bkt:.3f}")
print(f"With log correction: chi_F ~ N^{alpha_corr:.3f} / (ln N)^{-gamma:.3f}")
print(f"  (gamma={gamma:.3f}, expected ~-2 for BKT)")

print(f"\nPairwise alpha (BKT):")
for i in range(len(Ns)-1):
    N1, N2 = Ns[i], Ns[i+1]
    a = (log_chi_bkt[i+1] - log_chi_bkt[i]) / (log_N[i+1] - log_N[i])
    print(f"  ({N1},{N2}): {a:.3f}")

# MG scaling
coeffs_mg = np.polyfit(log_N, log_chi_mg, 1)
alpha_mg = coeffs_mg[0]
print(f"\nMG point (J2=0.498): alpha = {alpha_mg:.3f}")
print(f"Pairwise alpha (MG):")
for i in range(len(Ns)-1):
    N1, N2 = Ns[i], Ns[i+1]
    a = (log_chi_mg[i+1] - log_chi_mg[i]) / (log_N[i+1] - log_N[i])
    print(f"  ({N1},{N2}): {a:.3f}")

# Summary comparison table
print(f"\n{'=' * 65}")
print("COMPLETE COMPARISON: chi_F scaling exponents")
print("=" * 65)
print(f"{'Transition':>25}  {'alpha':>8}  {'alpha_eff(last)':>14}  {'Mechanism':>20}")
print(f"{'-'*75}")
print(f"{'J1-J2 BKT (J2_c)':>25}  {alpha_bkt:>8.3f}  {'—':>14}  {'Invisible at N≤20':>20}")
print(f"{'J1-J2 MG (J2=0.5)':>25}  {alpha_mg:>8.3f}  {'—':>14}  {'Saturating':>20}")
print(f"{'Potts q=2 (2nd order)':>25}  {'0.98':>8}  {'0.98':>14}  {'Power law (ν=1)':>20}")
print(f"{'Potts q=3 (2nd order)':>25}  {'1.38':>8}  {'1.38':>14}  {'Power law (ν=0.84)':>20}")
print(f"{'Potts q=5 (walking)':>25}  {'2.09':>8}  {'2.10':>14}  {'Super-1st-order':>20}")
print(f"{'Potts q=7 (walking)':>25}  {'2.65':>8}  {'2.67':>14}  {'Super-1st-order':>20}")

results['alpha_bkt_simple'] = float(alpha_bkt)
results['alpha_bkt_corrected'] = float(alpha_corr)
results['gamma_bkt'] = float(gamma)
results['alpha_mg'] = float(alpha_mg)

save()
print(f"\nResults saved.")
