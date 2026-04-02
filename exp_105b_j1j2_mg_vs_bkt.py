"""Sprint 105b: Resolve J1-J2 chi_F — MG peak vs BKT invisibility.

Two focused scans:
1. Fine scan near J2_c=0.2412 (BKT) — is chi_F featureless?
2. Fine scan near J2=0.5 (Majumdar-Ghosh) — resolve the true peak and FSS

Also: handle level crossings by computing gap to detect them.
"""
import numpy as np
from scipy.sparse import csr_matrix
from gpu_utils import eigsh
import json, time, os
from itertools import combinations
from db_utils import record

results = {
    'experiment': '105b_mg_vs_bkt',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'bkt_region': {},
    'mg_region': {},
}

def save():
    outpath = os.path.join(os.path.dirname(__file__), 'results', 'sprint_105b_mg_vs_bkt.json')
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

def gs_and_gap(H_NN, H_NNN, J2):
    H = H_NN + J2 * H_NNN
    evals, evecs = eigsh(H, k=2, which='SA')
    idx = np.argsort(evals)
    return evecs[:, idx[0]], float(evals[idx[1]] - evals[idx[0]])

def chi_F_point(H_NN, H_NNN, J2, dJ2, N):
    H1 = H_NN + J2 * H_NNN
    H2 = H_NN + (J2 + dJ2) * H_NNN
    _, v1 = eigsh(H1, k=1, which='SA')
    _, v2 = eigsh(H2, k=1, which='SA')
    overlap = abs(np.vdot(v1[:, 0], v2[:, 0]))
    overlap = min(overlap, 1.0 - 1e-15)
    return 2.0 * (1.0 - overlap) / (dJ2**2 * N)

dJ2 = 1e-4
sizes = [8, 10, 12, 14, 16, 18, 20]

# ============================================
# Part 1: BKT region scan (J2 = 0.15 to 0.35)
# ============================================
print("=" * 65)
print("Part 1: BKT region (J2 = 0.15 to 0.35)")
print("=" * 65)

J2_bkt = np.linspace(0.15, 0.35, 25)

for N in sizes:
    t0 = time.time()
    H_NN, H_NNN, basis = build_j1j2_parts(N)
    dim = H_NN.shape[0]

    chi_vals, gap_vals = [], []
    for J2 in J2_bkt:
        chi = chi_F_point(H_NN, H_NNN, J2, dJ2, N)
        _, gap = gs_and_gap(H_NN, H_NNN, J2)
        chi_vals.append(float(chi))
        gap_vals.append(float(gap))

    dt = time.time() - t0

    i_peak = np.argmax(chi_vals)
    i_min_gap = np.argmin(gap_vals)

    results['bkt_region'][str(N)] = {
        'dim': dim,
        'J2_vals': J2_bkt.tolist(),
        'chi_F': chi_vals,
        'gaps': gap_vals,
        'chi_peak': float(chi_vals[i_peak]),
        'J2_chi_peak': float(J2_bkt[i_peak]),
        'min_gap': float(gap_vals[i_min_gap]),
        'J2_min_gap': float(J2_bkt[i_min_gap]),
        'time_s': dt,
    }

    # Is chi_F monotonic in this region?
    monotonic = all(chi_vals[i+1] >= chi_vals[i] - 1e-6 for i in range(len(chi_vals)-1))

    print(f"N={N:2d} dim={dim:>6d}  chi_peak={chi_vals[i_peak]:.4f} at J2={J2_bkt[i_peak]:.3f}  "
          f"min_gap={gap_vals[i_min_gap]:.4f} at J2={J2_bkt[i_min_gap]:.3f}  "
          f"{'MONOTONIC' if monotonic else 'has peak'}  ({dt:.1f}s)")

    save()
    if dt > 40:
        print("  Skipping larger sizes")
        break

# ============================================
# Part 2: MG region scan (J2 = 0.40 to 0.55)
# ============================================
print(f"\n{'=' * 65}")
print("Part 2: MG region (J2 = 0.40 to 0.55)")
print("=" * 65)

J2_mg = np.linspace(0.40, 0.55, 31)

mg_peaks = {}

for N in sizes:
    t0 = time.time()
    H_NN, H_NNN, basis = build_j1j2_parts(N)
    dim = H_NN.shape[0]

    chi_vals, gap_vals = [], []
    for J2 in J2_mg:
        chi = chi_F_point(H_NN, H_NNN, J2, dJ2, N)
        _, gap = gs_and_gap(H_NN, H_NNN, J2)
        chi_vals.append(float(chi))
        gap_vals.append(float(gap))

    dt = time.time() - t0

    # Filter level crossings: flag if gap < 1e-6
    clean_chi = [(J2_mg[i], chi_vals[i]) for i in range(len(chi_vals))
                 if gap_vals[i] > 1e-4]
    if clean_chi:
        i_clean_peak = max(range(len(clean_chi)), key=lambda i: clean_chi[i][1])
        J2_peak_clean = clean_chi[i_clean_peak][0]
        chi_peak_clean = clean_chi[i_clean_peak][1]
    else:
        J2_peak_clean = float('nan')
        chi_peak_clean = float('nan')

    results['mg_region'][str(N)] = {
        'dim': dim,
        'J2_vals': J2_mg.tolist(),
        'chi_F': chi_vals,
        'gaps': gap_vals,
        'chi_peak_clean': float(chi_peak_clean),
        'J2_peak_clean': float(J2_peak_clean),
        'time_s': dt,
    }

    mg_peaks[N] = chi_peak_clean

    # Check for level crossings
    n_crossings = sum(1 for g in gap_vals if g < 1e-4)
    cross_flag = f"  [{n_crossings} level crossings]" if n_crossings > 0 else ""

    record(sprint=105, model='J1J2', q=0, n=N,
           quantity='chi_F_MG_clean', value=chi_peak_clean, method='MG_J2')

    print(f"N={N:2d} dim={dim:>6d}  chi_peak(clean)={chi_peak_clean:.4f} "
          f"at J2={J2_peak_clean:.4f}{cross_flag}  ({dt:.1f}s)")

    save()
    if dt > 40:
        print("  Skipping larger sizes")
        break

# FSS of MG peak
print(f"\n{'=' * 65}")
print("FSS of MG chi_F peak (level-crossing filtered)")
print("=" * 65)

Ns_mg = sorted(mg_peaks.keys())
if len(Ns_mg) >= 3:
    log_N = np.log([float(n) for n in Ns_mg])
    log_chi = np.log([mg_peaks[n] for n in Ns_mg])
    coeffs = np.polyfit(log_N, log_chi, 1)
    alpha_mg = coeffs[0]
    print(f"Global alpha (MG peak): {alpha_mg:.3f}")
    print("\nPairwise alpha:")
    for i in range(len(Ns_mg)-1):
        N1, N2 = Ns_mg[i], Ns_mg[i+1]
        a = (log_chi[i+1] - log_chi[i]) / (log_N[i+1] - log_N[i])
        print(f"  ({N1},{N2}): alpha = {a:.3f}")

    results['mg_alpha_global'] = float(alpha_mg)

# Compare with Potts
print(f"\n{'=' * 65}")
print("COMPARISON: chi_F alpha exponents")
print("=" * 65)
print(f"  J1-J2 BKT (J2_c=0.2412):  NO PEAK detectable at N<=20")
print(f"  J1-J2 MG  (J2~0.49):      alpha = {alpha_mg:.3f}" if len(Ns_mg) >= 3 else "")
print(f"  S_q Potts q=2 (2nd order): alpha = 0.98  (nu=1.009)")
print(f"  S_q Potts q=3 (2nd order): alpha = 1.38  (nu=0.841)")
print(f"  S_q Potts q=5 (walking):   alpha = 2.09  (super-1st-order)")
print(f"  S_q Potts q=7 (walking):   alpha = 2.65  (super-1st-order)")

save()
print(f"\nResults saved.")
