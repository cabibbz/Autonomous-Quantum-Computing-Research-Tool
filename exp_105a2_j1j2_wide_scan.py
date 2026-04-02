"""Sprint 105a (revised): Wide J2 scan for chi_F at J1-J2 BKT.

Previous scan had peak at edge (J2=0.32). BKT transitions have large FSS shifts.
Widen to J2 = 0.05 to 0.60 to find the true peak location.
"""
import numpy as np
from scipy.sparse import csr_matrix
from gpu_utils import eigsh
import json, time, os
from itertools import combinations
from db_utils import record

results = {
    'experiment': '105a2_j1j2_wide_scan',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'data': {},
}

def save():
    outpath = os.path.join(os.path.dirname(__file__), 'results', 'sprint_105a2_j1j2_wide.json')
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
    H_NN = build_bond_matrix(nn_bonds)
    H_NNN = build_bond_matrix(nnn_bonds)
    return H_NN, H_NNN, basis

def ground_state(H_NN, H_NNN, J2):
    H = H_NN + J2 * H_NNN
    _, vecs = eigsh(H, k=1, which='SA')
    return vecs[:, 0]

# Wide scan
J2_c = 0.2412
dJ2 = 1e-4
J2_vals = np.linspace(0.05, 0.60, 35)  # wide range

sizes = [8, 10, 12, 14, 16, 18, 20]

print("=" * 65)
print("Sprint 105a2: WIDE chi_F scan at J1-J2 BKT (J2=0.05 to 0.60)")
print("=" * 65)

for N in sizes:
    t0 = time.time()
    H_NN, H_NNN, basis = build_j1j2_parts(N)
    dim = H_NN.shape[0]
    t_build = time.time() - t0
    print(f"\nN={N} (dim={dim}, built in {t_build:.1f}s)", flush=True)

    chi_vals = []
    for J2 in J2_vals:
        psi1 = ground_state(H_NN, H_NNN, J2)
        psi2 = ground_state(H_NN, H_NNN, J2 + dJ2)
        overlap = abs(np.vdot(psi1, psi2))
        overlap = min(overlap, 1.0 - 1e-15)
        chi_F = 2.0 * (1.0 - overlap) / (dJ2**2 * N)
        chi_vals.append(float(chi_F))

    dt = time.time() - t0

    # Find peak
    i_peak = np.argmax(chi_vals)
    J2_peak = float(J2_vals[i_peak])
    chi_peak = float(chi_vals[i_peak])

    # FWHM
    half_max = chi_peak / 2.0
    above = np.array(chi_vals) > half_max
    if above.any():
        left = np.where(above)[0][0]
        right = np.where(above)[0][-1]
        width = float(J2_vals[right] - J2_vals[left])
    else:
        width = float('nan')

    # Is peak at edge?
    at_edge = (i_peak == 0) or (i_peak == len(J2_vals) - 1)

    results['data'][str(N)] = {
        'dim': dim,
        'J2_vals': J2_vals.tolist(),
        'chi_F': chi_vals,
        'J2_peak': J2_peak,
        'chi_peak': chi_peak,
        'width_FWHM': width,
        'at_edge': at_edge,
        'time_s': dt,
    }

    record(sprint=105, model='J1J2', q=0, n=N,
           quantity='chi_F_max_wide', value=chi_peak, method='BKT_J2')
    record(sprint=105, model='J1J2', q=0, n=N,
           quantity='J2_peak_wide', value=J2_peak, method='BKT_J2')

    edge_flag = " ** AT EDGE **" if at_edge else ""
    print(f"  J2_peak={J2_peak:.4f} (shift={J2_peak-J2_c:+.4f}){edge_flag}")
    print(f"  chi_peak={chi_peak:.4f}, FWHM={width:.4f}, time={dt:.1f}s")

    # Print profile around peak
    lo = max(0, i_peak - 3)
    hi = min(len(J2_vals), i_peak + 4)
    profile = " | ".join(f"{J2_vals[i]:.3f}:{chi_vals[i]:.4f}" for i in range(lo, hi))
    print(f"  Profile: {profile}")

    save()

    if dt > 55:
        print(f"  WARNING: {dt:.0f}s — skipping larger sizes")
        break

# Summary
print("\n" + "=" * 65)
print("SUMMARY: chi_F_max vs N (wide scan)")
print(f"{'N':>4} {'dim':>8} {'J2_peak':>8} {'shift':>8} {'chi_max':>10} {'FWHM':>8} {'edge':>5}")
print("-" * 60)
for N_str in sorted(results['data'].keys(), key=int):
    d = results['data'][N_str]
    shift = d['J2_peak'] - J2_c
    edge = "YES" if d['at_edge'] else ""
    print(f"{int(N_str):>4} {d['dim']:>8} {d['J2_peak']:>8.4f} {shift:>+8.4f} "
          f"{d['chi_peak']:>10.4f} {d['width_FWHM']:>8.4f} {edge:>5}")

# Alpha from peak values (only non-edge)
Ns_good = sorted([int(k) for k in results['data'].keys()
                   if not results['data'][k]['at_edge']])
if len(Ns_good) >= 2:
    log_N = np.log([float(n) for n in Ns_good])
    log_chi = np.log([results['data'][str(n)]['chi_peak'] for n in Ns_good])
    if len(Ns_good) >= 3:
        coeffs = np.polyfit(log_N, log_chi, 1)
        print(f"\nGlobal alpha (non-edge, {len(Ns_good)} pts): {coeffs[0]:.3f}")
    print("\nPairwise alpha (non-edge):")
    for i in range(len(Ns_good)-1):
        N1, N2 = Ns_good[i], Ns_good[i+1]
        a = (log_chi[i+1] - log_chi[i]) / (log_N[i+1] - log_N[i])
        print(f"  ({N1},{N2}): alpha = {a:.3f}")

save()
print(f"\nResults saved.")
