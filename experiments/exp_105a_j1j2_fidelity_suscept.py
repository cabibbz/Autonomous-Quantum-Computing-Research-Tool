"""Sprint 105a: Fidelity susceptibility at J1-J2 BKT transition.

chi_F(J2) = (2/N) * (1 - |<psi(J2)|psi(J2+dJ2)>|) / dJ2^2

Split: H(J2) = H_NN + J2 * H_NNN  (both in Sz=0 sector, periodic BC, Delta=1)
Scan J2 near J2_c=0.2412 for N=8,10,12,14,16,18,20.

Also measure chi_F at the Ising-like XX transition (Delta scan near Delta_c=1 for J2=0)
as a second-order control with known nu=1/2 (XXZ anisotropy transition).

Actually, for cleaner comparison: measure chi_F for the XX model varying Delta near 0
where the transition is known. But the J1-J2 BKT is the main target.
"""
import numpy as np
from scipy.sparse import lil_matrix, csr_matrix
from gpu_utils import eigsh
import json, time, os
from itertools import combinations
from db_utils import record

results = {
    'experiment': '105a_j1j2_fidelity_suscept',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'data': {},
}

def save():
    outpath = os.path.join(os.path.dirname(__file__), 'results', 'sprint_105a_j1j2_chi_F.json')
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    with open(outpath, 'w') as f:
        json.dump(results, f, indent=2, default=str)

# --- Sz=0 sector basis ---
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

# --- Build NN and NNN parts separately (Heisenberg, periodic BC) ---
def build_j1j2_parts(N):
    """H(J2) = H_NN + J2 * H_NNN where both are Heisenberg (Delta=1)."""
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
                # Sz_i Sz_j
                diag += (si - 0.5) * (sj - 0.5)
                # S+_i S-_j
                if si == 0 and sj == 1:
                    new_state = (state | (1 << i)) & ~(1 << j)
                    ni = idx.get(int(new_state))
                    if ni is not None:
                        rows.append(k); cols.append(ni); vals.append(0.5)
                # S-_i S+_j
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

# --- Main scan ---
J2_c = 0.2412
dJ2 = 1e-4  # finite difference step
J2_range = 0.08
n_g = 21  # number of scan points

sizes = [8, 10, 12, 14, 16, 18, 20]
J2_vals = np.linspace(J2_c - J2_range, J2_c + J2_range, n_g)

print("=" * 65)
print("Sprint 105a: chi_F at J1-J2 BKT transition (J2_c=0.2412)")
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

    results['data'][str(N)] = {
        'dim': dim,
        'J2_vals': J2_vals.tolist(),
        'chi_F': chi_vals,
        'J2_peak': J2_peak,
        'chi_peak': chi_peak,
        'width_FWHM': width,
        'time_s': dt,
    }

    record(sprint=105, model='J1J2', q=0, n=N,
           quantity='chi_F_max', value=chi_peak, method='BKT_J2')
    record(sprint=105, model='J1J2', q=0, n=N,
           quantity='J2_peak', value=J2_peak, method='BKT_J2')

    print(f"  J2_peak={J2_peak:.4f} (J2_c=0.2412, shift={J2_peak-J2_c:+.4f})")
    print(f"  chi_peak={chi_peak:.2f}, FWHM={width:.4f}, time={dt:.1f}s")

    save()

    # Time guard
    if dt > 55:
        print(f"  WARNING: {dt:.0f}s — skipping larger sizes")
        break

# Summary
print("\n" + "=" * 65)
print("SUMMARY: chi_F_max vs N at J1-J2 BKT")
print(f"{'N':>4} {'dim':>8} {'J2_peak':>8} {'chi_max':>10} {'FWHM':>8} {'time':>6}")
print("-" * 55)
for N_str in sorted(results['data'].keys(), key=int):
    d = results['data'][N_str]
    print(f"{N_str:>4} {d['dim']:>8} {d['J2_peak']:>8.4f} {d['chi_peak']:>10.2f} "
          f"{d['width_FWHM']:>8.4f} {d['time_s']:>5.1f}s")

# Quick alpha estimate from log-log
Ns_done = sorted([int(k) for k in results['data'].keys()])
if len(Ns_done) >= 2:
    log_N = np.log([float(n) for n in Ns_done])
    log_chi = np.log([results['data'][str(n)]['chi_peak'] for n in Ns_done])
    if len(Ns_done) >= 3:
        # Fit log(chi) = alpha*log(N) + const
        coeffs = np.polyfit(log_N, log_chi, 1)
        alpha = coeffs[0]
        print(f"\nGlobal alpha (log-log fit, {len(Ns_done)} pts): {alpha:.3f}")
    # Pairwise
    print("\nPairwise alpha:")
    for i in range(len(Ns_done)-1):
        N1, N2 = Ns_done[i], Ns_done[i+1]
        a = (log_chi[i+1] - log_chi[i]) / (log_N[i+1] - log_N[i])
        print(f"  ({N1},{N2}): alpha = {a:.3f}")

save()
print(f"\nResults saved.")
