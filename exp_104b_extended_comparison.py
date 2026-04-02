"""Sprint 104b: Extended sizes (N=22,24) + multi-size fit for velocity extraction.

Key issue from 104a: for BKT, v was extracted circularly from last pair.
Fix: fit E₀/N = a + b/N² + d/N⁴ to get vc=-6b/π from the fit. Then compare
convergence rates between models.

Also: compute v from S_z=1 sector gap (independent of Casimir formula).
"""
import numpy as np
from scipy.sparse import csr_matrix
from scipy.optimize import curve_fit
from gpu_utils import eigsh
import json, time, os
from itertools import combinations
from db_utils import record

# --- S_z sector basis ---
def sz_basis(N, n_up):
    """Return sorted list of basis states with exactly n_up up spins."""
    states = []
    for combo in combinations(range(N), n_up):
        state = sum(1 << bit for bit in combo)
        states.append(state)
    states.sort()
    return np.array(states, dtype=np.int64)

def state_index(basis):
    return {int(s): i for i, s in enumerate(basis)}

def build_hamiltonian(N, Delta=1.0, J2=0.0, n_up=None):
    """Build XXZ + J2 Hamiltonian in specified S_z sector."""
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
            j1 = (i + 1) % N
            bonds.append((i, j1, 1.0, Delta))
            if J2 != 0.0:
                j2 = (i + 2) % N
                bonds.append((i, j2, J2, J2 * Delta))

        for i, j, Jxy, Jzz in bonds:
            si = (state >> i) & 1
            sj = (state >> j) & 1
            diag += Jzz * (si - 0.5) * (sj - 0.5)

            if si == 0 and sj == 1:
                new_state = (state | (1 << i)) & ~(1 << j)
                ni = idx.get(int(new_state))
                if ni is not None:
                    rows.append(k); cols.append(ni); vals.append(0.5 * Jxy)
            if si == 1 and sj == 0:
                new_state = (state & ~(1 << i)) | (1 << j)
                ni = idx.get(int(new_state))
                if ni is not None:
                    rows.append(k); cols.append(ni); vals.append(0.5 * Jxy)

        rows.append(k); cols.append(k); vals.append(diag)

    H = csr_matrix((vals, (rows, cols)), shape=(dim, dim))
    return H, basis

def midchain_entropy(psi, basis, N):
    nA = N // 2
    mask_A = (1 << nA) - 1
    from collections import defaultdict

    a_groups = defaultdict(list)
    for k, state in enumerate(basis):
        a_config = state & mask_A
        b_config = state >> nA
        a_groups[int(a_config)].append((int(b_config), psi[k]))

    a_configs = sorted(a_groups.keys())
    a_idx = {a: i for i, a in enumerate(a_configs)}

    b_set = set()
    for a in a_configs:
        for b, _ in a_groups[a]:
            b_set.add(b)
    b_configs = sorted(b_set)
    b_idx = {b: i for i, b in enumerate(b_configs)}

    M = np.zeros((len(a_configs), len(b_configs)), dtype=complex)
    for a in a_configs:
        ai = a_idx[a]
        for b, amp in a_groups[a]:
            M[ai, b_idx[b]] = amp

    rho_A = M @ M.conj().T
    eigvals = np.linalg.eigvalsh(rho_A)
    eigvals = eigvals[eigvals > 1e-15]
    return float(-np.sum(eigvals * np.log(eigvals)))

# --- Velocity from S_z=1 sector ---
def extract_velocity_sz1(N, Delta=1.0, J2=0.0, E0_sz0=None):
    """Compute velocity from S_z=0 -> S_z=1 gap.
    Gap = 2*pi*v*x_spin / N, where x_spin depends on the CFT.
    """
    n_up = N // 2 - 1  # S_z = +1 sector: one fewer up spin... wait
    # Actually S_z = total_Sz. For up=1, down=0:
    # S_z = n_up - N/2. For n_up = N/2: S_z=0. For n_up = N/2+1: S_z=+1.
    # But we need n_up <= N, and N/2+1 <= N for N>=2.
    n_up_sz1 = N // 2 + 1  # Wait: if up=1, then n_up means number of sites with spin=1 (up)
    # S_z = sum(s_i - 1/2) = n_up - N/2
    # For S_z=1: n_up = N/2 + 1

    H_sz1, _ = build_hamiltonian(N, Delta=Delta, J2=J2, n_up=n_up_sz1)
    evals_sz1 = eigsh(H_sz1, k=1, which='SA', return_eigenvectors=False)
    E0_sz1 = float(np.min(evals_sz1))

    gap_sz = E0_sz1 - E0_sz0
    return gap_sz

# --- Main ---
def run():
    models = {
        'XX': {'Delta': 0.0, 'J2': 0.0, 'v_exact': 1.0},
        'Heisenberg': {'Delta': 1.0, 'J2': 0.0, 'v_exact': np.pi/2},
        'BKT': {'Delta': 1.0, 'J2': 0.2412, 'v_exact': None},
    }

    sizes = [8, 10, 12, 14, 16, 18, 20, 22]
    results = {}

    for name, params in models.items():
        print(f"\n{'='*60}")
        print(f"  {name}  (Delta={params['Delta']}, J2={params['J2']})")
        print(f"{'='*60}")

        energies, entropies, gaps_sz0, gaps_sz01 = {}, {}, {}, {}

        for N in sizes:
            t0 = time.time()

            # S_z=0 ground state
            H, basis = build_hamiltonian(N, Delta=params['Delta'], J2=params['J2'])
            dim = H.shape[0]

            evals, evecs = eigsh(H, k=2, which='SA')
            idx_sort = np.argsort(evals)
            evals = evals[idx_sort]
            evecs = evecs[:, idx_sort]

            E0 = float(evals[0])
            gap = float(evals[1] - evals[0])
            psi = evecs[:, 0]
            S = midchain_entropy(psi, basis, N)

            # S_z=0 -> S_z=1 gap for velocity
            gap_01 = extract_velocity_sz1(N, Delta=params['Delta'],
                                          J2=params['J2'], E0_sz0=E0)

            dt = time.time() - t0

            energies[N] = E0
            entropies[N] = S
            gaps_sz0[N] = gap
            gaps_sz01[N] = gap_01

            print(f"  N={N:2d}  dim(Sz0)={dim:>8d}  E0/N={E0/N:.8f}  "
                  f"gap_sz0={gap:.6f}  gap_01={gap_01:.6f}  S={S:.6f}  ({dt:.1f}s)")

        # --- Velocity extraction from S_z=0→1 gap ---
        # For SU(2)_1 WZW: gap_01 = 2πv·x_spin/N where x_spin=1/2
        # So v = gap_01·N / (2π·(1/2)) = gap_01·N / π
        print(f"\n  Velocity from Sz=0→1 gap (assuming x_spin=1/2):")
        v_from_gap = {}
        for N in sizes:
            v_g = gaps_sz01[N] * N / np.pi
            v_from_gap[N] = v_g

        # Extrapolate v from largest sizes using pairwise
        Ns = sorted(sizes)
        v_last = v_from_gap[Ns[-1]]
        v_second = v_from_gap[Ns[-2]]
        # Linear extrapolation in 1/N
        v_extrap = v_last + (v_last - v_second) * (1.0/Ns[-1]) / (1.0/Ns[-2] - 1.0/Ns[-1])

        if params['v_exact'] is not None:
            v_use = params['v_exact']
            print(f"  v_exact = {v_use:.6f}")
            for N in sorted(v_from_gap.keys()):
                print(f"    N={N:2d}: v_gap={v_from_gap[N]:.6f}  "
                      f"ratio={v_from_gap[N]/v_use:.4f}")
        else:
            v_use = v_last  # use largest-N estimate
            print(f"  v_gap(N={Ns[-1]}) = {v_last:.6f}  (used as v)")
            for N in sorted(v_from_gap.keys()):
                print(f"    N={N:2d}: v_gap={v_from_gap[N]:.6f}")

        # --- Multi-size fit: E0/N = a + b/N^2 + d/N^4 ---
        print(f"\n  Multi-size fit: E0/N = a + b/N^2 + d/N^4")
        x_data = np.array([1.0/N**2 for N in Ns])
        y_data = np.array([energies[N]/N for N in Ns])

        def casimir_model(x, a, b, d):
            return a + b * x + d * x**2

        popt, pcov = curve_fit(casimir_model, x_data, y_data)
        a_fit, b_fit, d_fit = popt
        vc_fit = -6.0 * b_fit / np.pi
        c_fit = vc_fit / v_use

        residuals = y_data - casimir_model(x_data, *popt)
        print(f"  eps_inf = {a_fit:.10f}")
        print(f"  vc_fit = {vc_fit:.6f}")
        print(f"  c_Cas(fit) = {c_fit:.6f}")
        print(f"  max residual = {np.max(np.abs(residuals)):.2e}")

        # --- Multi-size fit: S(N) = c/3 * ln(N/pi) + c' ---
        print(f"\n  Entropy fit: S = c/3 * ln(N/pi) + c'")
        x_ent = np.array([np.log(N / np.pi) for N in Ns])
        y_ent = np.array([entropies[N] for N in Ns])

        # Include 1/N^2 correction: S = c/3 * ln(N/pi) + c' + e/N^2
        def entropy_model(x_N, c_over3, cp, e):
            ln_N_pi, inv_N2 = x_N
            return c_over3 * ln_N_pi + cp + e * inv_N2

        x_ent_2d = np.vstack([x_ent, [1.0/N**2 for N in Ns]])
        popt_ent, _ = curve_fit(entropy_model, x_ent_2d, y_ent)
        c_eff_fit = 3.0 * popt_ent[0]
        residuals_ent = y_ent - entropy_model(x_ent_2d, *popt_ent)
        print(f"  c_eff(fit) = {c_eff_fit:.6f}")
        print(f"  max residual = {np.max(np.abs(residuals_ent)):.2e}")

        # --- Pairwise analysis ---
        print(f"\n  {'N1':>4} {'N2':>4}  {'c_eff':>8}  {'c_Cas':>8}  "
              f"{'|Δc_eff|':>10}  {'|Δc_Cas|':>10}  {'ratio':>6}")

        pw_ceff, pw_ccas = [], []
        for i in range(len(Ns)-1):
            N1, N2 = Ns[i], Ns[i+1]
            c_eff = 3.0 * (entropies[N2] - entropies[N1]) / np.log(N2/N1)
            eN1 = energies[N1] / N1
            eN2 = energies[N2] / N2
            vc = 6.0 * (eN2 - eN1) / (np.pi * (1.0/N1**2 - 1.0/N2**2))
            c_cas = vc / v_use

            pw_ceff.append((N1, N2, float(c_eff)))
            pw_ccas.append((N1, N2, float(c_cas)))

            dev_eff = abs(c_eff - 1.0)
            dev_cas = abs(c_cas - 1.0)
            ratio = dev_eff / dev_cas if dev_cas > 1e-8 else float('inf')
            print(f"  {N1:4d} {N2:4d}  {c_eff:8.5f}  {c_cas:8.5f}  "
                  f"{dev_eff:10.6f}  {dev_cas:10.6f}  {ratio:6.1f}×")

        # Convergence rate: |c_pair - c_last_pair| for last 3 pairs
        print(f"\n  Convergence (last-pair deviation from c=1):")
        c_eff_last = pw_ceff[-1][2]
        c_cas_last = pw_ccas[-1][2]
        dev_eff = abs(c_eff_last - 1.0)
        dev_cas = abs(c_cas_last - 1.0)
        ratio = dev_eff / dev_cas if dev_cas > 1e-8 else float('inf')
        print(f"  c_eff = {c_eff_last:.6f}  |Δ| = {dev_eff:.6f}")
        print(f"  c_Cas = {c_cas_last:.6f}  |Δ| = {dev_cas:.6f}")
        print(f"  ratio (eff/Cas) = {ratio:.2f}")

        results[name] = {
            'params': params,
            'sizes': Ns,
            'energies': {str(k): v for k, v in energies.items()},
            'entropies': {str(k): v for k, v in entropies.items()},
            'gaps_sz0': {str(k): v for k, v in gaps_sz0.items()},
            'gaps_sz01': {str(k): v for k, v in gaps_sz01.items()},
            'v_use': float(v_use),
            'v_from_gap': {str(k): v for k, v in v_from_gap.items()},
            'c_Cas_fit': float(c_fit),
            'c_eff_fit': float(c_eff_fit),
            'pairwise_ceff': pw_ceff,
            'pairwise_ccas': pw_ccas,
        }

    # --- FINAL SUMMARY ---
    print(f"\n{'='*70}")
    print(f"  FINAL SUMMARY")
    print(f"{'='*70}")
    print(f"\n  Fits (all sizes):")
    print(f"  {'Model':>15}  {'c_eff(fit)':>10}  {'c_Cas(fit)':>10}  {'v_used':>8}")
    for name in ['XX', 'Heisenberg', 'BKT']:
        r = results[name]
        print(f"  {name:>15}  {r['c_eff_fit']:10.5f}  {r['c_Cas_fit']:10.5f}  {r['v_use']:8.4f}")

    print(f"\n  Last-pair (20,22) deviations:")
    print(f"  {'Model':>15}  {'|Δc_eff|':>10}  {'|Δc_Cas|':>10}  {'winner':>10}  {'margin':>8}")
    for name in ['XX', 'Heisenberg', 'BKT']:
        r = results[name]
        ce = r['pairwise_ceff'][-1][2]
        cc = r['pairwise_ccas'][-1][2]
        de = abs(ce - 1.0)
        dc = abs(cc - 1.0)
        winner = "entropy" if de < dc else "Casimir"
        margin = dc/de if de > 0 and de < dc else de/dc if dc > 0 else 1.0
        print(f"  {name:>15}  {de:10.6f}  {dc:10.6f}  {winner:>10}  {margin:8.1f}×")

    # Save
    outpath = os.path.join(os.path.dirname(__file__), 'results', 'sprint_104b_extended.json')
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    with open(outpath, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nResults saved to {outpath}")

if __name__ == '__main__':
    run()
