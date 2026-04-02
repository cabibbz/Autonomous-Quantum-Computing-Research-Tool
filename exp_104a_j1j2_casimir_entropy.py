"""Sprint 104a: J1-J2 Heisenberg chain — Casimir energy vs entropy.

Compare c_Cas vs c_eff convergence to c=1 across three regimes:
  1. XX model (Delta=0): c=1, v=1, no log corrections (clean control)
  2. Heisenberg (Delta=1, J2=0): c=1, v=pi/2, mild log corrections (marginal operator)
  3. J1-J2 at BKT (J2c=0.2412): c=1, strong log^3 corrections

Uses S_z=0 sector, periodic BC, exact diag with GPU for N>=18.

Casimir formula: E0(N) = eps_inf*N - pi*v*c/(6N) + O(1/N^3)
=> E0(N)/N = eps_inf - pi*v*c/(6*N^2)
Pairwise: vc = 6*(E0(N2)/N2 - E0(N1)/N1) / (pi*(1/N1^2 - 1/N2^2))  [N2 > N1]

Entropy formula: S(N/2, N) = c/3 * ln(N/pi) + c'
Pairwise: c_eff = 3*(S(N2) - S(N1)) / ln(N2/N1)
"""
import numpy as np
from scipy.sparse import csr_matrix
from gpu_utils import eigsh
import json, time, os
from itertools import combinations
from db_utils import record

# --- S_z=0 sector basis ---
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

# --- Build XXZ + J2 Hamiltonian in S_z=0 sector (periodic BC) ---
def build_hamiltonian(N, Delta=1.0, J2=0.0):
    """H = sum_i [0.5*(S+_i S-_{i+1} + h.c.) + Delta*Sz_i*Sz_{i+1}]
         + J2 * sum_i [0.5*(S+_i S-_{i+2} + h.c.) + Delta*Sz_i*Sz_{i+2}]
    """
    basis = sz0_basis(N)
    dim = len(basis)
    idx = state_index(basis)

    rows, cols, vals = [], [], []

    for k, state in enumerate(basis):
        diag = 0.0

        # Build bond list: (site_i, site_j, coupling_xy, coupling_zz)
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

            # Sz_i Sz_j
            diag += Jzz * (si - 0.5) * (sj - 0.5)

            # S+_i S-_j: need si=0 (down->up), sj=1 (up->down)
            if si == 0 and sj == 1:
                new_state = (state | (1 << i)) & ~(1 << j)
                ni = idx.get(int(new_state))
                if ni is not None:
                    rows.append(k); cols.append(ni); vals.append(0.5 * Jxy)

            # S-_i S+_j: need si=1 (up->down), sj=0 (down->up)
            if si == 1 and sj == 0:
                new_state = (state & ~(1 << i)) | (1 << j)
                ni = idx.get(int(new_state))
                if ni is not None:
                    rows.append(k); cols.append(ni); vals.append(0.5 * Jxy)

        rows.append(k); cols.append(k); vals.append(diag)

    H = csr_matrix((vals, (rows, cols)), shape=(dim, dim))
    return H, basis

# --- Midchain entanglement entropy ---
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
    dim_eff = len(a_configs)

    b_set = set()
    for a in a_configs:
        for b, _ in a_groups[a]:
            b_set.add(b)
    b_configs = sorted(b_set)
    b_idx = {b: i for i, b in enumerate(b_configs)}
    dim_b = len(b_configs)

    M = np.zeros((dim_eff, dim_b), dtype=complex)
    for a in a_configs:
        ai = a_idx[a]
        for b, amp in a_groups[a]:
            M[ai, b_idx[b]] = amp

    rho_A = M @ M.conj().T
    eigvals = np.linalg.eigvalsh(rho_A)
    eigvals = eigvals[eigvals > 1e-15]
    return float(-np.sum(eigvals * np.log(eigvals)))

# --- Main ---
def run():
    models = {
        'XX': {'Delta': 0.0, 'J2': 0.0, 'v_exact': 1.0, 'label': 'XX (Δ=0)'},
        'Heisenberg': {'Delta': 1.0, 'J2': 0.0, 'v_exact': np.pi/2, 'label': 'Heisenberg (Δ=1)'},
        'BKT': {'Delta': 1.0, 'J2': 0.2412, 'v_exact': None, 'label': 'J1-J2 BKT (J2=0.2412)'},
    }

    sizes = [8, 10, 12, 14, 16, 18, 20]
    results = {}

    for name, params in models.items():
        print(f"\n{'='*60}")
        print(f"  {params['label']}  (Delta={params['Delta']}, J2={params['J2']})")
        print(f"{'='*60}")

        energies, entropies, gaps = {}, {}, {}

        for N in sizes:
            t0 = time.time()
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
            dt = time.time() - t0

            energies[N] = E0
            entropies[N] = S
            gaps[N] = gap

            print(f"  N={N:2d}  dim={dim:>8d}  E0/N={E0/N:.8f}  gap={gap:.6f}  "
                  f"S={S:.6f}  ({dt:.1f}s)")

            record(sprint=104, model='J1J2', q=0, n=N,
                   quantity='E0', value=E0, method=name)
            record(sprint=104, model='J1J2', q=0, n=N,
                   quantity='S_mid', value=S, method=name)

        # --- Pairwise extraction ---
        Ns = sorted(energies.keys())

        # Velocity: use exact if known, else extract from largest-pair vc
        v_exact = params['v_exact']

        pw_vc = []
        pw_ceff = []

        for i in range(len(Ns) - 1):
            N1, N2 = Ns[i], Ns[i+1]

            # c_eff from entropy
            c_eff = 3.0 * (entropies[N2] - entropies[N1]) / np.log(N2 / N1)

            # vc from Casimir: E0(N)/N = eps_inf - pi*v*c/(6*N^2)
            # vc = 6*(E0(N2)/N2 - E0(N1)/N1) / (pi * (1/N1^2 - 1/N2^2))
            eN1 = energies[N1] / N1
            eN2 = energies[N2] / N2
            denom = np.pi * (1.0/N1**2 - 1.0/N2**2)
            vc = 6.0 * (eN2 - eN1) / denom  # Note: eN2 > eN1 (less negative), denom > 0

            pw_vc.append((N1, N2, float(vc)))
            pw_ceff.append((N1, N2, float(c_eff)))

        # For BKT: estimate v from the converged vc (since c=1, v=vc)
        if v_exact is None:
            v_est = pw_vc[-1][2]  # use last pair as best estimate
            print(f"\n  v estimated from last pair: {v_est:.6f}")
            v_use = v_est
        else:
            v_use = v_exact
            print(f"\n  v (exact): {v_use:.6f}")

        print(f"\n  {'N1':>4} {'N2':>4}  {'c_eff':>8}  {'c_Cas':>8}  "
              f"{'|c_eff-1|':>10}  {'|c_Cas-1|':>10}  {'ratio':>6}")
        print(f"  {'-'*62}")

        pw_ccas = []
        for i in range(len(pw_vc)):
            N1, N2, vc = pw_vc[i]
            c_cas = vc / v_use
            c_eff = pw_ceff[i][2]

            dev_cas = abs(c_cas - 1.0)
            dev_eff = abs(c_eff - 1.0)
            ratio = dev_eff / dev_cas if dev_cas > 1e-10 else float('inf')

            pw_ccas.append((N1, N2, float(c_cas)))
            print(f"  {N1:4d} {N2:4d}  {c_eff:8.5f}  {c_cas:8.5f}  "
                  f"{dev_eff:10.6f}  {dev_cas:10.6f}  {ratio:6.1f}×")

        # Store
        results[name] = {
            'params': {k: v for k, v in params.items() if k != 'label'},
            'label': params['label'],
            'sizes': Ns,
            'energies': {str(k): v for k, v in energies.items()},
            'entropies': {str(k): v for k, v in entropies.items()},
            'gaps': {str(k): v for k, v in gaps.items()},
            'v_use': float(v_use),
            'pairwise_ceff': pw_ceff,
            'pairwise_ccas': pw_ccas,
            'pairwise_vc': pw_vc,
        }

    # --- Summary comparison ---
    print(f"\n{'='*70}")
    print(f"  SUMMARY: Last-pair deviations from c=1")
    print(f"{'='*70}")
    print(f"  {'Model':>20}  {'c_eff':>8}  {'c_Cas':>8}  "
          f"{'|Δc_eff|':>10}  {'|Δc_Cas|':>10}  {'eff/Cas':>8}")
    for name in ['XX', 'Heisenberg', 'BKT']:
        r = results[name]
        c_eff = r['pairwise_ceff'][-1][2]
        c_cas = r['pairwise_ccas'][-1][2]
        dev_eff = abs(c_eff - 1.0)
        dev_cas = abs(c_cas - 1.0)
        ratio = dev_eff / dev_cas if dev_cas > 1e-10 else float('inf')
        print(f"  {r['label']:>20}  {c_eff:8.5f}  {c_cas:8.5f}  "
              f"{dev_eff:10.6f}  {dev_cas:10.6f}  {ratio:8.1f}×")

    # Save
    outpath = os.path.join(os.path.dirname(__file__), 'results', 'sprint_104a_j1j2.json')
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    with open(outpath, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nResults saved to {outpath}")

if __name__ == '__main__':
    run()
