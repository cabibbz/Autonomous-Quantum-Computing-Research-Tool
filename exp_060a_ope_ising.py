"""Sprint 060a: OPE coefficient extraction validated on q=2 Ising CFT.

Method: On a periodic chain of N sites, CFT eigenstates |a> correspond to
primary operators. The matrix element of a local operator phi(0) between
eigenstates gives the OPE coefficient:

  <a|phi_local(0)|b> = C_{a,phi,b} * (2*pi/N)^{x_phi}

For the Ising model:
  - |0> = identity (ground state)
  - |sigma> = spin field (1st excited, momentum 0)
  - |epsilon> = energy field

Local operators:
  - sigma_local = sigma_z at site 0 (spin field)
  - epsilon_local = sigma_z(0)*sigma_z(1) (energy density)

Known exact OPE coefficients (Ising CFT):
  C_{sigma,sigma,epsilon} = 1/2
  C_{sigma,sigma,I} = 1 (normalization)

We extract C from: C = <a|O_local(0)|b> * N^{x_O} / (2*pi)^{x_O}
with appropriate normalization.
"""
import numpy as np
from scipy.sparse import csr_matrix, kron as sp_kron, eye as sp_eye
from scipy.sparse.linalg import eigsh
import json, time

def potts_hamiltonian_periodic(n, q, g):
    """Build H = -sum delta(s_i, s_j) - g * sum (X + X^dag) for q-state Potts."""
    dim = q**n
    eye_q = np.eye(q)
    X = np.zeros((q, q))
    for s in range(q):
        X[(s+1) % q, s] = 1.0
    XpXd = X + X.T
    delta_2 = np.zeros((q**2, q**2))
    for s in range(q):
        delta_2[s*q+s, s*q+s] = 1.0

    H = csr_matrix((dim, dim))
    for i in range(n):
        j = (i + 1) % n
        # Two-site delta interaction
        if j == i + 1:
            left_dim = q**i if i > 0 else 1
            right_dim = q**(n-j-1) if n-j-1 > 0 else 1
            op = csr_matrix(delta_2)
            if i > 0:
                op = sp_kron(csr_matrix(np.eye(left_dim)), op, format='csr')
            if n-j-1 > 0:
                op = sp_kron(op, csr_matrix(np.eye(right_dim)), format='csr')
            H = H - op
        else:
            # Wrap-around bond (i=n-1, j=0)
            for s in range(q):
                proj = np.zeros((q, q)); proj[s, s] = 1.0
                ops = [eye_q] * n; ops[0] = proj; ops[n-1] = proj
                H = H - tensor_op(ops, q)
    for i in range(n):
        ops = [eye_q] * n; ops[i] = XpXd
        H = H - g * tensor_op(ops, q)
    return H

def tensor_op(ops_list, q):
    result = csr_matrix(ops_list[0])
    for op in ops_list[1:]:
        result = sp_kron(result, csr_matrix(op), format='csr')
    return result

def build_translation_operator(n, q):
    dim = q**n
    rows, cols = [], []
    for idx in range(dim):
        state = []
        tmp = idx
        for _ in range(n):
            state.append(tmp % q)
            tmp //= q
        new_state = state[1:] + [state[0]]
        new_idx = sum(s * q**k for k, s in enumerate(new_state))
        rows.append(new_idx)
        cols.append(idx)
    return csr_matrix((np.ones(dim), (rows, cols)), shape=(dim, dim))

def build_local_clock_op(n, q, site=0):
    """Build Z = diag(1, omega, omega^2, ...) at given site. This is the
    local spin field operator (order parameter) for the q-state Potts model."""
    omega = np.exp(2j * np.pi / q)
    Z = np.diag([omega**s for s in range(q)])
    eye_q = np.eye(q)
    ops = [eye_q] * n
    ops[site] = Z
    return tensor_op(ops, q)

def build_local_energy_op(n, q, site=0):
    """Build the local energy density operator: delta(s_site, s_{site+1}) at one bond."""
    eye_q = np.eye(q)
    delta_2 = np.zeros((q**2, q**2))
    for s in range(q):
        delta_2[s*q+s, s*q+s] = 1.0

    j = (site + 1) % n
    if j == site + 1:
        left_dim = q**site if site > 0 else 1
        right_dim = q**(n-j-1) if n-j-1 > 0 else 1
        op = csr_matrix(delta_2)
        if site > 0:
            op = sp_kron(csr_matrix(np.eye(left_dim)), op, format='csr')
        if n-j-1 > 0:
            op = sp_kron(op, csr_matrix(np.eye(right_dim)), format='csr')
        return op
    else:
        # Wrap-around
        result = csr_matrix((q**n, q**n))
        for s in range(q):
            proj = np.zeros((q, q)); proj[s, s] = 1.0
            ops = [eye_q] * n; ops[0] = proj; ops[n-1] = proj
            result = result + tensor_op(ops, q)
        return result

def get_momentum_states(evals, evecs, T, n):
    """Assign momentum quantum numbers to eigenstates."""
    n_states = len(evals)
    momenta = np.zeros(n_states, dtype=int)
    tol = 1e-8 * abs(evals[-1] - evals[0]) if len(evals) > 1 else 1e-10

    i = 0
    while i < n_states:
        group = [i]
        j = i + 1
        while j < n_states and abs(evals[j] - evals[i]) < tol:
            group.append(j)
            j += 1
        sub = evecs[:, group]
        T_sub = (sub.conj().T) @ (T @ sub)
        t_evals = np.linalg.eigvals(T_sub)
        for idx_g, g_idx in enumerate(group):
            phase = np.angle(t_evals[idx_g])
            k = phase / (2 * np.pi / n)
            momenta[g_idx] = int(round(k)) % n
        i = j
    return momenta

def extract_ope(n_sites, q, g_c, x_sigma, x_epsilon, n_eig=30):
    """Extract OPE coefficients from finite-size matrix elements.

    Method: In CFT on cylinder of circumference L=N:
      <a|phi(0)|b> = C_{a,phi,b} * (2*pi/N)^{Delta_phi}

    For normalized CFT states, the OPE coefficient is:
      C_{sigma,sigma,epsilon} = <epsilon|sigma_local|sigma> * N / <0|sigma_local|sigma>

    We use the ratio method to cancel unknown normalizations.
    """
    t0 = time.time()
    dim = q**n_sites
    print(f"\n  q={q}, n={n_sites} (dim={dim})...", flush=True)

    H = potts_hamiltonian_periodic(n_sites, q, g_c)
    T = build_translation_operator(n_sites, q)

    n_eig = min(n_eig, dim - 2)
    evals, evecs = eigsh(H, k=n_eig, which='SA')
    order = np.argsort(evals)
    evals = evals[order]
    evecs = evecs[:, order]

    momenta = get_momentum_states(evals, evecs, T, n_sites)
    gaps = evals - evals[0]
    delta1 = gaps[1]
    ratios = gaps / delta1 if delta1 > 0 else gaps

    dt = time.time() - t0
    print(f"  Diag: {dt:.1f}s. E0={evals[0]:.6f}, Delta1={delta1:.6f}")

    # Build local clock operator Z at site 0
    Z_local = build_local_clock_op(n_sites, q, site=0)

    # Build local energy density at bond (0,1)
    E_local = build_local_energy_op(n_sites, q, site=0)

    # Identify states by momentum and ratio
    # Ground state: index 0, k=0
    # Sigma field: first k=0 excited state (R~1)
    # Epsilon field: look for R matching known ratio

    print(f"\n  Spectrum (first 15 levels):")
    print(f"  {'i':>3} {'R':>8} {'k':>3} {'gap':>10}")
    for i in range(min(15, len(ratios))):
        k = int(momenta[i])
        s = k if k <= n_sites//2 else k - n_sites
        print(f"  {i:>3} {ratios[i]:>8.4f} {s:>3} {gaps[i]:>10.6f}")

    # Find sigma: first k=0 state above ground (R near 1)
    sigma_indices = []
    for i in range(1, len(ratios)):
        k = int(momenta[i])
        s = k if k <= n_sites//2 else k - n_sites
        if s == 0 and 0.5 < ratios[i] < 1.5:
            sigma_indices.append(i)

    if not sigma_indices:
        print("  WARNING: No sigma state found!")
        return None

    # For q=2 Ising: sigma is at k=0 (but actually can be at k=pi for some conventions)
    # Let's also check non-zero momentum
    all_low = [(i, ratios[i], int(momenta[i])) for i in range(1, min(10, len(ratios)))]
    print(f"\n  Low-lying states: {all_low}")

    # Sigma state(s): the first excited level (may be degenerate)
    sigma_idx = 1  # First excited state

    # Epsilon: for Ising, R_epsilon/R_sigma = 8 → look for R~8
    # General: use known R_epsilon from Sprint 057
    R_epsilon_expected = {2: 8.0, 3: 6.0, 4: 6.6, 5: 7.1, 7: 8.1, 10: 8.3}
    R_eps = R_epsilon_expected.get(q, 6.0)

    # Find epsilon: k=0 state near R_eps
    epsilon_candidates = []
    for i in range(1, len(ratios)):
        k = int(momenta[i])
        s = k if k <= n_sites//2 else k - n_sites
        if s == 0 and abs(ratios[i] - R_eps) < 2.0:
            epsilon_candidates.append((i, ratios[i]))

    print(f"\n  Sigma at index {sigma_idx}, R={ratios[sigma_idx]:.4f}, k={int(momenta[sigma_idx])}")
    print(f"  Epsilon candidates (k=0, R near {R_eps}): {epsilon_candidates}")

    # Compute matrix elements of Z_local (clock operator = spin field)
    # <0|Z|sigma>, <epsilon|Z|sigma>, <0|Z|0>

    psi_0 = evecs[:, 0]  # Ground state
    psi_sigma = evecs[:, sigma_idx]

    # Z|sigma>
    Z_sigma = Z_local @ psi_sigma

    # <0|Z|sigma>
    me_0_Z_sigma = np.dot(psi_0.conj(), Z_sigma)
    print(f"\n  <0|Z|sigma> = {me_0_Z_sigma:.6f} (complex: {me_0_Z_sigma})")

    results = {
        "q": q, "n": n_sites, "dim": dim, "g_c": g_c,
        "E0": float(evals[0]),
        "delta1": float(delta1),
        "sigma_idx": sigma_idx,
        "R_sigma": float(ratios[sigma_idx]),
        "k_sigma": int(momenta[sigma_idx]),
        "me_0_Z_sigma": {"real": float(np.real(me_0_Z_sigma)),
                          "imag": float(np.imag(me_0_Z_sigma)),
                          "abs": float(np.abs(me_0_Z_sigma))},
        "epsilon_candidates": [],
        "ope_coefficients": {},
        "time_s": round(time.time() - t0, 2)
    }

    # For each epsilon candidate, compute C_{sigma,sigma,epsilon}
    # The OPE ratio method:
    #   C_{sigma,sigma,epsilon} / C_{I,sigma,sigma} = <epsilon|Z|sigma> / <0|Z|sigma> * N^{x_eps - x_sigma}...
    #
    # Actually, the clean formula from Cardy (1986):
    # On cylinder of circumference L, with states normalized <a|a>=1:
    #   <a|phi(0)|b> = C_{a,phi,b} * (2*pi/L)^{x_phi}  [for phi a primary at origin]
    #
    # But phi(0) is summed over the chain for translation eigenstates.
    # For a LOCAL operator at site 0:
    #   <a|phi_0|b> = C_{a,phi,b} * (2*pi/N)^{x_phi} * (1/sqrt(N)) for momentum eigenstates
    #
    # More carefully: if |sigma> is a momentum eigenstate with k=0, then
    #   <epsilon|Z_0|sigma> ~ C_{eps,sigma,sigma} * N^{-x_Z}
    # where Z_0 is the LOCAL clock operator at site 0.
    #
    # The RATIO cancels normalization issues:
    #   C_{sigma,sigma,epsilon} / C_{I,sigma,sigma} = |<epsilon|Z_0|sigma>| / |<0|Z_0|sigma>| * N^{x_epsilon - 0}
    # ... but this doesn't simplify nicely.
    #
    # Simplest reliable method: the DERIVATIVE approach.
    # dE_n/dg|_{g_c} = -<n|V|n> where V = sum_i (X+X^dag)_i (extensive)
    # The mixed matrix element <m|V_local|n> for V_local = (X+X^dag) at site 0:
    #   <m|V_0|n> = C_{m,epsilon,n} * (2*pi/N)^{x_epsilon}
    # because V_local is proportional to the energy density operator near g_c.
    #
    # Let's use BOTH methods and compare.

    # Method 1: Clock operator Z matrix elements
    # The properly normalized OPE coefficient from matrix elements:
    # C_{sigma,sigma,epsilon}^2 = |<epsilon|Z_0|sigma>|^2 / |<0|Z_0|sigma>|^2
    #                              * (ratios[eps]/1)^{2*x_sigma} ... complicated
    #
    # Method 2: Use the standard Cardy-Calabrese normalization.
    # For a primary phi with scaling dim x_phi on a periodic chain of N sites:
    #   <0|phi_0|phi> = A_phi * N^{-x_phi}
    # where A_phi = (2*pi)^{x_phi} is the normalization constant.
    # Then:
    #   C_{phi,phi,psi} = <psi|phi_0|phi> / <0|phi_0|phi> * (2*pi/N)^{x_psi - x_phi}
    #
    # For Ising: x_sigma = 1/8, x_epsilon = 1
    # C_{sigma,sigma,epsilon} = <epsilon|Z_0|sigma> / <0|Z_0|sigma> * (2*pi/N)^{1 - 1/8}
    #
    # Actually, the simplest and most reliable: use the RATIO of matrix elements
    # with appropriate scaling powers, and check N-independence.

    # Let me just compute ALL matrix elements and fit the scaling
    print(f"\n  === Matrix elements of Z_local (clock operator) ===")
    me_Z = np.zeros((min(15, n_eig), min(15, n_eig)), dtype=complex)
    for a in range(min(15, n_eig)):
        for b in range(min(15, n_eig)):
            me_Z[a, b] = np.dot(evecs[:, a].conj(), Z_local @ evecs[:, b])

    print(f"  {'a':>3} {'b':>3} {'|<a|Z|b>|':>12} {'R_a':>8} {'R_b':>8} {'k_a':>4} {'k_b':>4}")
    significant = []
    for a in range(min(15, n_eig)):
        for b in range(min(15, n_eig)):
            val = np.abs(me_Z[a, b])
            if val > 1e-6:
                ka = int(momenta[a]); sa = ka if ka <= n_sites//2 else ka - n_sites
                kb = int(momenta[b]); sb = kb if kb <= n_sites//2 else kb - n_sites
                significant.append((a, b, val, ratios[a], ratios[b], sa, sb))
                if a < 8 and b < 8:
                    print(f"  {a:>3} {b:>3} {val:>12.6f} {ratios[a]:>8.4f} {ratios[b]:>8.4f} {sa:>4} {sb:>4}")

    results["matrix_elements_Z"] = [
        {"a": int(a), "b": int(b), "abs_me": round(float(v), 8),
         "R_a": round(float(ra), 6), "R_b": round(float(rb), 6),
         "k_a": int(ka), "k_b": int(kb)}
        for a, b, v, ra, rb, ka, kb in significant if a < 15 and b < 15
    ]

    # Also compute matrix elements of the energy density operator
    print(f"\n  === Matrix elements of E_local (energy density) ===")
    me_E = np.zeros((min(10, n_eig), min(10, n_eig)), dtype=complex)
    for a in range(min(10, n_eig)):
        for b in range(min(10, n_eig)):
            me_E[a, b] = np.dot(evecs[:, a].conj(), E_local @ evecs[:, b])

    for a in range(min(10, n_eig)):
        for b in range(a, min(10, n_eig)):
            val = np.abs(me_E[a, b] - (me_E[a, b].real if abs(me_E[a,b].imag) < 1e-10 else me_E[a,b]))
            val_real = me_E[a, b].real
            val_abs = np.abs(me_E[a, b])
            if val_abs > 1e-6 and a != b:
                ka = int(momenta[a]); sa = ka if ka <= n_sites//2 else ka - n_sites
                kb = int(momenta[b]); sb = kb if kb <= n_sites//2 else kb - n_sites
                print(f"  {a:>3} {b:>3} {val_abs:>12.6f} {ratios[a]:>8.4f} {ratios[b]:>8.4f} {sa:>4} {sb:>4}")

    return results


# ===== ISING (q=2) VALIDATION =====
print("=" * 70)
print("SPRINT 060a: OPE COEFFICIENT EXTRACTION — q=2 ISING VALIDATION")
print("=" * 70)
print("Known exact: C_{sigma,sigma,epsilon} = 1/2")
print("x_sigma = 1/8, x_epsilon = 1, c = 1/2")
print()

all_results = {}
# q=2: can go to n=12 (dim=4096)
for n in [6, 8, 10, 12]:
    result = extract_ope(n, q=2, g_c=0.250, x_sigma=0.125, x_epsilon=1.0)
    if result:
        all_results[f"n={n}"] = result

# ===== SCALING ANALYSIS =====
print("\n\n" + "=" * 70)
print("SCALING ANALYSIS: Extracting C_{sigma,sigma,epsilon}")
print("=" * 70)

# For Ising, sigma is at k=0 or k=pi/a depending on boundary conditions.
# On our periodic chain with Z_2 symmetry:
#   - Ground state is in even sector
#   - sigma field connects even/odd sectors
# The clock operator Z = sigma_z maps between Z_2 sectors.
# So <0|Z|sigma> is nonzero only if |sigma> is in the odd sector.
#
# The OPE coefficient extraction:
# |<0|Z_0|sigma>|^2 = (2*pi)^{2*x_sigma} * N^{-2*x_sigma} * |C_{I,sigma,sigma}|^2  (=1 by normalization)
# So: |<0|Z_0|sigma>| = (2*pi)^{x_sigma} * N^{-x_sigma}
# Check: |<0|Z_0|sigma>| * N^{x_sigma} should be constant = (2*pi)^{x_sigma}

print("\nChecking <0|Z_0|sigma> scaling:")
for key, res in all_results.items():
    n = res["n"]
    me = res["me_0_Z_sigma"]["abs"]
    scaled = me * n**0.125  # x_sigma = 1/8
    print(f"  {key}: |<0|Z_0|sigma>| = {me:.6f}, * N^(1/8) = {scaled:.6f}")

# For epsilon:
# |<epsilon|Z_0|sigma>| = (2*pi)^{x_sigma} * N^{-x_sigma} * C_{sigma,sigma,epsilon}
# So: C_{sigma,sigma,epsilon} = |<epsilon|Z_0|sigma>| / |<0|Z_0|sigma>|
#
# This is the CLEAN ratio that doesn't need any scaling powers!
# C_{sigma,sigma,epsilon} = |<epsilon|Z_0|sigma>| / |<0|Z_0|sigma>|

print("\n\nC_{sigma,sigma,epsilon} from ratio method:")
print("C = |<epsilon|Z_0|sigma>| / |<0|Z_0|sigma>|")
print("(Exact Ising value: 0.5)")

for key, res in all_results.items():
    me_list = res.get("matrix_elements_Z", [])
    n = res["n"]
    # Find <0|Z|sigma> (a=0, b=sigma_idx)
    me_0_sigma = res["me_0_Z_sigma"]["abs"]

    # Find epsilon state: k=0, R near 8
    # <epsilon|Z|sigma>: look for matrix element with R_a near 8
    for me in me_list:
        if me["k_a"] == 0 and 6 < me["R_a"] < 10 and me["b"] == res["sigma_idx"]:
            C_ratio = me["abs_me"] / me_0_sigma if me_0_sigma > 0 else 0
            print(f"  {key}: <eps(R={me['R_a']:.2f})|Z|sigma> = {me['abs_me']:.6f}, "
                  f"C_ratio = {C_ratio:.4f}")

with open("results/sprint_060a_ope_ising.json", "w") as f:
    json.dump(all_results, f, indent=2, default=str)
print("\nResults saved to results/sprint_060a_ope_ising.json")
