"""
Sprint 024a: Rotated Surface Code [[9,1,3]] — Entanglement Structure

Build the d=3 rotated surface code (9 data qubits) and compare its
entanglement structure to Shor [[9,1,3]] and [[5,1,3]].

Both surface and Shor are [[9,1,3]] — same n, same k, same d.
How does the entanglement topology differ?

Measures: pairwise MI, tripartite information (I3), negativity spectrum.
"""

import numpy as np
import json, time
from itertools import combinations

start = time.time()

I2 = np.eye(2)
X = np.array([[0, 1], [1, 0]])
Y = np.array([[0, -1j], [1j, 0]])
Z = np.array([[1, 0], [0, -1]])

def multi_pauli(ops_dict, n):
    result = np.eye(1)
    for q in range(n):
        result = np.kron(result, ops_dict.get(q, I2))
    return result

def partial_trace(rho, keep, n):
    """Trace out all qubits not in 'keep' list."""
    dim = 2**n
    keep = sorted(keep)
    trace_out = sorted(set(range(n)) - set(keep))
    rho_reshaped = rho.reshape([2]*2*n)
    # Contract traced-out indices
    for i, q in enumerate(sorted(trace_out, reverse=True)):
        rho_reshaped = np.trace(rho_reshaped, axis1=q, axis2=q+n-i)
    return rho_reshaped.reshape(2**len(keep), 2**len(keep))

def von_neumann_entropy(rho):
    evals = np.linalg.eigvalsh(rho)
    evals = evals[evals > 1e-12]
    return -np.sum(evals * np.log2(evals))

def mutual_info(rho, i, j, n):
    rho_i = partial_trace(rho, [i], n)
    rho_j = partial_trace(rho, [j], n)
    rho_ij = partial_trace(rho, [i, j], n)
    return von_neumann_entropy(rho_i) + von_neumann_entropy(rho_j) - von_neumann_entropy(rho_ij)

def tripartite_info(rho, i, j, k, n):
    S_i = von_neumann_entropy(partial_trace(rho, [i], n))
    S_j = von_neumann_entropy(partial_trace(rho, [j], n))
    S_k = von_neumann_entropy(partial_trace(rho, [k], n))
    S_ij = von_neumann_entropy(partial_trace(rho, [i, j], n))
    S_ik = von_neumann_entropy(partial_trace(rho, [i, k], n))
    S_jk = von_neumann_entropy(partial_trace(rho, [j, k], n))
    S_ijk = von_neumann_entropy(partial_trace(rho, [i, j, k], n))
    return S_i + S_j + S_k - S_ij - S_ik - S_jk + S_ijk

def negativity(rho, subsys_a, n):
    """Compute negativity across bipartition A|B."""
    keep = list(range(n))
    n_a = len(subsys_a)
    n_b = n - n_a
    dim_a = 2**n_a
    dim_b = 2**n_b
    # Reorder qubits: A first, then B
    subsys_b = sorted(set(range(n)) - set(subsys_a))
    order = sorted(subsys_a) + subsys_b
    # Permute rho
    rho_perm = rho.reshape([2]*2*n)
    axes = order + [o + n for o in order]
    rho_perm = rho_perm.transpose(axes)
    rho_perm = rho_perm.reshape(dim_a * dim_b, dim_a * dim_b)
    # Partial transpose over A
    rho_pt = rho_perm.reshape(dim_a, dim_b, dim_a, dim_b)
    rho_pt = rho_pt.transpose(2, 1, 0, 3)
    rho_pt = rho_pt.reshape(dim_a * dim_b, dim_a * dim_b)
    evals = np.linalg.eigvalsh(rho_pt)
    return float(np.sum(np.abs(evals[evals < -1e-12])))

def density_matrix(sv):
    sv = sv.reshape(-1, 1)
    return sv @ sv.conj().T

# ============================================================
# Build Rotated Surface Code [[9,1,3]]
# ============================================================
# Layout (d=3 rotated surface code):
#   0 - 1 - 2
#   |   |   |
#   3 - 4 - 5
#   |   |   |
#   6 - 7 - 8
#
# X stabilizers:
#   X0 X1          (weight-2, top boundary)
#   X1 X2 X4 X5    (weight-4)
#   X3 X4 X6 X7    (weight-4)
#   X7 X8          (weight-2, bottom boundary)
#
# Z stabilizers:
#   Z3 Z6          (weight-2, left boundary)
#   Z0 Z1 Z3 Z4    (weight-4)
#   Z4 Z5 Z7 Z8    (weight-4)
#   Z2 Z5          (weight-2, right boundary)
#
# Logical operators:
#   X_L = X0 X3 X6 (left column)
#   Z_L = Z0 Z1 Z2 (top row)

n_surf = 9
dim_surf = 2**n_surf

# Stabilizer generators
x_stabs = [
    {0: X, 1: X},                     # weight-2 top
    {1: X, 2: X, 4: X, 5: X},         # weight-4
    {3: X, 4: X, 6: X, 7: X},         # weight-4
    {7: X, 8: X},                     # weight-2 bottom
]
z_stabs = [
    {3: Z, 6: Z},                     # weight-2 left
    {0: Z, 1: Z, 3: Z, 4: Z},         # weight-4
    {4: Z, 5: Z, 7: Z, 8: Z},         # weight-4
    {2: Z, 5: Z},                     # weight-2 right
]

# Logical operators
X_L_surf = multi_pauli({0: X, 3: X, 6: X}, n_surf)
Z_L_surf = multi_pauli({0: Z, 1: Z, 2: Z}, n_surf)

# Build projector onto code space
all_stabs = x_stabs + z_stabs
print("Building surface code projector...")
t0 = time.time()
proj = np.eye(dim_surf, dtype=complex)
for s in all_stabs:
    S = multi_pauli(s, n_surf)
    proj = proj @ (np.eye(dim_surf) + S) / 2
print(f"  Projector built in {time.time()-t0:.2f}s, rank={np.linalg.matrix_rank(proj, tol=1e-8)}")

# Extract logical |0> and |1>
# |0_L>: +1 eigenstate of Z_L within code space
# Project |00...0> into code space
state_0 = np.zeros(dim_surf)
state_0[0] = 1.0
sv_0z_surf = proj @ state_0
norm = np.linalg.norm(sv_0z_surf)
if norm < 1e-10:
    # Try a different seed state
    for seed in range(dim_surf):
        state_s = np.zeros(dim_surf)
        state_s[seed] = 1.0
        sv_0z_surf = proj @ state_s
        norm = np.linalg.norm(sv_0z_surf)
        if norm > 1e-10:
            print(f"  Used seed state |{seed:09b}> (norm={norm:.4f})")
            break

sv_0z_surf /= np.linalg.norm(sv_0z_surf)

# Check Z_L eigenvalue
zl_expect = np.real(sv_0z_surf.conj() @ Z_L_surf @ sv_0z_surf)
print(f"  Z_L expectation of projected state: {zl_expect:.4f}")

if zl_expect < 0:
    # This is |1_L>, flip
    sv_1z_surf = sv_0z_surf.copy()
    sv_0z_surf = X_L_surf @ sv_1z_surf
    sv_0z_surf /= np.linalg.norm(sv_0z_surf)
else:
    sv_1z_surf = X_L_surf @ sv_0z_surf
    sv_1z_surf /= np.linalg.norm(sv_1z_surf)

# Verify
assert abs(np.dot(sv_0z_surf.conj(), sv_1z_surf)) < 1e-10, "Logical states not orthogonal!"
print(f"  Z_L|0_L> = {np.real(sv_0z_surf.conj() @ Z_L_surf @ sv_0z_surf):.4f}")
print(f"  Z_L|1_L> = {np.real(sv_1z_surf.conj() @ Z_L_surf @ sv_1z_surf):.4f}")

# ============================================================
# Build Shor [[9,1,3]] for comparison
# ============================================================
block_plus = np.zeros(8)
block_plus[0b000] = 1.0; block_plus[0b111] = 1.0
block_plus /= np.linalg.norm(block_plus)

block_minus = np.zeros(8)
block_minus[0b000] = 1.0; block_minus[0b111] = -1.0
block_minus /= np.linalg.norm(block_minus)

sv_0z_shor = np.kron(np.kron(block_plus, block_plus), block_plus)
sv_1z_shor = np.kron(np.kron(block_minus, block_minus), block_minus)

# ============================================================
# Build [[5,1,3]] for comparison
# ============================================================
n5 = 5; dim5 = 2**n5
stabs_513 = [
    {0: X, 1: Z, 2: Z, 3: X},
    {1: X, 2: Z, 3: Z, 4: X},
    {0: X, 2: X, 3: Z, 4: Z},
    {0: Z, 1: X, 3: X, 4: Z},
]
proj_513 = np.eye(dim5)
for s in stabs_513:
    S = multi_pauli(s, n5)
    proj_513 = proj_513 @ (np.eye(dim5) + S) / 2
state5_0 = np.zeros(dim5); state5_0[0] = 1.0
sv_0z_513 = proj_513 @ state5_0
sv_0z_513 /= np.linalg.norm(sv_0z_513)
X_L_513 = multi_pauli({q: X for q in range(5)}, 5)
sv_1z_513 = X_L_513 @ sv_0z_513
sv_1z_513 /= np.linalg.norm(sv_1z_513)

print(f"\nAll code states built in {time.time()-start:.1f}s")

# ============================================================
# Compute entanglement structure for each code
# ============================================================

codes = {
    'surface-9': {'n': 9, 'sv0': sv_0z_surf, 'sv1': sv_1z_surf},
    'shor-9':    {'n': 9, 'sv0': sv_0z_shor, 'sv1': sv_1z_shor},
    '[[5,1,3]]': {'n': 5, 'sv0': sv_0z_513,  'sv1': sv_1z_513},
}

results = {}

for code_name, code_data in codes.items():
    n = code_data['n']
    # Use |0_L> for entanglement structure (code is symmetric)
    rho = density_matrix(code_data['sv0'])

    print(f"\n=== {code_name} (n={n}) ===")

    # Time test for partial_trace
    t0 = time.time()
    _ = partial_trace(rho, [0], n)
    t_pt = time.time() - t0
    print(f"  partial_trace timing: {t_pt:.3f}s")

    # Pairwise MI
    print("  Computing pairwise MI...")
    mi_pairs = {}
    total_mi = 0
    for i, j in combinations(range(n), 2):
        mi = mutual_info(rho, i, j, n)
        mi_pairs[f"({i},{j})"] = round(float(mi), 4)
        total_mi += mi

    nonzero_mi = {k: v for k, v in mi_pairs.items() if abs(v) > 0.01}
    print(f"  Total MI: {total_mi:.4f}")
    print(f"  Non-zero MI pairs ({len(nonzero_mi)}/{len(mi_pairs)}): {nonzero_mi}")

    # I3 for all triples (only for n <= 9)
    print("  Computing tripartite information...")
    i3_triples = {}
    n_neg = 0
    n_pos = 0
    n_zero = 0
    i3_min = float('inf')
    i3_max = float('-inf')

    t0 = time.time()
    triple_list = list(combinations(range(n), 3))
    # Time test on first triple
    if len(triple_list) > 0:
        i3_test = tripartite_info(rho, *triple_list[0], n)
        t_i3 = time.time() - t0
        print(f"  I3 timing per triple: {t_i3:.3f}s, total triples: {len(triple_list)}")
        est = t_i3 * len(triple_list)
        if est > 45:
            print(f"  WARNING: I3 would take ~{est:.0f}s, skipping")
            i3_triples = "SKIPPED"
        else:
            for triple in triple_list:
                i3 = tripartite_info(rho, *triple, n)
                key = str(triple)
                i3_triples[key] = round(float(i3), 4)
                if i3 < -0.01: n_neg += 1
                elif i3 > 0.01: n_pos += 1
                else: n_zero += 1
                i3_min = min(i3_min, i3)
                i3_max = max(i3_max, i3)

            print(f"  I3: neg={n_neg}, pos={n_pos}, zero={n_zero}")
            print(f"  I3 range: [{i3_min:.4f}, {i3_max:.4f}]")

    # Negativity spectrum (all bipartitions of size 1,2,3,4)
    print("  Computing negativity spectrum...")
    neg_spectrum = {}
    for size_a in range(1, n//2 + 1):
        negs = []
        for sub_a in combinations(range(n), size_a):
            neg = negativity(rho, list(sub_a), n)
            negs.append(neg)
        neg_spectrum[f"|A|={size_a}"] = {
            'min': round(min(negs), 4),
            'max': round(max(negs), 4),
            'mean': round(float(np.mean(negs)), 4),
            'spread': round(max(negs) - min(negs), 4),
        }
        print(f"  Neg |A|={size_a}: min={min(negs):.4f} max={max(negs):.4f} mean={np.mean(negs):.4f}")

    results[code_name] = {
        'n': n,
        'total_mi': round(float(total_mi), 4),
        'mi_pairs': mi_pairs,
        'nonzero_mi_count': len(nonzero_mi),
        'i3_triples': i3_triples,
        'i3_stats': {
            'negative': n_neg, 'positive': n_pos, 'zero': n_zero,
            'min': round(float(i3_min), 4) if i3_min != float('inf') else None,
            'max': round(float(i3_max), 4) if i3_max != float('-inf') else None,
        } if i3_triples != "SKIPPED" else "SKIPPED",
        'negativity_spectrum': neg_spectrum,
    }

# ============================================================
# Comparative summary
# ============================================================
print("\n" + "="*60)
print("COMPARATIVE SUMMARY")
print("="*60)
print(f"{'Metric':30s} {'Surface-9':>12s} {'Shor-9':>12s} {'[[5,1,3]]':>12s}")
print("-"*66)
print(f"{'Total pairwise MI':30s} {results['surface-9']['total_mi']:>12.4f} {results['shor-9']['total_mi']:>12.4f} {results['[[5,1,3]]']['total_mi']:>12.4f}")
print(f"{'Non-zero MI pairs':30s} {results['surface-9']['nonzero_mi_count']:>12d} {results['shor-9']['nonzero_mi_count']:>12d} {results['[[5,1,3]]']['nonzero_mi_count']:>12d}")

for code in ['surface-9', 'shor-9', '[[5,1,3]]']:
    i3s = results[code]['i3_stats']
    if i3s != "SKIPPED":
        print(f"  {code} I3: {i3s['negative']} neg, {i3s['positive']} pos, {i3s['zero']} zero, range [{i3s['min']}, {i3s['max']}]")

for size_key in ['|A|=1', '|A|=2']:
    for code in ['surface-9', 'shor-9', '[[5,1,3]]']:
        ns = results[code]['negativity_spectrum'].get(size_key, {})
        if ns:
            print(f"  {code} Neg {size_key}: {ns['min']:.4f}-{ns['max']:.4f} (spread={ns['spread']:.4f})")

elapsed = time.time() - start
print(f"\nTotal time: {elapsed:.1f}s")

# Save results
output = {
    'sprint': '024a',
    'description': 'Surface code [[9,1,3]] entanglement structure vs Shor and [[5,1,3]]',
    'results': results,
    'elapsed_seconds': round(elapsed, 1),
}

with open('results/sprint_024a_surface_code_structure.json', 'w') as f:
    json.dump(output, f, indent=2)
print("Results saved to results/sprint_024a_surface_code_structure.json")
