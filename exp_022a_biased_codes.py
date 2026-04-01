"""
Sprint 022a: Specialized (bias-aware) codes vs isotropic [[5,1,3]] under biased noise

Key question: Under Z-biased noise (dephasing >> relaxation, as in real hardware),
does a phase-flip repetition code (specialized for Z errors) outperform [[5,1,3]]?

Codes:
- 3-qubit bit-flip code: |0>_L=|000>, |1>_L=|111>, stabilizers XXI, IXX (protects X errors)
- 3-qubit phase-flip code: |0>_L=|+++>, |1>_L=|--->, stabilizers ZZI, IZZ (protects Z errors)
- [[5,1,3]] perfect code: protects all single-qubit errors (isotropic)
- Uncoded: single qubit baseline

Sweep: fix gamma=0.05 (mild T1), sweep lambda from 0 to 0.5 (increasing T2 dominance).
Test all three logical bases (Z, X, Y) for fair comparison.
"""

import numpy as np
import json, time

start = time.time()

I2 = np.eye(2)
X = np.array([[0, 1], [1, 0]])
Y = np.array([[0, -1j], [1j, 0]])
Z = np.array([[1, 0], [0, -1]])
H_gate = np.array([[1, 1], [1, -1]]) / np.sqrt(2)

def density_matrix(sv):
    sv = sv.reshape(-1, 1)
    return sv @ sv.conj().T

def von_neumann_entropy(rho):
    evals = np.linalg.eigvalsh(rho)
    evals = evals[evals > 1e-12]
    return -np.sum(evals * np.log2(evals))

def multi_pauli(ops_dict, n):
    result = np.eye(1)
    for q in range(n):
        result = np.kron(result, ops_dict.get(q, I2))
    return result

def apply_amplitude_damping(rho, gamma, n):
    K0 = np.array([[1, 0], [0, np.sqrt(1 - gamma)]])
    K1 = np.array([[0, np.sqrt(gamma)], [0, 0]])
    result = rho.copy().astype(complex)
    for q in range(n):
        rho_new = np.zeros_like(result)
        for K in [K0, K1]:
            K_full = np.eye(1, dtype=complex)
            for i in range(n):
                K_full = np.kron(K_full, K if i == q else I2)
            rho_new += K_full @ result @ K_full.conj().T
        result = rho_new
    return result

def apply_phase_damping(rho, lam, n):
    K0 = np.array([[1, 0], [0, np.sqrt(1 - lam)]])
    K1 = np.array([[0, 0], [0, np.sqrt(lam)]])
    result = rho.copy().astype(complex)
    for q in range(n):
        rho_new = np.zeros_like(result)
        for K in [K0, K1]:
            K_full = np.eye(1, dtype=complex)
            for i in range(n):
                K_full = np.kron(K_full, K if i == q else I2)
            rho_new += K_full @ result @ K_full.conj().T
        result = rho_new
    return result

def apply_combined_noise(rho, gamma, lam, n):
    rho = apply_amplitude_damping(rho, gamma, n)
    rho = apply_phase_damping(rho, lam, n)
    return rho

def holevo_info(rho_0, rho_1):
    rho_avg = (rho_0 + rho_1) / 2
    return von_neumann_entropy(rho_avg) - (von_neumann_entropy(rho_0) + von_neumann_entropy(rho_1)) / 2

# ============================================================
# Build code states
# ============================================================

# --- 3-qubit bit-flip code ---
n3 = 3
dim3 = 2**n3
sv_0z_bf = np.zeros(dim3); sv_0z_bf[0] = 1.0   # |000>
sv_1z_bf = np.zeros(dim3); sv_1z_bf[7] = 1.0   # |111>
sv_0x_bf = (sv_0z_bf + sv_1z_bf) / np.sqrt(2)
sv_1x_bf = (sv_0z_bf - sv_1z_bf) / np.sqrt(2)
sv_0y_bf = (sv_0z_bf + 1j * sv_1z_bf) / np.sqrt(2)
sv_1y_bf = (sv_0z_bf - 1j * sv_1z_bf) / np.sqrt(2)

# --- 3-qubit phase-flip code ---
# |0>_L = |+++> = H^3 |000>
# |1>_L = |---> = H^3 |111>
# Stabilizers: ZZI, IZZ (in the Hadamard-rotated basis)
plus = np.array([1, 1]) / np.sqrt(2)
minus = np.array([1, -1]) / np.sqrt(2)

def tensor_product(*states):
    result = states[0]
    for s in states[1:]:
        result = np.kron(result, s)
    return result

sv_0z_pf = tensor_product(plus, plus, plus)    # |+++>
sv_1z_pf = tensor_product(minus, minus, minus)  # |--->
sv_0x_pf = (sv_0z_pf + sv_1z_pf) / np.sqrt(2)
sv_1x_pf = (sv_0z_pf - sv_1z_pf) / np.sqrt(2)
sv_0y_pf = (sv_0z_pf + 1j * sv_1z_pf) / np.sqrt(2)
sv_1y_pf = (sv_0z_pf - 1j * sv_1z_pf) / np.sqrt(2)

# --- [[5,1,3]] code ---
n5 = 5
dim5 = 2**n5
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

state_0 = np.zeros(dim5); state_0[0] = 1.0
sv_0z_513 = proj_513 @ state_0
sv_0z_513 = sv_0z_513 / np.linalg.norm(sv_0z_513)
X_L_513 = multi_pauli({q: X for q in range(5)}, 5)
sv_1z_513 = X_L_513 @ sv_0z_513
sv_1z_513 = sv_1z_513 / np.linalg.norm(sv_1z_513)
sv_0x_513 = (sv_0z_513 + sv_1z_513) / np.sqrt(2)
sv_1x_513 = (sv_0z_513 - sv_1z_513) / np.sqrt(2)
sv_0y_513 = (sv_0z_513 + 1j * sv_1z_513) / np.sqrt(2)
sv_1y_513 = (sv_0z_513 - 1j * sv_1z_513) / np.sqrt(2)

# --- Uncoded single qubit ---
sv_0z_unc = np.array([1.0, 0.0])
sv_1z_unc = np.array([0.0, 1.0])
sv_0x_unc = np.array([1.0, 1.0]) / np.sqrt(2)
sv_1x_unc = np.array([1.0, -1.0]) / np.sqrt(2)
sv_0y_unc = np.array([1.0, 1j]) / np.sqrt(2)
sv_1y_unc = np.array([1.0, -1j]) / np.sqrt(2)

print(f"Code states built in {time.time()-start:.1f}s")

# ============================================================
# Experiment: Sweep lambda at fixed gamma, all bases
# ============================================================

gamma_fixed = 0.05  # mild T1 relaxation
lambdas = np.linspace(0.0, 0.50, 11)  # 0 to 0.5 dephasing

codes = {
    'uncoded': {
        'n': 1,
        'Z': (sv_0z_unc, sv_1z_unc),
        'X': (sv_0x_unc, sv_1x_unc),
        'Y': (sv_0y_unc, sv_1y_unc),
    },
    '3q-bitflip': {
        'n': n3,
        'Z': (sv_0z_bf, sv_1z_bf),
        'X': (sv_0x_bf, sv_1x_bf),
        'Y': (sv_0y_bf, sv_1y_bf),
    },
    '3q-phaseflip': {
        'n': n3,
        'Z': (sv_0z_pf, sv_1z_pf),
        'X': (sv_0x_pf, sv_1x_pf),
        'Y': (sv_0y_pf, sv_1y_pf),
    },
    '[[5,1,3]]': {
        'n': n5,
        'Z': (sv_0z_513, sv_1z_513),
        'X': (sv_0x_513, sv_1x_513),
        'Y': (sv_0y_513, sv_1y_513),
    },
}

print(f"\n=== Sweep λ at γ={gamma_fixed} — Holevo by basis ===\n")

all_results = {}
for lam in lambdas:
    key = f"lam={lam:.2f}"
    all_results[key] = {}
    print(f"\nλ={lam:.2f} (bias ratio λ/γ = {lam/gamma_fixed:.1f}):")

    for code_name, code_data in codes.items():
        n = code_data['n']
        basis_holevo = {}
        for basis in ['Z', 'X', 'Y']:
            sv0, sv1 = code_data[basis]
            rho_0 = density_matrix(sv0)
            rho_1 = density_matrix(sv1)
            rho_0n = apply_combined_noise(rho_0, gamma_fixed, lam, n)
            rho_1n = apply_combined_noise(rho_1, gamma_fixed, lam, n)
            h = holevo_info(rho_0n, rho_1n)
            basis_holevo[basis] = round(float(h), 4)

        avg = round(float(np.mean([basis_holevo[b] for b in ['Z', 'X', 'Y']])), 4)
        asym = round(float(max(basis_holevo.values()) - min(basis_holevo.values())), 4)

        all_results[key][code_name] = {
            **basis_holevo,
            'average': avg,
            'asymmetry': asym,
        }
        print(f"  {code_name:14s}: Z={basis_holevo['Z']:.3f}  X={basis_holevo['X']:.3f}  Y={basis_holevo['Y']:.3f}  avg={avg:.3f}  asym={asym:.3f}")

# ============================================================
# Also sweep gamma at fixed lambda (X-biased noise)
# ============================================================
print(f"\n\n=== Sweep γ at λ=0.05 (X-biased noise) — Holevo by basis ===\n")

lambda_fixed = 0.05
gammas = np.linspace(0.0, 0.50, 11)

x_bias_results = {}
for gamma in gammas:
    key = f"gam={gamma:.2f}"
    x_bias_results[key] = {}
    print(f"\nγ={gamma:.2f} (bias ratio γ/λ = {gamma/lambda_fixed:.1f}):")

    for code_name, code_data in codes.items():
        n = code_data['n']
        basis_holevo = {}
        for basis in ['Z', 'X', 'Y']:
            sv0, sv1 = code_data[basis]
            rho_0 = density_matrix(sv0)
            rho_1 = density_matrix(sv1)
            rho_0n = apply_combined_noise(rho_0, gamma, lambda_fixed, n)
            rho_1n = apply_combined_noise(rho_1, gamma, lambda_fixed, n)
            h = holevo_info(rho_0n, rho_1n)
            basis_holevo[basis] = round(float(h), 4)

        avg = round(float(np.mean([basis_holevo[b] for b in ['Z', 'X', 'Y']])), 4)
        asym = round(float(max(basis_holevo.values()) - min(basis_holevo.values())), 4)

        x_bias_results[key][code_name] = {
            **basis_holevo,
            'average': avg,
            'asymmetry': asym,
        }
        print(f"  {code_name:14s}: Z={basis_holevo['Z']:.3f}  X={basis_holevo['X']:.3f}  Y={basis_holevo['Y']:.3f}  avg={avg:.3f}  asym={asym:.3f}")

# ============================================================
# Summary: when does specialization beat isotropy?
# ============================================================
print("\n\n=== KEY ANALYSIS: Basis-averaged Holevo winners ===\n")

print("Z-BIASED NOISE (γ=0.05, varying λ):")
print(f"{'λ':>6s} {'λ/γ':>6s} | {'uncoded':>8s} {'3q-bf':>8s} {'3q-pf':>8s} {'[[5,1,3]]':>10s} | WINNER")
print("-" * 70)
for lam in lambdas:
    key = f"lam={lam:.2f}"
    avgs = {name: all_results[key][name]['average'] for name in codes}
    winner = max(avgs, key=avgs.get)
    ratio = lam / gamma_fixed
    print(f"  {lam:.2f} {ratio:6.1f} | {avgs['uncoded']:8.3f} {avgs['3q-bitflip']:8.3f} {avgs['3q-phaseflip']:8.3f} {avgs['[[5,1,3]]']:10.3f} | {winner}")

print("\nX-BIASED NOISE (λ=0.05, varying γ):")
print(f"{'γ':>6s} {'γ/λ':>6s} | {'uncoded':>8s} {'3q-bf':>8s} {'3q-pf':>8s} {'[[5,1,3]]':>10s} | WINNER")
print("-" * 70)
for gamma in gammas:
    key = f"gam={gamma:.2f}"
    avgs = {name: x_bias_results[key][name]['average'] for name in codes}
    winner = max(avgs, key=avgs.get)
    ratio = gamma / lambda_fixed
    print(f"  {gamma:.2f} {ratio:6.1f} | {avgs['uncoded']:8.3f} {avgs['3q-bitflip']:8.3f} {avgs['3q-phaseflip']:8.3f} {avgs['[[5,1,3]]']:10.3f} | {winner}")

elapsed = time.time() - start
print(f"\nTotal time: {elapsed:.1f}s")

# Save results
results = {
    'sprint': '022a',
    'description': 'Specialized (bias-aware) codes vs isotropic [[5,1,3]] under biased noise',
    'gamma_fixed': gamma_fixed,
    'lambda_fixed': lambda_fixed,
    'z_biased_results': all_results,
    'x_biased_results': x_bias_results,
    'elapsed_seconds': round(elapsed, 1)
}

with open('results/sprint_022a_biased_codes.json', 'w') as f:
    json.dump(results, f, indent=2)
print("Results saved to results/sprint_022a_biased_codes.json")
