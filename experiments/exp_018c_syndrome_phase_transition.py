"""
Sprint 018c: Syndrome MI at the Phase Transition.

From 18a: syndrome is a lossy compression of the error.
From 18b: logical information (Holevo) undergoes a phase transition.
Question: does syndrome information undergo its own transition?

Key idea: the syndrome contains two types of information:
  1. "Useful" info — tells you which correctable error occurred
  2. "Useless" info — tells you something happened, but not what (ambiguous)

As noise increases, the syndrome transitions from mostly-useful to mostly-useless.
Does this transition align with the logical information transition?

Test: 3-qubit bit-flip code (clean threshold at p=0.5) and compare syndrome
information decomposition to Holevo information.
Also: concatenated codes to see if the syndrome transition sharpens.
"""

import numpy as np
import json
from itertools import product

def h_binary(x):
    if x <= 0 or x >= 1:
        return 0.0
    return -x * np.log2(x) - (1 - x) * np.log2(1 - x)

def entropy_from_dist(dist):
    return -sum(p * np.log2(p) for p in dist if p > 1e-30)


# ============================================================
# 3-QUBIT BIT-FLIP CODE: Detailed syndrome decomposition
# ============================================================
def three_qubit_detailed(p):
    """Full syndrome decomposition for 3-qubit code."""
    syndrome_table = {
        (0,0,0): (0,0), (1,0,0): (1,0), (0,1,0): (1,1), (0,0,1): (0,1),
        (1,1,0): (0,1), (1,0,1): (1,1), (0,1,1): (1,0), (1,1,1): (0,0),
    }
    correction = {(0,0): None, (1,0): 0, (1,1): 1, (0,1): 2}

    # Build joint distribution P(syndrome, correctable)
    joint = {}  # (syndrome, correct?) -> probability
    syndrome_probs = {}
    total_useful = 0.0  # P(syndrome uniquely identifies correctable error)

    for f0, f1, f2 in product([0,1], repeat=3):
        n_flips = f0 + f1 + f2
        prob = (p ** n_flips) * ((1-p) ** (3 - n_flips))
        error = (f0, f1, f2)
        syn = syndrome_table[error]

        # After correction, is the result correct?
        corr_qubit = correction[syn]
        final_flips = list(error)
        if corr_qubit is not None:
            final_flips[corr_qubit] ^= 1
        correct = (sum(final_flips) <= 1)

        key = (syn, correct)
        joint[key] = joint.get(key, 0) + prob
        syndrome_probs[syn] = syndrome_probs.get(syn, 0) + prob

    # Decompose syndrome info
    h_syndrome = entropy_from_dist(list(syndrome_probs.values()))

    # "Useful" syndrome info = MI(S : correctable_error)
    # For the 3-qubit code, each syndrome maps to exactly one correctable error
    # and one uncorrectable error pattern. So knowing S tells you:
    #   - which error correction is applied (always known)
    #   - whether that correction is right or wrong (NOT known from S alone)

    # Fraction of syndrome that's "disambiguating" vs "ambiguous"
    # For each syndrome, compute P(correct | syndrome)
    p_correct_total = 0.0
    p_correct_given_s = {}
    for syn in syndrome_probs:
        p_corr = joint.get((syn, True), 0) / syndrome_probs[syn]
        p_correct_given_s[syn] = p_corr
        p_correct_total += joint.get((syn, True), 0)

    # H(correct | S)
    h_correct_given_s = sum(
        syndrome_probs[syn] * h_binary(p_correct_given_s[syn])
        for syn in syndrome_probs
    )

    # H(correct)
    h_correct = h_binary(p_correct_total)

    # MI(S : correct) = H(correct) - H(correct | S)
    mi_s_correct = h_correct - h_correct_given_s

    # Holevo information (analytical for bit-flip on repetition code)
    # |0_L> = |000>, |1_L> = |111>
    # After bit-flip noise, both are diagonal in computational basis
    dim = 8
    rho0 = np.zeros(dim)
    rho1 = np.zeros(dim)
    for x in range(dim):
        bits = [(x >> i) & 1 for i in range(3)]
        n_ones = sum(bits)
        rho0[x] = (p ** n_ones) * ((1-p) ** (3 - n_ones))
        rho1[x] = (p ** (3 - n_ones)) * ((1-p) ** n_ones)

    rho_avg = (rho0 + rho1) / 2
    h_avg = entropy_from_dist(rho_avg)
    h0 = entropy_from_dist(rho0)
    h1 = entropy_from_dist(rho1)
    holevo = h_avg - 0.5 * h0 - 0.5 * h1

    return {
        'h_syndrome': round(h_syndrome, 6),
        'p_correct': round(p_correct_total, 6),
        'h_correct': round(h_correct, 6),
        'h_correct_given_s': round(h_correct_given_s, 6),
        'mi_s_correct': round(mi_s_correct, 6),
        'holevo': round(holevo, 6),
        'p_correct_per_syndrome': {str(k): round(v, 6) for k, v in p_correct_given_s.items()},
    }


# ============================================================
# 9-QUBIT CONCATENATED CODE: syndrome at two levels
# ============================================================
def nine_qubit_detailed(p):
    """Syndrome analysis for 2-level concatenated 3-qubit bit-flip code."""
    # Each qubit flips independently with prob p
    # Inner code: 3 blocks of 3 qubits, each with 2 syndrome bits
    # Outer code: 3 "super-qubits" (inner-decoded), with 2 syndrome bits

    # Inner decode: majority vote on each block
    # P(inner block correct) = (1-p)^3 + 3p(1-p)^2 = 1 - 3p^2 + 2p^3
    p_inner_correct = (1-p)**3 + 3*p*(1-p)**2

    # Inner syndrome: same as 3-qubit code at rate p
    inner_h_syn = 0
    inner_mi = 0

    # For each inner block, compute syndrome distribution
    # This is just the 3-qubit analysis at rate p
    inner = three_qubit_detailed(p)
    inner_h_syn_per_block = inner['h_syndrome']

    # Total inner syndrome entropy (3 independent blocks)
    total_inner_h_syn = 3 * inner_h_syn_per_block

    # Outer code: 3 "super-qubits", each with error rate p_inner_error = 1 - p_inner_correct
    p_outer = 1 - p_inner_correct  # = 3p^2 - 2p^3
    outer = three_qubit_detailed(p_outer)
    outer_h_syn = outer['h_syndrome']

    # Total syndrome entropy
    total_h_syn = total_inner_h_syn + outer_h_syn

    # Total correction probability
    # Inner corrects each block independently, then outer corrects the super-qubits
    # P(overall correct) = p_inner_correct applied at the outer level
    p_outer_correct = (1 - p_outer)**3 + 3*p_outer*(1-p_outer)**2
    # Which simplifies to the concatenated formula

    # Holevo for 9-qubit code
    dim = 2**9
    rho0 = np.zeros(dim)
    rho1 = np.zeros(dim)
    for x in range(dim):
        bits = [(x >> i) & 1 for i in range(9)]
        n_ones = sum(bits)
        rho0[x] = (p ** n_ones) * ((1-p) ** (9 - n_ones))
        rho1[x] = (p ** (9 - n_ones)) * ((1-p) ** n_ones)

    rho_avg = (rho0 + rho1) / 2
    holevo = entropy_from_dist(rho_avg) - 0.5*entropy_from_dist(rho0) - 0.5*entropy_from_dist(rho1)

    return {
        'total_h_syndrome': round(total_h_syn, 6),
        'inner_h_syndrome': round(total_inner_h_syn, 6),
        'outer_h_syndrome': round(outer_h_syn, 6),
        'p_inner_correct': round(p_inner_correct, 6),
        'p_outer_error': round(p_outer, 6),
        'p_overall_correct': round(p_outer_correct, 6),
        'holevo': round(holevo, 6),
    }


# ============================================================
# MAIN: Fine-grained sweep around threshold
# ============================================================
p_values = np.concatenate([
    np.arange(0, 0.10, 0.01),
    np.arange(0.10, 0.50, 0.02),
    np.arange(0.50, 0.55, 0.01),
    [0.60, 0.70, 0.80, 0.90, 1.0]
])
p_values = sorted(set(np.round(p_values, 4)))

print("=== SYNDROME INFORMATION AT THE PHASE TRANSITION ===\n")
print(f"{'p':>5} | {'H(S)3q':>7} | {'MI(S:L)':>7} | {'P(corr)':>7} | {'Holevo':>7} | {'H(S)9q':>7} | {'P(c)9q':>7} | {'Hol_9q':>7}")
print("-" * 80)

results_3q = []
results_9q = []

for p in p_values:
    r3 = three_qubit_detailed(p)
    r9 = nine_qubit_detailed(p)
    results_3q.append(r3)
    results_9q.append(r9)

    print(f"{p:5.2f} | {r3['h_syndrome']:7.4f} | {r3['mi_s_correct']:7.4f} | {r3['p_correct']:7.4f} | "
          f"{r3['holevo']:7.4f} | {r9['total_h_syndrome']:7.4f} | {r9['p_overall_correct']:7.4f} | {r9['holevo']:7.4f}")


# ============================================================
# Analysis: Syndrome-Holevo correlation
# ============================================================
print("\n\n=== SYNDROME vs HOLEVO: PHASE TRANSITION COMPARISON ===\n")
print("Do syndrome measures show the same transition as Holevo?")
print()
print(f"{'p':>5} | {'H(S)/Hmax':>9} | {'Holevo':>7} | {'MI(S:L)/H(L)':>12} | {'P(correct)':>10}")
print("-" * 55)
for i, p in enumerate(p_values):
    r = results_3q[i]
    h_norm = r['h_syndrome'] / 2.0  # max syndrome entropy = 2 bits
    pred = r['mi_s_correct'] / r['h_correct'] if r['h_correct'] > 0.001 else 0
    print(f"{p:5.2f} | {h_norm:9.4f} | {r['holevo']:7.4f} | {pred:12.4f} | {r['p_correct']:10.4f}")

# Key crossover points
print("\n=== KEY TRANSITION POINTS ===")
for i in range(1, len(p_values)):
    p = p_values[i]
    r = results_3q[i]
    r_prev = results_3q[i-1]
    # Holevo = 0.5 (half information retained)
    if r['holevo'] <= 0.5 and r_prev['holevo'] > 0.5:
        print(f"Holevo = 0.5 (half info retained) at p ≈ {p:.3f}")
    # P(correct) = 0.75
    if r['p_correct'] <= 0.75 and r_prev['p_correct'] > 0.75:
        print(f"P(correct) = 0.75 at p ≈ {p:.3f}")
    # H(S) = 1.5 (75% of max)
    if r['h_syndrome'] >= 1.5 and r_prev['h_syndrome'] < 1.5:
        print(f"H(S) = 1.5 (75% of max) at p ≈ {p:.3f}")

# Concatenation sharpening
print("\n=== CONCATENATION SHARPENING ===")
print(f"{'p':>5} | {'3q Holevo':>9} | {'9q Holevo':>9} | {'Ratio':>6} | {'3q P(c)':>7} | {'9q P(c)':>7}")
print("-" * 55)
for i, p in enumerate(p_values):
    r3 = results_3q[i]
    r9 = results_9q[i]
    ratio = r9['holevo'] / r3['holevo'] if r3['holevo'] > 0.001 else 0
    print(f"{p:5.2f} | {r3['holevo']:9.4f} | {r9['holevo']:9.4f} | {ratio:6.2f} | {r3['p_correct']:7.4f} | {r9['p_overall_correct']:7.4f}")

# Syndrome hierarchy at concatenation
print("\n=== SYNDROME HIERARCHY IN 9-QUBIT CODE ===")
print(f"{'p':>5} | {'Inner H(S)':>10} | {'Outer H(S)':>10} | {'Total H(S)':>10} | {'Outer/Total':>11}")
print("-" * 55)
for i, p in enumerate(p_values):
    r9 = results_9q[i]
    total = r9['total_h_syndrome']
    ratio = r9['outer_h_syndrome'] / total if total > 0.001 else 0
    print(f"{p:5.2f} | {r9['inner_h_syndrome']:10.4f} | {r9['outer_h_syndrome']:10.4f} | {total:10.4f} | {ratio:11.4f}")


# Save
all_results = {
    'p_values': [round(float(p), 4) for p in p_values],
    'three_qubit': results_3q,
    'nine_qubit': results_9q,
    'description': 'Syndrome information at the phase transition for 3-qubit and 9-qubit repetition codes',
}

with open('results/sprint_018c_syndrome_phase_transition.json', 'w') as f:
    json.dump(all_results, f, indent=2)

print("\nResults saved to results/sprint_018c_syndrome_phase_transition.json")
