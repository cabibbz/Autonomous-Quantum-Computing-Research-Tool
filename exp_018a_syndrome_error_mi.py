"""
Sprint 018a: Syndrome-Error Mutual Information.

For each QEC code, the syndrome is a deterministic function of the error.
Key quantities:
  H(S) = syndrome entropy (total info carried by syndrome)
  H(E|S) = residual error ambiguity (error uncertainty after observing syndrome)
  P(correct) = probability that syndrome-based correction succeeds
  MI(S : logical_outcome) = how much syndrome predicts success vs failure

Tests 3-qubit bit-flip code and [[5,1,3]] perfect code under their
respective noise models at varying error rates.
"""

import numpy as np
import json
from itertools import product

# Pauli matrices
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)

def kron_list(mats):
    result = mats[0]
    for m in mats[1:]:
        result = np.kron(result, m)
    return result

def pauli_string(paulis, n):
    pmap = {'I': I2, 'X': X, 'Y': Y, 'Z': Z}
    return kron_list([pmap[c] for c in paulis])

def h_binary(x):
    """Binary entropy."""
    if x <= 0 or x >= 1:
        return 0.0
    return -x * np.log2(x) - (1 - x) * np.log2(1 - x)

def entropy_from_dist(dist):
    """Shannon entropy of a probability distribution (bits)."""
    return -sum(p * np.log2(p) for p in dist if p > 1e-30)


# ============================================================
# 3-QUBIT BIT-FLIP CODE
# ============================================================
def three_qubit_analysis(p):
    """Analyze syndrome information for 3-qubit bit-flip code at error rate p."""
    n = 3
    # All 2^3 = 8 error patterns (each qubit flips or not)
    # Error pattern = (f0, f1, f2) where fi=1 means qubit i flips

    # Syndromes: s0 = Z0Z1 parity, s1 = Z1Z2 parity
    # No error:  s = (0,0) -> correct = no flip
    # Flip q0:   s = (1,0) -> correct = flip q0
    # Flip q1:   s = (1,1) -> correct = flip q1
    # Flip q2:   s = (0,1) -> correct = flip q2

    syndrome_table = {
        (0,0,0): (0,0),  # no error
        (1,0,0): (1,0),  # flip q0
        (0,1,0): (1,1),  # flip q1
        (0,0,1): (0,1),  # flip q2
        (1,1,0): (0,1),  # flips q0,q1 -> syndrome same as flip q2 -> MISCORRECT
        (1,0,1): (1,1),  # flips q0,q2 -> syndrome same as flip q1 -> MISCORRECT
        (0,1,1): (1,0),  # flips q1,q2 -> syndrome same as flip q0 -> MISCORRECT
        (1,1,1): (0,0),  # all flip -> syndrome 00 (no correction) -> MISCORRECT
    }

    # Correction table: syndrome -> which qubit to flip
    correction = {(0,0): None, (1,0): 0, (1,1): 1, (0,1): 2}

    # After correction, does the logical state survive?
    # Logical state = majority vote. If total flips (original + correction) is even -> correct
    error_probs = {}
    syndrome_probs = {(0,0): 0, (1,0): 0, (1,1): 0, (0,1): 0}
    p_correct = 0.0

    for f0, f1, f2 in product([0,1], repeat=3):
        n_flips = f0 + f1 + f2
        prob = (p ** n_flips) * ((1-p) ** (n - n_flips))
        error = (f0, f1, f2)
        syn = syndrome_table[error]

        error_probs[error] = prob
        syndrome_probs[syn] += prob

        # After correction
        corr_qubit = correction[syn]
        final_flips = list(error)
        if corr_qubit is not None:
            final_flips[corr_qubit] ^= 1

        # Logical correct if majority vote gives original (even total flips)
        if sum(final_flips) <= 1:  # majority unflipped
            p_correct += prob

    # Compute H(S)
    h_syndrome = entropy_from_dist(list(syndrome_probs.values()))

    # Compute H(E)
    h_error = entropy_from_dist(list(error_probs.values()))

    # H(E|S) = H(E) - H(S) since S = f(E)
    h_error_given_syndrome = h_error - h_syndrome

    # MI(S : logical_outcome)
    # logical_outcome is binary: correct or not
    p_wrong = 1.0 - p_correct
    h_logical = h_binary(p_correct)

    # P(outcome | syndrome) for each syndrome
    h_logical_given_s = 0.0
    for syn in syndrome_probs:
        p_s = syndrome_probs[syn]
        if p_s < 1e-30:
            continue
        # Which errors map to this syndrome?
        p_correct_given_s = 0.0
        for error, s in syndrome_table.items():
            if s == syn:
                p_e = error_probs[error]
                # Does correction work for this error?
                corr_qubit = correction[syn]
                final_flips = list(error)
                if corr_qubit is not None:
                    final_flips[corr_qubit] ^= 1
                if sum(final_flips) <= 1:
                    p_correct_given_s += p_e / p_s
        h_logical_given_s += p_s * h_binary(p_correct_given_s)

    mi_syndrome_logical = h_logical - h_logical_given_s

    return {
        'h_syndrome': round(h_syndrome, 6),
        'h_error': round(h_error, 6),
        'h_error_given_syndrome': round(max(0, h_error_given_syndrome), 6),
        'p_correct': round(p_correct, 6),
        'h_logical_outcome': round(h_logical, 6),
        'mi_syndrome_logical': round(mi_syndrome_logical, 6),
        'syndrome_probs': {str(k): round(v, 6) for k, v in syndrome_probs.items()},
        'n_syndromes': 4,
        'max_h_syndrome': 2.0,  # log2(4)
    }


# ============================================================
# [[5,1,3]] CODE
# ============================================================
def five_qubit_analysis(p):
    """Analyze syndrome information for [[5,1,3]] code under depolarizing noise."""
    n = 5
    stabilizers = ['XZZXI', 'IXZZX', 'XIXZZ', 'ZXIXZ']
    paulis = [I2, X, Y, Z]
    pauli_labels = ['I', 'X', 'Y', 'Z']

    # Build stabilizer matrices
    stab_mats = [pauli_string(s, n) for s in stabilizers]

    # For each single-qubit error (I,X,Y,Z on each qubit):
    # Compute syndrome by checking commutation with stabilizers
    # S_k anticommutes with error E => syndrome bit k = 1

    # All 4^5 = 1024 error patterns
    syndrome_probs = {}
    error_probs = {}
    correction_success = {}  # syndrome -> {error_pattern: prob, ...}
    p_correct_total = 0.0

    # For depolarizing: each qubit independently gets I (prob 1-p), X (prob p/3), Y (prob p/3), Z (prob p/3)
    for error_indices in product(range(4), repeat=n):
        # Probability
        n_errors = sum(1 for e in error_indices if e != 0)
        prob = ((1 - p) ** (n - n_errors)) * ((p / 3) ** n_errors)

        error_label = ''.join(pauli_labels[e] for e in error_indices)
        error_probs[error_label] = prob

        # Build error operator
        error_ops = [paulis[e] for e in error_indices]
        E = kron_list(error_ops)

        # Compute syndrome
        syn = 0
        for k, S in enumerate(stab_mats):
            # Check if E anticommutes with S: SE = -ES => SES^dag = -E
            comm = S @ E + E @ S  # anticommutator for anticommutation
            # Actually: if {S, E} = 0 (anticommute), syndrome bit = 1
            # Check: S @ E vs E @ S
            diff = S @ E - E @ S
            if np.linalg.norm(diff) > 1e-10:  # anticommute
                syn |= (1 << k)

        if syn not in syndrome_probs:
            syndrome_probs[syn] = 0.0
            correction_success[syn] = {'correct': 0.0, 'incorrect': 0.0}
        syndrome_probs[syn] += prob

        # Does standard correction work?
        # Standard correction: for syndrome s, apply the minimum-weight error with that syndrome
        # For [[5,1,3]], all weight-1 errors give distinct syndromes
        # Weight-0: syndrome 0 -> no correction -> correct
        # Weight-1: unique syndrome -> correct correction -> correct
        # Weight >= 2: may give same syndrome as different weight-1 error -> miscorrection

        if n_errors == 0:
            p_correct_total += prob
            correction_success[syn]['correct'] += prob
        elif n_errors == 1:
            # Single-qubit error: uniquely identified by syndrome -> correct
            p_correct_total += prob
            correction_success[syn]['correct'] += prob
        else:
            # Multi-qubit error: may or may not be correctable
            # The correction applies a weight-1 error matching the syndrome
            # Need to check if E * correction is in the stabilizer group (logical identity)
            # or is a logical operator (logical error)

            # Find the weight-1 error with this syndrome
            corr_error = None
            if syn == 0:
                # No correction applied. Multi-qubit error with syndrome 0
                # is either in stabilizer group (correct) or logical operator (error)
                corr_error = np.eye(2**n, dtype=complex)
            else:
                for q in range(n):
                    for pi, P in enumerate([X, Y, Z]):
                        test_ops = [I2] * n
                        test_ops[q] = P
                        test_E = kron_list(test_ops)
                        test_syn = 0
                        for k, S in enumerate(stab_mats):
                            diff = S @ test_E - test_E @ S
                            if np.linalg.norm(diff) > 1e-10:
                                test_syn |= (1 << k)
                        if test_syn == syn:
                            corr_error = test_E
                            break
                    if corr_error is not None:
                        break

            if corr_error is None:
                corr_error = np.eye(2**n, dtype=complex)

            # Net effect: corr_error @ E
            # This should be in stabilizer group (correct) or logical op (error)
            net = corr_error @ E

            # Check if net commutes with all stabilizers (-> in normalizer)
            # and check if it's in the stabilizer group or a logical operator
            in_normalizer = True
            for S in stab_mats:
                diff = S @ net - net @ S
                if np.linalg.norm(diff) > 1e-10:
                    in_normalizer = False
                    break

            if in_normalizer:
                # In normalizer: either stabilizer (correct) or logical (error)
                # Check against logical operators X_L = XXXXX, Z_L = ZZZZZ
                # If net is proportional to identity in codespace -> correct
                # Build codespace projector
                proj = np.eye(2**n, dtype=complex)
                for S in stab_mats:
                    proj = (np.eye(2**n) + S) / 2 @ proj

                # Check: proj @ net @ proj == c * proj for some c
                restricted = proj @ net @ proj
                # If net is in stabilizer group, restricted = proj (up to phase)
                # If net is logical op, restricted has off-diagonal blocks

                # Simple check: is net @ |0_L> = ±|0_L>?
                state_0 = np.zeros(2**n, dtype=complex)
                state_0[0] = 1.0
                for S in stab_mats:
                    state_0 = (np.eye(2**n) + S) / 2 @ state_0
                    norm = np.linalg.norm(state_0)
                    if norm > 1e-15:
                        state_0 = state_0 / norm

                result_state = net @ state_0
                overlap = abs(np.dot(state_0.conj(), result_state))

                if overlap > 0.99:
                    p_correct_total += prob
                    correction_success[syn]['correct'] += prob
                else:
                    correction_success[syn]['incorrect'] += prob
            else:
                # Not in normalizer -> this shouldn't happen for [[5,1,3]]
                # but treat as incorrect
                correction_success[syn]['incorrect'] += prob

    # Compute information measures
    h_syndrome = entropy_from_dist(list(syndrome_probs.values()))
    h_error = entropy_from_dist(list(error_probs.values()))
    h_error_given_syndrome = h_error - h_syndrome

    p_wrong = 1.0 - p_correct_total
    h_logical = h_binary(p_correct_total)

    # MI(S : logical_outcome)
    h_logical_given_s = 0.0
    for syn, probs in correction_success.items():
        p_s = syndrome_probs[syn]
        if p_s < 1e-30:
            continue
        p_corr_given_s = probs['correct'] / p_s
        h_logical_given_s += p_s * h_binary(p_corr_given_s)

    mi_syndrome_logical = h_logical - h_logical_given_s

    return {
        'h_syndrome': round(h_syndrome, 6),
        'h_error': round(h_error, 6),
        'h_error_given_syndrome': round(max(0, h_error_given_syndrome), 6),
        'p_correct': round(p_correct_total, 6),
        'h_logical_outcome': round(h_logical, 6),
        'mi_syndrome_logical': round(mi_syndrome_logical, 6),
        'n_syndromes': len(syndrome_probs),
        'max_h_syndrome': round(np.log2(16), 6),  # 4 syndrome bits -> 16 outcomes
    }


# ============================================================
# MAIN EXPERIMENT
# ============================================================
p_values = [0.0, 0.01, 0.02, 0.05, 0.08, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50]

print("=== 3-QUBIT BIT-FLIP CODE: SYNDROME INFORMATION ===\n")
print(f"{'p':>5} | {'H(S)':>6} | {'H(E)':>6} | {'H(E|S)':>6} | {'P(corr)':>7} | {'H(L)':>6} | {'MI(S:L)':>7}")
print("-" * 62)

results_3q = []
for p in p_values:
    r = three_qubit_analysis(p)
    results_3q.append(r)
    print(f"{p:5.2f} | {r['h_syndrome']:6.3f} | {r['h_error']:6.3f} | {r['h_error_given_syndrome']:6.3f} | "
          f"{r['p_correct']:7.4f} | {r['h_logical_outcome']:6.3f} | {r['mi_syndrome_logical']:7.4f}")

print(f"\n=== [[5,1,3]] CODE: SYNDROME INFORMATION (depolarizing) ===\n")
print(f"{'p':>5} | {'H(S)':>6} | {'H(E)':>6} | {'H(E|S)':>6} | {'P(corr)':>7} | {'H(L)':>6} | {'MI(S:L)':>7}")
print("-" * 62)

results_5q = []
for p in p_values:
    r = five_qubit_analysis(p)
    results_5q.append(r)
    print(f"{p:5.2f} | {r['h_syndrome']:6.3f} | {r['h_error']:6.3f} | {r['h_error_given_syndrome']:6.3f} | "
          f"{r['p_correct']:7.4f} | {r['h_logical_outcome']:6.3f} | {r['mi_syndrome_logical']:7.4f}")


# ============================================================
# Analysis: Syndrome efficiency
# ============================================================
print("\n=== SYNDROME EFFICIENCY ANALYSIS ===\n")
print("Syndrome efficiency = H(S) / H(E) = fraction of error info captured by syndrome")
print()
for i, p in enumerate(p_values):
    if p == 0:
        continue
    r3 = results_3q[i]
    r5 = results_5q[i]
    eff3 = r3['h_syndrome'] / r3['h_error'] if r3['h_error'] > 0 else 0
    eff5 = r5['h_syndrome'] / r5['h_error'] if r5['h_error'] > 0 else 0
    print(f"p={p:.2f}: 3-qubit eff={eff3:.3f} ({r3['h_syndrome']:.3f}/{r3['h_error']:.3f}), "
          f"[[5,1,3]] eff={eff5:.3f} ({r5['h_syndrome']:.3f}/{r5['h_error']:.3f})")

print("\nSyndrome predictive power = MI(S:logical_outcome) / H(logical_outcome)")
print("  = how much knowing the syndrome tells you about whether correction will succeed")
print()
for i, p in enumerate(p_values):
    if p == 0:
        continue
    r3 = results_3q[i]
    r5 = results_5q[i]
    pred3 = r3['mi_syndrome_logical'] / r3['h_logical_outcome'] if r3['h_logical_outcome'] > 0.001 else 0
    pred5 = r5['mi_syndrome_logical'] / r5['h_logical_outcome'] if r5['h_logical_outcome'] > 0.001 else 0
    print(f"p={p:.2f}: 3-qubit pred={pred3:.3f}, [[5,1,3]] pred={pred5:.3f}")


# Save results
all_results = {
    'p_values': p_values,
    'three_qubit_bitflip': results_3q,
    'five_qubit_depolarizing': results_5q,
    'description': 'Syndrome-error mutual information analysis for QEC codes',
    'measures': {
        'h_syndrome': 'Shannon entropy of syndrome distribution H(S)',
        'h_error': 'Shannon entropy of error distribution H(E)',
        'h_error_given_syndrome': 'Residual error entropy H(E|S) = H(E) - H(S)',
        'p_correct': 'Probability of successful correction',
        'h_logical_outcome': 'Binary entropy of correction success H(L)',
        'mi_syndrome_logical': 'MI between syndrome and correction outcome MI(S:L)',
    }
}

with open('results/sprint_018a_syndrome_error_mi.json', 'w') as f:
    json.dump(all_results, f, indent=2)

print("\nResults saved to results/sprint_018a_syndrome_error_mi.json")
