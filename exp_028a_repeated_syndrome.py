"""
Sprint 028a: Repeated syndrome extraction with majority vote for [[5,1,3]].

Sprints 026-027 proved that single-round syndrome extraction (bare or flag-FT)
cannot beat passive encoding. The root cause: one noisy syndrome measurement
provides insufficient information for reliable correction.

This experiment builds a REPEATED syndrome extraction circuit:
- Extract syndrome d=3 times using fresh ancillas each round
- Take majority vote across rounds for each stabilizer
- Use majority syndrome for correction

Key design choices:
- Fresh ancillas each round (reset would require mid-circuit reset, not universally available)
- 3 rounds for distance-3 code (standard FT protocol)
- Total: 5 data + 4×3=12 ancilla = 17 qubits (within simulator limits)
- Total 2Q gates: 16 per round × 3 = 48

This experiment:
1. Builds the repeated syndrome circuit
2. Verifies majority vote gives correct syndrome in noiseless case
3. Tests with single-qubit errors to confirm all corrections work
4. Tests with single ancilla errors to confirm majority vote corrects them
"""

import numpy as np
import json, time
from qiskit.circuit import QuantumCircuit
from qiskit.quantum_info import Clifford, Statevector
from qiskit_aer import AerSimulator

start = time.time()

# ============================================================
# [[5,1,3]] Encoding circuit (from Sprint 025a)
# ============================================================
n = 5
destab_x = np.array([[1,1,1,1,1],[1,0,1,1,1],[0,1,1,0,0],[1,0,0,0,1],[0,1,1,1,1]], dtype=bool)
destab_z = np.zeros((5,5), dtype=bool)
stab_x = np.array([[0,0,0,0,0],[1,0,0,1,0],[0,1,0,0,1],[1,0,1,0,0],[0,1,0,1,0]], dtype=bool)
stab_z = np.array([[1,1,1,1,1],[0,1,1,0,0],[0,0,1,1,0],[0,0,0,1,1],[1,0,0,0,1]], dtype=bool)
tableau = np.zeros((10, 11), dtype=bool)
tableau[:5, :5] = destab_x; tableau[:5, 5:10] = destab_z
tableau[5:, :5] = stab_x;  tableau[5:, 5:10] = stab_z
enc_513 = Clifford(tableau).to_circuit()

# ============================================================
# Stabilizer definitions
# ============================================================
stabilizers = [
    [(0,'X'), (1,'Z'), (2,'Z'), (3,'X')],   # g1 = XZZXI
    [(1,'X'), (2,'Z'), (3,'Z'), (4,'X')],   # g2 = IXZZX
    [(0,'X'), (2,'X'), (3,'Z'), (4,'Z')],   # g3 = XIXZZ
    [(0,'Z'), (1,'X'), (3,'X'), (4,'Z')],   # g4 = ZXIXZ
]

I2 = np.eye(2, dtype=complex)
Xm = np.array([[0,1],[1,0]], dtype=complex)
Ym = np.array([[0,-1j],[1j,0]], dtype=complex)
Zm = np.array([[1,0],[0,-1]], dtype=complex)
pauli_map = {'I': I2, 'X': Xm, 'Y': Ym, 'Z': Zm}

def multi_kron(ops, nq):
    r = np.eye(1, dtype=complex)
    for q in range(nq):
        r = np.kron(r, ops.get(q, I2))
    return r

stab_strings = ['XZZXI', 'IXZZX', 'XIXZZ', 'ZXIXZ']
stab_matrices = []
for ss in stab_strings:
    ops = {q: pauli_map[ss[q]] for q in range(5) if ss[q] != 'I'}
    stab_matrices.append(multi_kron(ops, 5))

def compute_syndrome(error_op):
    syndrome = []
    for S in stab_matrices:
        commutator = error_op @ S - S @ error_op
        anticommutes = np.linalg.norm(commutator) > 1e-10
        syndrome.append(1 if anticommutes else 0)
    return tuple(syndrome)

# Build syndrome -> correction lookup
syndrome_table = {(0,0,0,0): None}
for q in range(5):
    for p_label, p_mat in [('X', Xm), ('Y', Ym), ('Z', Zm)]:
        E = multi_kron({q: p_mat}, 5)
        syn = compute_syndrome(E)
        syndrome_table[syn] = (q, p_label)

# ============================================================
# Build repeated syndrome extraction circuit
# ============================================================
def build_single_round_syndrome(n_data, anc_start):
    """Build one round of syndrome extraction.
    Data qubits: 0..n_data-1. Ancillas: anc_start..anc_start+3."""
    qc = QuantumCircuit(anc_start + 4)
    for s_idx, stab in enumerate(stabilizers):
        anc = anc_start + s_idx
        qc.h(anc)
        for data_q, pauli in stab:
            if pauli == 'X': qc.cx(anc, data_q)
            elif pauli == 'Z': qc.cz(anc, data_q)
        qc.h(anc)
    return qc

def build_repeated_syndrome(n_rounds=3):
    """Build repeated syndrome extraction with fresh ancillas per round.
    Total qubits: 5 data + 4*n_rounds ancilla.
    Returns (circuit, list of ancilla qubit indices per round)."""
    n_total = 5 + 4 * n_rounds
    qc = QuantumCircuit(n_total)

    round_ancillas = []
    for r in range(n_rounds):
        anc_start = 5 + 4 * r
        round_ancillas.append(list(range(anc_start, anc_start + 4)))

        # Add syndrome extraction for this round
        for s_idx, stab in enumerate(stabilizers):
            anc = anc_start + s_idx
            qc.h(anc)
            for data_q, pauli in stab:
                if pauli == 'X': qc.cx(anc, data_q)
                elif pauli == 'Z': qc.cz(anc, data_q)
            qc.h(anc)

        qc.barrier()  # Visual separation between rounds

    return qc, round_ancillas

# Build 3-round circuit
n_rounds = 3
syn_circuit, round_ancillas = build_repeated_syndrome(n_rounds)
n_total_qubits = 5 + 4 * n_rounds

print(f"Repeated syndrome circuit ({n_rounds} rounds):")
print(f"  Total qubits: {n_total_qubits}")
print(f"  Gates: {dict(syn_circuit.count_ops())}")
print(f"  Depth: {syn_circuit.depth()}")
cx_count = syn_circuit.count_ops().get('cx', 0)
cz_count = syn_circuit.count_ops().get('cz', 0)
print(f"  Total 2Q gates: {cx_count + cz_count}")
print(f"  Ancilla layout: {round_ancillas}")

# ============================================================
# Noiseless verification: single-qubit errors
# ============================================================
print("\n" + "="*60)
print("VERIFICATION 1: Single-qubit errors, noiseless")
print("="*60)

sim = AerSimulator(method='statevector')
shots = 10000

verification_results = []
n_correct = 0
n_total = 0

for q_err in range(5):
    for p_label in ['X', 'Y', 'Z']:
        for logical in ['0', '1']:
            # Build full circuit
            n_qb = n_total_qubits
            n_cl = 4 * n_rounds  # measure all ancillas
            qc = QuantumCircuit(n_qb, n_cl)

            # Prepare logical state
            if logical == '1': qc.x(0)

            # Encode
            qc.compose(enc_513, qubits=range(5), inplace=True)

            # Inject error
            if p_label == 'X': qc.x(q_err)
            elif p_label == 'Y': qc.y(q_err)
            elif p_label == 'Z': qc.z(q_err)

            # Repeated syndrome extraction
            qc.compose(syn_circuit, qubits=range(n_qb), inplace=True)

            # Measure all ancillas
            for r in range(n_rounds):
                for s in range(4):
                    cl_bit = r * 4 + s
                    qc.measure(round_ancillas[r][s], cl_bit)

            # Run
            counts = sim.run(qc, shots=shots).result().get_counts()

            # Should be deterministic in noiseless case
            dominant = max(counts, key=counts.get)
            dominant_frac = counts[dominant] / shots

            # Extract per-round syndromes and majority vote
            bits = dominant[::-1]  # LSB first
            round_syndromes = []
            for r in range(n_rounds):
                syn_r = tuple(int(bits[r*4 + s]) for s in range(4))
                round_syndromes.append(syn_r)

            # Majority vote per stabilizer
            majority_syn = []
            for s in range(4):
                votes = [round_syndromes[r][s] for r in range(n_rounds)]
                majority_syn.append(1 if sum(votes) > n_rounds / 2 else 0)
            majority_syn = tuple(majority_syn)

            # Expected syndrome
            E = multi_kron({q_err: pauli_map[p_label]}, 5)
            expected = compute_syndrome(E)

            match = majority_syn == expected
            all_rounds_match = all(rs == expected for rs in round_syndromes)

            n_total += 1
            if match: n_correct += 1

            verification_results.append({
                'error': f"{p_label} on q{q_err}",
                'logical': logical,
                'round_syndromes': [list(rs) for rs in round_syndromes],
                'majority_syndrome': list(majority_syn),
                'expected_syndrome': list(expected),
                'majority_match': match,
                'all_rounds_match': all_rounds_match,
                'dominant_frac': round(float(dominant_frac), 4),
            })

print(f"Majority vote correct: {n_correct}/{n_total}")
all_rounds_correct = sum(1 for r in verification_results if r['all_rounds_match'])
print(f"All rounds individually correct: {all_rounds_correct}/{n_total}")

# ============================================================
# Verification 2: Single ancilla errors (test majority vote)
# ============================================================
print("\n" + "="*60)
print("VERIFICATION 2: Single ancilla errors — majority vote should override")
print("="*60)

# For a codeword with no data error, inject a single X error on one ancilla
# Majority vote should still give (0,0,0,0)
ancilla_error_results = []
n_ancilla_correct = 0
n_ancilla_total = 0

for target_round in range(n_rounds):
    for target_stab in range(4):
        target_anc = round_ancillas[target_round][target_stab]

        for logical in ['0']:  # Just test |0>_L
            qc = QuantumCircuit(n_total_qubits, 4 * n_rounds)

            # Prepare and encode
            if logical == '1': qc.x(0)
            qc.compose(enc_513, qubits=range(5), inplace=True)

            # We need to inject the ancilla error at the right time
            # The ancilla error should happen DURING syndrome extraction
            # For simplicity, we'll inject it AFTER the H gate but before measurement
            # This is equivalent to a bit-flip on the ancilla readout

            # Build syndrome round by round, injecting error in the target round
            for r in range(n_rounds):
                anc_start = 5 + 4 * r
                for s_idx, stab in enumerate(stabilizers):
                    anc = anc_start + s_idx
                    qc.h(anc)
                    for data_q, pauli in stab:
                        if pauli == 'X': qc.cx(anc, data_q)
                        elif pauli == 'Z': qc.cz(anc, data_q)
                    qc.h(anc)

                # Inject ancilla error after this round's extraction
                if r == target_round:
                    qc.x(target_anc)  # Flip the ancilla readout

                qc.barrier()

            # Measure all ancillas
            for r in range(n_rounds):
                for s in range(4):
                    qc.measure(round_ancillas[r][s], r * 4 + s)

            # Run
            counts = sim.run(qc, shots=1000).result().get_counts()
            dominant = max(counts, key=counts.get)
            bits = dominant[::-1]

            # Extract syndromes per round
            round_syndromes = []
            for r in range(n_rounds):
                syn_r = tuple(int(bits[r*4 + s]) for s in range(4))
                round_syndromes.append(syn_r)

            # Majority vote
            majority_syn = []
            for s in range(4):
                votes = [round_syndromes[r][s] for r in range(n_rounds)]
                majority_syn.append(1 if sum(votes) > n_rounds / 2 else 0)
            majority_syn = tuple(majority_syn)

            # Expected: no data error, so syndrome should be (0,0,0,0)
            expected = (0, 0, 0, 0)
            match = majority_syn == expected

            # Check that the targeted round has the flipped bit
            target_round_syn = round_syndromes[target_round]
            flipped = target_round_syn[target_stab] == 1
            other_rounds_clean = all(
                round_syndromes[r] == (0,0,0,0)
                for r in range(n_rounds) if r != target_round
            )

            n_ancilla_total += 1
            if match: n_ancilla_correct += 1

            ancilla_error_results.append({
                'target_round': target_round,
                'target_stabilizer': target_stab,
                'target_ancilla': target_anc,
                'round_syndromes': [list(rs) for rs in round_syndromes],
                'majority_syndrome': list(majority_syn),
                'majority_correct': match,
                'target_bit_flipped': flipped,
                'other_rounds_clean': other_rounds_clean,
            })

print(f"Majority vote overrides ancilla error: {n_ancilla_correct}/{n_ancilla_total}")
n_flipped = sum(1 for r in ancilla_error_results if r['target_bit_flipped'])
print(f"Target bit correctly flipped in error round: {n_flipped}/{n_ancilla_total}")
n_clean = sum(1 for r in ancilla_error_results if r['other_rounds_clean'])
print(f"Other rounds clean: {n_clean}/{n_ancilla_total}")

# Print details for any failures
failures = [r for r in ancilla_error_results if not r['majority_correct']]
if failures:
    print(f"\nFAILURES ({len(failures)}):")
    for f in failures:
        print(f"  Round {f['target_round']}, stab {f['target_stabilizer']}: "
              f"syndromes={f['round_syndromes']}, majority={f['majority_syndrome']}")
else:
    print("ALL ancilla errors correctly overridden by majority vote!")

# ============================================================
# Verification 3: Ancilla error DURING syndrome extraction (mid-CNOT)
# ============================================================
print("\n" + "="*60)
print("VERIFICATION 3: Ancilla X error propagation through CX gates")
print("="*60)

# The real concern: an ancilla X error BEFORE the CX gates propagates to data
# With CX(anc, data), an X error on anc propagates X to data
# This creates a weight-2+ data error that may be uncorrectable
# But with majority vote, the other 2 rounds should give the correct syndrome

# We test: inject X on ancilla BEFORE the CX/CZ gates in one round
# The other rounds should still extract the correct (data) syndrome
# Even if one round's syndrome is corrupted, majority should work

propagation_results = []
n_prop_correct = 0

for target_round in range(n_rounds):
    for target_stab in range(4):
        target_anc = round_ancillas[target_round][target_stab]

        qc = QuantumCircuit(n_total_qubits, 4 * n_rounds)

        # Encode |0>_L
        qc.compose(enc_513, qubits=range(5), inplace=True)

        # Build syndrome round by round
        for r in range(n_rounds):
            anc_start = 5 + 4 * r
            for s_idx, stab in enumerate(stabilizers):
                anc = anc_start + s_idx
                qc.h(anc)

                # Inject error BEFORE CX/CZ gates (after H)
                if r == target_round and s_idx == target_stab:
                    qc.x(anc)  # This will propagate through subsequent CX gates

                for data_q, pauli in stab:
                    if pauli == 'X': qc.cx(anc, data_q)
                    elif pauli == 'Z': qc.cz(anc, data_q)
                qc.h(anc)
            qc.barrier()

        # Measure all ancillas
        for r in range(n_rounds):
            for s in range(4):
                qc.measure(round_ancillas[r][s], r * 4 + s)

        # Run
        counts = sim.run(qc, shots=1000).result().get_counts()
        dominant = max(counts, key=counts.get)
        dominant_frac = counts[dominant] / 1000
        bits = dominant[::-1]

        round_syndromes = []
        for r in range(n_rounds):
            syn_r = tuple(int(bits[r*4 + s]) for s in range(4))
            round_syndromes.append(syn_r)

        # Majority vote
        majority_syn = []
        for s in range(4):
            votes = [round_syndromes[r][s] for r in range(n_rounds)]
            majority_syn.append(1 if sum(votes) > n_rounds / 2 else 0)
        majority_syn = tuple(majority_syn)

        # The ancilla X error propagates through CX gates to DATA qubits
        # This means the data state is now corrupted with a multi-qubit error
        # The corrupted round will see the wrong syndrome (ancilla measurement also affected)
        # But subsequent rounds will see the syndrome of the DATA error caused by propagation
        # Majority vote may or may NOT give the right correction here

        # Key insight: the ancilla error propagates X to some data qubits,
        # changing the data state. So later rounds will detect a DATA error syndrome,
        # which may differ from the "correct" (0,0,0,0) for the original state.

        propagation_results.append({
            'target_round': target_round,
            'target_stabilizer': target_stab,
            'round_syndromes': [list(rs) for rs in round_syndromes],
            'majority_syndrome': list(majority_syn),
            'majority_is_trivial': majority_syn == (0,0,0,0),
            'dominant_frac': round(float(dominant_frac), 4),
        })

        if majority_syn != (0,0,0,0):
            # Check if the majority syndrome corresponds to a correctable error
            correction = syndrome_table.get(majority_syn)
            propagation_results[-1]['correction'] = str(correction)

        if majority_syn == (0,0,0,0):
            n_prop_correct += 1

print(f"Majority vote gives trivial syndrome (ideal): {n_prop_correct}/{len(propagation_results)}")

# Show what happens
for r in propagation_results:
    status = "OK" if r['majority_is_trivial'] else f"WRONG: majority={r['majority_syndrome']}"
    print(f"  Round {r['target_round']}, stab {r['target_stabilizer']}: "
          f"syndromes={r['round_syndromes']} -> {status}")

# ============================================================
# Save results
# ============================================================
elapsed = time.time() - start
print(f"\nTotal time: {elapsed:.1f}s")

output = {
    'sprint': '028a',
    'description': 'Repeated syndrome extraction circuit — build and verify',
    'circuit': {
        'n_rounds': n_rounds,
        'total_qubits': n_total_qubits,
        'gates': dict(syn_circuit.count_ops()),
        'depth': syn_circuit.depth(),
        'total_2q_gates': cx_count + cz_count,
        'ancilla_layout': round_ancillas,
    },
    'verification_data_errors': {
        'majority_correct': n_correct,
        'total': n_total,
        'all_rounds_correct': all_rounds_correct,
    },
    'verification_ancilla_errors': {
        'majority_overrides': n_ancilla_correct,
        'total': n_ancilla_total,
        'details': ancilla_error_results,
    },
    'verification_propagation': {
        'majority_trivial': n_prop_correct,
        'total': len(propagation_results),
        'details': propagation_results,
    },
    'elapsed_seconds': round(elapsed, 1),
}

with open('results/sprint_028a_repeated_syndrome.json', 'w') as f:
    json.dump(output, f, indent=2)
print("Saved to results/sprint_028a_repeated_syndrome.json")
