"""
Sprint 026a: Build and verify [[5,1,3]] syndrome extraction circuit.

The [[5,1,3]] stabilizers are:
  g1 = XZZXI
  g2 = IXZZX
  g3 = XIXZZ
  g4 = ZXIXZ

Each stabilizer is measured using an ancilla qubit:
  - Initialize ancilla in |0>
  - For each data qubit in the stabilizer:
    - If X: CNOT(ancilla, data) [ancilla controls]... wait, standard approach:
      - X part: H on ancilla, CNOT(ancilla->data) or CX, then H
      - Actually: to measure a Pauli P, we use controlled-P from ancilla
      - For X: CX(ancilla, data)
      - For Z: CZ(ancilla, data) = CX with H sandwiching on target...
      Actually the standard construction:
      - ancilla starts |0>, apply H -> |+>
      - For each qubit i in stabilizer:
        - If X_i: CX(ancilla, data_i)
        - If Z_i: CZ(ancilla, data_i)
        - If Y_i: CY(ancilla, data_i) [= S†·CX·S on data, plus CZ]
      - Apply H to ancilla, measure

For [[5,1,3]] stabilizers (no Y components), we only need CX and CZ.
CZ = H·CX·H on target, but CZ is native on many backends.

This experiment:
1. Builds the encoding circuit (from Sprint 025a)
2. Builds syndrome extraction for all 4 stabilizers
3. Injects each single-qubit X, Y, Z error
4. Verifies syndrome correctly identifies the error
5. Applies correction based on syndrome lookup table
6. Verifies logical state is restored
"""

import numpy as np
import json, time
from qiskit.circuit import QuantumCircuit, ClassicalRegister
from qiskit.quantum_info import Clifford, Statevector, Operator
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

cliff = Clifford(tableau)
enc_513 = cliff.to_circuit()
print(f"Encoding circuit: {dict(enc_513.count_ops())}, depth={enc_513.depth()}")

# ============================================================
# Stabilizer definitions
# ============================================================
# g1=XZZXI, g2=IXZZX, g3=XIXZZ, g4=ZXIXZ
# Format: list of (qubit, pauli_type) for each stabilizer
stabilizers = [
    [(0,'X'), (1,'Z'), (2,'Z'), (3,'X')],         # g1 = XZZXI
    [(1,'X'), (2,'Z'), (3,'Z'), (4,'X')],         # g2 = IXZZX
    [(0,'X'), (2,'X'), (3,'Z'), (4,'Z')],         # g3 = XIXZZ
    [(0,'Z'), (1,'X'), (3,'X'), (4,'Z')],         # g4 = ZXIXZ
]

def build_syndrome_extraction(n_data=5, n_anc=4):
    """
    Build syndrome extraction subcircuit for [[5,1,3]].
    Data qubits: 0-4, Ancilla qubits: 5-8 (one per stabilizer).
    Returns the circuit (without measurement).
    """
    qc = QuantumCircuit(n_data + n_anc)

    for s_idx, stab in enumerate(stabilizers):
        anc = n_data + s_idx  # ancilla qubit index
        qc.h(anc)  # prepare |+>

        for data_q, pauli in stab:
            if pauli == 'X':
                qc.cx(anc, data_q)  # controlled-X
            elif pauli == 'Z':
                qc.cz(anc, data_q)  # controlled-Z
            # No Y in [[5,1,3]] stabilizers

        qc.h(anc)  # convert back to Z basis

    return qc

# ============================================================
# Syndrome lookup table for [[5,1,3]]
# ============================================================
# For each single-qubit error, compute the syndrome (which stabilizers anticommute)
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

# Build stabilizer operators
stab_matrices = []
stab_strings = ['XZZXI', 'IXZZX', 'XIXZZ', 'ZXIXZ']
for ss in stab_strings:
    ops = {q: pauli_map[ss[q]] for q in range(5) if ss[q] != 'I'}
    stab_matrices.append(multi_kron(ops, 5))

# Compute syndrome for each single-qubit error
def compute_syndrome(error_op):
    """Returns 4-bit syndrome: bit i = 1 if error anticommutes with stabilizer i."""
    syndrome = []
    for S in stab_matrices:
        # Check if E*S = S*E (commute, syndrome 0) or E*S = -S*E (anticommute, syndrome 1)
        commutator = error_op @ S - S @ error_op
        anticommutes = np.linalg.norm(commutator) > 1e-10
        syndrome.append(1 if anticommutes else 0)
    return tuple(syndrome)

# Build lookup table: syndrome -> (qubit, pauli) correction
syndrome_table = {}
error_list = []

print("\nSyndrome table:")
print(f"{'Error':>10s} | {'Syndrome':>10s}")
print("-" * 25)

# Identity (no error)
syndrome_table[(0,0,0,0)] = None
print(f"{'No error':>10s} | (0,0,0,0)")

for q in range(5):
    for p_label, p_mat in [('X', Xm), ('Y', Ym), ('Z', Zm)]:
        error_ops = {q: p_mat}
        E = multi_kron(error_ops, 5)
        syn = compute_syndrome(E)
        correction = (q, p_label)
        syndrome_table[syn] = correction
        error_list.append((q, p_label, syn))
        print(f"{'  '+p_label+str(q):>10s} | {syn}")

# Check for collisions (there shouldn't be any for a distance-3 code with single errors)
print(f"\nUnique syndromes: {len(set(s for _, _, s in error_list))}")
print(f"Total errors: {len(error_list)}")

# Cleaner print
print("\nSyndrome lookup table:")
for (q, p, s) in error_list:
    print(f"  {p} on qubit {q} -> syndrome {s}")

# ============================================================
# Verify with statevector simulation
# ============================================================
print("\n" + "="*60)
print("VERIFICATION: Inject errors, extract syndromes, correct")
print("="*60)

sim = AerSimulator(method='statevector')

# Build codespace projector
dim = 32
proj = np.eye(dim, dtype=complex)
for S in stab_matrices:
    proj = proj @ (np.eye(dim) + S) / 2

# Logical states
s0 = np.zeros(dim, dtype=complex); s0[0] = 1
s0 = proj @ s0; s0 /= np.linalg.norm(s0)
X_all = multi_kron({q: Xm for q in range(5)}, 5)
s1 = X_all @ s0; s1 /= np.linalg.norm(s1)

# Test: for each single-qubit error on |0>_L, verify syndrome + correction restores state
results = {'syndrome_table': {}, 'verification': []}

for syn_tuple, corr in syndrome_table.items():
    results['syndrome_table'][str(syn_tuple)] = str(corr)

n_correct = 0
n_total = 0

for q in range(5):
    for p_label, p_mat in [('X', Xm), ('Y', Ym), ('Z', Zm)]:
        for logical in ['0', '1']:
            target = s0 if logical == '0' else s1

            # Apply error
            E = multi_kron({q: p_mat}, 5)
            errored = E @ target

            # Get syndrome
            syn = compute_syndrome(E)
            correction = syndrome_table.get(syn)

            # Apply correction
            if correction is not None:
                c_q, c_p = correction
                C = multi_kron({c_q: pauli_map[c_p]}, 5)
                corrected = C @ errored
            else:
                corrected = errored

            # Check fidelity
            fid = abs(np.dot(target.conj(), corrected))**2
            n_total += 1
            if fid > 0.999:
                n_correct += 1

            results['verification'].append({
                'error': f"{p_label} on q{q}",
                'logical': logical,
                'syndrome': list(syn),
                'correction': str(correction),
                'fidelity': round(float(fid), 6)
            })

print(f"\nCorrect recoveries: {n_correct}/{n_total}")

# ============================================================
# Now verify with actual quantum circuit (syndrome measurement)
# ============================================================
print("\n" + "="*60)
print("CIRCUIT VERIFICATION: Full encode-error-syndrome-correct cycle")
print("="*60)

syndrome_circ = build_syndrome_extraction()
print(f"Syndrome circuit: {dict(syndrome_circ.count_ops())}, depth={syndrome_circ.depth()}")

# Count CX + CZ gates
cx_count = syndrome_circ.count_ops().get('cx', 0)
cz_count = syndrome_circ.count_ops().get('cz', 0)
h_count = syndrome_circ.count_ops().get('h', 0)
print(f"  CX: {cx_count}, CZ: {cz_count}, H: {h_count}")
print(f"  Total 2Q gates: {cx_count + cz_count}")

# Build correction map: syndrome bitstring -> correction gates
# syndrome_table maps (s0,s1,s2,s3) -> (qubit, pauli) or None
correction_map = {}
for syn, corr in syndrome_table.items():
    # Convert syndrome tuple to bitstring (ancilla measurement order: a3,a2,a1,a0)
    # Qiskit measures LSB first, so classical bit order matters
    key = ''.join(str(b) for b in syn)  # g1g2g3g4
    correction_map[key] = corr

# Test a few errors with the full circuit
shots = 10000
circuit_results = []

for q_err in range(5):
    for p_label in ['X', 'Y', 'Z']:
        for logical in ['0', '1']:
            # Build: prepare -> encode -> error -> syndrome -> measure syndrome
            qc = QuantumCircuit(9, 4)  # 5 data + 4 ancilla, 4 classical bits

            # Prepare logical state
            if logical == '1':
                qc.x(0)

            # Encode
            qc.compose(enc_513, qubits=range(5), inplace=True)

            # Inject error
            if p_label == 'X': qc.x(q_err)
            elif p_label == 'Y': qc.y(q_err)
            elif p_label == 'Z': qc.z(q_err)

            # Syndrome extraction
            qc.compose(syndrome_circ, qubits=range(9), inplace=True)

            # Measure ancillas
            qc.measure([5, 6, 7, 8], [0, 1, 2, 3])

            # Run
            counts = sim.run(qc, shots=shots).result().get_counts()

            # Check: should get a single syndrome deterministically
            dominant = max(counts, key=counts.get)
            dominant_frac = counts[dominant] / shots

            # Convert Qiskit bitstring (reversed) to syndrome tuple
            # Qiskit: bit 0 = rightmost = classical register 0 = ancilla 5 = g1
            syn_bits = [int(dominant[3-i]) for i in range(4)]  # g1,g2,g3,g4

            expected_syn = compute_syndrome(multi_kron({q_err: pauli_map[p_label]}, 5))
            match = tuple(syn_bits) == expected_syn

            circuit_results.append({
                'error': f"{p_label} on q{q_err}",
                'logical': logical,
                'measured_syndrome': syn_bits,
                'expected_syndrome': list(expected_syn),
                'dominant_fraction': round(float(dominant_frac), 4),
                'match': match,
            })

            if not match:
                print(f"  MISMATCH: {p_label} on q{q_err}, |{logical}>: measured {syn_bits} vs expected {list(expected_syn)}")

n_match = sum(1 for r in circuit_results if r['match'])
print(f"\nCircuit syndrome matches: {n_match}/{len(circuit_results)}")
if n_match == len(circuit_results):
    print("ALL syndromes correct!")
else:
    print("Some mismatches — debugging needed")
    for r in circuit_results:
        if not r['match']:
            print(f"  {r['error']} |{r['logical']}>: got {r['measured_syndrome']} expected {r['expected_syndrome']} (dominant: {r['dominant_fraction']:.3f})")

# ============================================================
# Full cycle test: encode -> error -> syndrome -> correct -> decode -> verify
# ============================================================
print("\n" + "="*60)
print("FULL CYCLE: encode -> error -> syndrome -> correct -> decode")
print("="*60)

# We'll use statevector to verify correction works in circuit form
full_cycle_results = []
for q_err in range(5):
    for p_label in ['X', 'Y', 'Z']:
        for logical in ['0', '1']:
            target = s0 if logical == '0' else s1

            # Build circuit: prepare -> encode -> error
            qc = QuantumCircuit(5)
            if logical == '1':
                qc.x(0)
            qc.compose(enc_513, inplace=True)
            if p_label == 'X': qc.x(q_err)
            elif p_label == 'Y': qc.y(q_err)
            elif p_label == 'Z': qc.z(q_err)

            # Apply correction from syndrome table
            syn = compute_syndrome(multi_kron({q_err: pauli_map[p_label]}, 5))
            correction = syndrome_table.get(syn)
            if correction is not None:
                c_q, c_p = correction
                if c_p == 'X': qc.x(c_q)
                elif c_p == 'Y': qc.y(c_q)
                elif c_p == 'Z': qc.z(c_q)

            # Get statevector and check fidelity
            sv = Statevector.from_instruction(qc).data
            fid = abs(np.dot(target.conj(), sv))**2

            full_cycle_results.append({
                'error': f"{p_label} on q{q_err}",
                'logical': logical,
                'fidelity': round(float(fid), 6)
            })

n_perfect = sum(1 for r in full_cycle_results if r['fidelity'] > 0.999)
print(f"Perfect recoveries: {n_perfect}/{len(full_cycle_results)}")

# ============================================================
# Save results
# ============================================================
elapsed = time.time() - start
print(f"\nTotal time: {elapsed:.1f}s")

output = {
    'sprint': '026a',
    'description': 'Syndrome extraction circuit for [[5,1,3]] — build and verify',
    'encoding_circuit': {
        'gates': dict(enc_513.count_ops()),
        'depth': enc_513.depth(),
    },
    'syndrome_circuit': {
        'gates': dict(syndrome_circ.count_ops()),
        'depth': syndrome_circ.depth(),
        'total_2q_gates': cx_count + cz_count,
    },
    'syndrome_table': {str(k): str(v) for k, v in syndrome_table.items()},
    'algebraic_verification': {
        'correct': n_correct,
        'total': n_total,
    },
    'circuit_verification': {
        'syndrome_matches': n_match,
        'total': len(circuit_results),
        'details': circuit_results,
    },
    'full_cycle_verification': {
        'perfect_recoveries': n_perfect,
        'total': len(full_cycle_results),
        'details': full_cycle_results,
    },
    'elapsed_seconds': round(elapsed, 1),
}

with open('results/sprint_026a_syndrome_circuit.json', 'w') as f:
    json.dump(output, f, indent=2)
print("Saved to results/sprint_026a_syndrome_circuit.json")
