"""
Sprint 027a: Build flag-FT syndrome extraction circuit for [[5,1,3]].

Flag qubit approach (Chao & Reichardt 2018):
- For each weight-4 stabilizer, add 1 flag qubit
- Flag CNOT brackets the middle 2 data gates
- If ancilla X error propagates to weight-2 data error, flag detects it
- Modified decoder uses flag info to apply correct weight-2 correction

Circuit layout: 5 data + 4 ancilla + 4 flag = 13 qubits
  data: 0-4, ancilla: 5-8, flag: 9-12

For each stabilizer with gates [g0, g1, g2, g3]:
  H(anc), g0, CX(anc,flag), g1, g2, CX(anc,flag), g3, H(anc)

Weight-2 propagation errors (from ancilla X between g1 and g2):
  g1=XZZXI: Z2X3, syndrome (0,1,0,0), normally decoded as Z4
  g2=IXZZX: Z3X4, syndrome (0,1,1,0), normally decoded as X3
  g3=XIXZZ: Z3Z4, syndrome (1,1,0,1), normally decoded as Y1
  g4=ZXIXZ: X3Z4, syndrome (0,0,1,0), normally decoded as Z2
"""

import numpy as np
import json, time
from qiskit.circuit import QuantumCircuit
from qiskit.quantum_info import Clifford, Statevector
from qiskit_aer import AerSimulator

start = time.time()

# ============================================================
# [[5,1,3]] Encoding (from Sprint 025a)
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
print(f"Encoding circuit: {dict(enc_513.count_ops())}, depth={enc_513.depth()}")

# ============================================================
# Stabilizer definitions
# ============================================================
stabilizers = [
    [(0,'X'), (1,'Z'), (2,'Z'), (3,'X')],   # g1 = XZZXI
    [(1,'X'), (2,'Z'), (3,'Z'), (4,'X')],   # g2 = IXZZX
    [(0,'X'), (2,'X'), (3,'Z'), (4,'Z')],   # g3 = XIXZZ
    [(0,'Z'), (1,'X'), (3,'X'), (4,'Z')],   # g4 = ZXIXZ
]

# Pauli matrices for syndrome computation
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

# Standard syndrome -> correction lookup
syndrome_table = {(0,0,0,0): None}
for q in range(5):
    for p_label, p_mat in [('X', Xm), ('Y', Ym), ('Z', Zm)]:
        E = multi_kron({q: p_mat}, 5)
        syn = compute_syndrome(E)
        syndrome_table[syn] = (q, p_label)

print("\nStandard syndrome table built: 16 entries (1 trivial + 15 errors)")

# ============================================================
# Flag circuit: identify weight-2 propagation errors
# ============================================================
print("\n" + "="*60)
print("WEIGHT-2 PROPAGATION ERROR ANALYSIS")
print("="*60)

# For each stabilizer, the weight-2 error from ancilla X between gates 2 and 3
# (middle of the 4-gate sequence)
flag_weight2_errors = {}

for s_idx, stab in enumerate(stabilizers):
    # Gates are stab[0], stab[1], stab[2], stab[3]
    # Weight-2 error: propagation through gates 2 and 3 only
    # gate[2]: if X -> X on data, if Z -> Z on data (from CZ with X on ancilla)
    # gate[3]: if X -> X on data, if Z -> Z on data
    q2, p2 = stab[2]
    q3, p3 = stab[3]

    # X error on ancilla propagates through CX as X on target, through CZ as Z on target
    err_on_q2 = 'X' if p2 == 'X' else 'Z'  # CX: X->X, CZ: X->Z
    err_on_q3 = 'X' if p3 == 'X' else 'Z'

    # Build the weight-2 error operator
    err_ops = {q2: pauli_map[err_on_q2], q3: pauli_map[err_on_q3]}
    E2 = multi_kron(err_ops, 5)
    syn2 = compute_syndrome(E2)

    # What would standard decoder do?
    std_correction = syndrome_table.get(syn2)

    flag_weight2_errors[s_idx] = {
        'error_desc': f"{err_on_q2}{q2} {err_on_q3}{q3}",
        'error_qubits': [(q2, err_on_q2), (q3, err_on_q3)],
        'syndrome': syn2,
        'std_correction': std_correction,
    }

    print(f"\n  g{s_idx+1} = {''.join(p for _, p in stab).ljust(4, 'I')}")
    print(f"    Weight-2 error: {err_on_q2} on q{q2}, {err_on_q3} on q{q3}")
    print(f"    Syndrome: {syn2}")
    print(f"    Standard correction: {std_correction}")
    print(f"    Standard correction WRONG: would create weight-3 = logical error!")

# ============================================================
# Build flag syndrome extraction circuit
# ============================================================
print("\n" + "="*60)
print("FLAG SYNDROME EXTRACTION CIRCUIT")
print("="*60)

def build_flag_syndrome_extraction(n_data=5, n_anc=4, n_flag=4):
    """
    Build flag-FT syndrome extraction for [[5,1,3]].
    Qubits: 0-4 data, 5-8 ancilla, 9-12 flag.

    For each stabilizer with 4 data gates [g0,g1,g2,g3]:
      H(anc), g0, CX(anc,flag), g1, g2, CX(anc,flag), g3, H(anc)
    """
    total_qubits = n_data + n_anc + n_flag
    qc = QuantumCircuit(total_qubits)

    for s_idx, stab in enumerate(stabilizers):
        anc = n_data + s_idx       # ancilla qubit: 5,6,7,8
        flg = n_data + n_anc + s_idx  # flag qubit: 9,10,11,12

        # Prepare ancilla in |+>
        qc.h(anc)

        # First data gate
        data_q, pauli = stab[0]
        if pauli == 'X': qc.cx(anc, data_q)
        elif pauli == 'Z': qc.cz(anc, data_q)

        # First flag CNOT
        qc.cx(anc, flg)

        # Middle two data gates
        for gate_idx in [1, 2]:
            data_q, pauli = stab[gate_idx]
            if pauli == 'X': qc.cx(anc, data_q)
            elif pauli == 'Z': qc.cz(anc, data_q)

        # Second flag CNOT
        qc.cx(anc, flg)

        # Last data gate
        data_q, pauli = stab[3]
        if pauli == 'X': qc.cx(anc, data_q)
        elif pauli == 'Z': qc.cz(anc, data_q)

        # Convert ancilla back to Z basis
        qc.h(anc)

    return qc

flag_syndrome = build_flag_syndrome_extraction()
gate_counts = dict(flag_syndrome.count_ops())
cx_count = gate_counts.get('cx', 0)
cz_count = gate_counts.get('cz', 0)
h_count = gate_counts.get('h', 0)
total_2q = cx_count + cz_count

print(f"Flag syndrome circuit: {gate_counts}")
print(f"  CX: {cx_count}, CZ: {cz_count}, H: {h_count}")
print(f"  Total 2Q gates: {total_2q} (bare had 16, flag adds {total_2q - 16})")
print(f"  Depth: {flag_syndrome.depth()}")

# ============================================================
# Noiseless verification: inject ancilla errors, check flag
# ============================================================
print("\n" + "="*60)
print("NOISELESS VERIFICATION: Ancilla Error Propagation")
print("="*60)

sim = AerSimulator(method='statevector')
shots = 1000

# Build codespace projector and logical states
dim = 32
proj = np.eye(dim, dtype=complex)
for S in stab_matrices:
    proj = proj @ (np.eye(dim) + S) / 2
s0 = np.zeros(dim, dtype=complex); s0[0] = 1
s0 = proj @ s0; s0 /= np.linalg.norm(s0)
X_all = multi_kron({q: Xm for q in range(5)}, 5)
s1 = X_all @ s0; s1 /= np.linalg.norm(s1)

propagation_tests = []

for s_idx, stab in enumerate(stabilizers):
    anc = 5 + s_idx
    flg = 9 + s_idx

    print(f"\n--- Stabilizer g{s_idx+1}: testing ancilla X error at each position ---")

    # Test X error on ancilla at 5 positions:
    # pos 0: before gate 0 (after H)
    # pos 1: between gate 0 and flag CNOT 1
    # pos 2: between flag CNOT 1 and gate 1 (= between gates in flag bracket)
    # pos 3: between gate 2 and flag CNOT 2
    # pos 4: between flag CNOT 2 and gate 3

    for err_pos in range(5):
        # Build circuit manually: encode |0>_L, then syndrome with ancilla X at specific position
        qc = QuantumCircuit(13, 8)  # 13 qubits, 8 classical (4 syndrome + 4 flag)

        # Encode
        qc.compose(enc_513, qubits=range(5), inplace=True)

        # Syndrome extraction for stabilizer s_idx with X error injection
        a = anc
        f = flg

        qc.h(a)

        gate_counter = 0
        # Position 0: before gate 0
        if err_pos == 0:
            qc.x(a)

        # Gate 0
        dq, pp = stab[0]
        if pp == 'X': qc.cx(a, dq)
        elif pp == 'Z': qc.cz(a, dq)

        # Position 1: between gate 0 and flag CNOT 1
        if err_pos == 1:
            qc.x(a)

        # Flag CNOT 1
        qc.cx(a, f)

        # Position 2: between flag CNOT 1 and gate 1
        if err_pos == 2:
            qc.x(a)

        # Gate 1
        dq, pp = stab[1]
        if pp == 'X': qc.cx(a, dq)
        elif pp == 'Z': qc.cz(a, dq)

        # Gate 2
        dq, pp = stab[2]
        if pp == 'X': qc.cx(a, dq)
        elif pp == 'Z': qc.cz(a, dq)

        # Position 3: between gate 2 and flag CNOT 2
        if err_pos == 3:
            qc.x(a)

        # Flag CNOT 2
        qc.cx(a, f)

        # Position 4: between flag CNOT 2 and gate 3
        if err_pos == 4:
            qc.x(a)

        # Gate 3
        dq, pp = stab[3]
        if pp == 'X': qc.cx(a, dq)
        elif pp == 'Z': qc.cz(a, dq)

        qc.h(a)

        # Also run the other 3 stabilizers (no errors) to get full syndrome
        for s2_idx, stab2 in enumerate(stabilizers):
            if s2_idx == s_idx:
                continue
            a2 = 5 + s2_idx
            f2 = 9 + s2_idx
            qc.h(a2)
            dq, pp = stab2[0]
            if pp == 'X': qc.cx(a2, dq)
            elif pp == 'Z': qc.cz(a2, dq)
            qc.cx(a2, f2)
            dq, pp = stab2[1]
            if pp == 'X': qc.cx(a2, dq)
            elif pp == 'Z': qc.cz(a2, dq)
            dq, pp = stab2[2]
            if pp == 'X': qc.cx(a2, dq)
            elif pp == 'Z': qc.cz(a2, dq)
            qc.cx(a2, f2)
            dq, pp = stab2[3]
            if pp == 'X': qc.cx(a2, dq)
            elif pp == 'Z': qc.cz(a2, dq)
            qc.h(a2)

        # Measure syndrome (ancillas 5-8) and flags (9-12)
        qc.measure([5, 6, 7, 8], [0, 1, 2, 3])  # syndrome bits
        qc.measure([9, 10, 11, 12], [4, 5, 6, 7])  # flag bits

        counts = sim.run(qc, shots=shots).result().get_counts()
        dominant = max(counts, key=counts.get)

        # Parse: Qiskit bitstring is reversed
        bits = dominant[::-1]  # bits[i] = classical register i
        syn_bits = tuple(int(bits[i]) for i in range(4))
        flag_bits = tuple(int(bits[i]) for i in range(4, 8))

        # Determine expected data error weight
        # Error positions and their data error weights:
        pos_labels = [
            "before gate 0 (≡ stabilizer, weight 0 effective)",
            "between gate 0 and flag CNOT 1",
            "between flag CNOT 1 and gate 1 (in bracket)",
            "between gate 2 and flag CNOT 2 (in bracket)",
            "between flag CNOT 2 and gate 3",
        ]

        flag_triggered = flag_bits[s_idx]

        result = {
            'stabilizer': s_idx,
            'error_position': err_pos,
            'position_label': pos_labels[err_pos],
            'syndrome': list(syn_bits),
            'flags': list(flag_bits),
            'flag_triggered': bool(flag_triggered),
        }
        propagation_tests.append(result)

        print(f"  pos {err_pos}: syn={syn_bits} flags={flag_bits} flag_g{s_idx+1}={'YES' if flag_triggered else 'no '} | {pos_labels[err_pos]}")

# ============================================================
# Verify flagged correction: weight-2 errors
# ============================================================
print("\n" + "="*60)
print("FLAGGED CORRECTION VERIFICATION")
print("="*60)

# Build flagged correction lookup
# When flag s triggers and syndrome matches: use weight-2 correction
flag_correction_lookup = {}
for s_idx, info in flag_weight2_errors.items():
    syn = info['syndrome']
    flag_correction_lookup[(s_idx, syn)] = info['error_qubits']
    print(f"  Flag g{s_idx+1} + syndrome {syn} -> correct with {info['error_desc']}")

# Test: for each weight-2 error, verify flagged correction restores state
print("\nTesting weight-2 correction:")
n_correct_flagged = 0
n_correct_standard = 0
n_total_w2 = 0

for s_idx, info in flag_weight2_errors.items():
    for logical_label, target in [('0', s0), ('1', s1)]:
        # Apply the weight-2 error
        err_qubits = info['error_qubits']
        err_ops = {q: pauli_map[p] for q, p in err_qubits}
        E2 = multi_kron(err_ops, 5)
        errored = E2 @ target

        # Standard correction (WRONG)
        std_corr = info['std_correction']
        if std_corr is not None:
            cq, cp = std_corr
            C_std = multi_kron({cq: pauli_map[cp]}, 5)
            std_result = C_std @ errored
        else:
            std_result = errored
        std_fid = abs(np.dot(target.conj(), std_result))**2

        # Flagged correction (RIGHT)
        flag_ops = {q: pauli_map[p] for q, p in err_qubits}
        C_flag = multi_kron(flag_ops, 5)
        flag_result = C_flag @ errored
        flag_fid = abs(np.dot(target.conj(), flag_result))**2

        n_total_w2 += 1
        if std_fid > 0.999: n_correct_standard += 1
        if flag_fid > 0.999: n_correct_flagged += 1

        print(f"  g{s_idx+1} {info['error_desc']} |{logical_label}>: std_fid={std_fid:.4f} flag_fid={flag_fid:.4f}")

print(f"\nStandard correction: {n_correct_standard}/{n_total_w2} correct")
print(f"Flagged correction:  {n_correct_flagged}/{n_total_w2} correct")

# ============================================================
# Full verification: all single-qubit errors (no flag should trigger)
# ============================================================
print("\n" + "="*60)
print("SINGLE-QUBIT ERROR VERIFICATION (flag should NOT trigger)")
print("="*60)

single_qubit_tests = []
for q in range(5):
    for p_label in ['X', 'Y', 'Z']:
        # Encode |0>_L, apply error, run flag syndrome
        qc = QuantumCircuit(13, 8)
        qc.compose(enc_513, qubits=range(5), inplace=True)

        if p_label == 'X': qc.x(q)
        elif p_label == 'Y': qc.y(q)
        elif p_label == 'Z': qc.z(q)

        qc.compose(flag_syndrome, qubits=range(13), inplace=True)
        qc.measure([5,6,7,8], [0,1,2,3])
        qc.measure([9,10,11,12], [4,5,6,7])

        counts = sim.run(qc, shots=shots).result().get_counts()
        dominant = max(counts, key=counts.get)
        bits = dominant[::-1]
        syn_bits = tuple(int(bits[i]) for i in range(4))
        flag_bits = tuple(int(bits[i]) for i in range(4, 8))
        any_flag = any(f == 1 for f in flag_bits)

        # Check syndrome matches expected
        E = multi_kron({q: pauli_map[p_label]}, 5)
        expected_syn = compute_syndrome(E)
        syn_match = syn_bits == expected_syn

        single_qubit_tests.append({
            'error': f"{p_label}{q}",
            'syndrome': list(syn_bits),
            'expected': list(expected_syn),
            'syn_match': syn_match,
            'flags': list(flag_bits),
            'any_flag': any_flag,
        })

n_syn_ok = sum(1 for t in single_qubit_tests if t['syn_match'])
n_flag_ok = sum(1 for t in single_qubit_tests if not t['any_flag'])
print(f"Syndrome matches: {n_syn_ok}/{len(single_qubit_tests)}")
print(f"No flag triggered: {n_flag_ok}/{len(single_qubit_tests)}")
if n_syn_ok < len(single_qubit_tests):
    for t in single_qubit_tests:
        if not t['syn_match']:
            print(f"  MISMATCH: {t['error']}: got {t['syndrome']} expected {t['expected']}")
if n_flag_ok < len(single_qubit_tests):
    for t in single_qubit_tests:
        if t['any_flag']:
            print(f"  FALSE FLAG: {t['error']}: flags={t['flags']}")

# ============================================================
# Save results
# ============================================================
elapsed = time.time() - start
print(f"\nTotal time: {elapsed:.1f}s")

output = {
    'sprint': '027a',
    'description': 'Flag-FT syndrome circuit for [[5,1,3]] — build and verify',
    'flag_circuit': {
        'gates': gate_counts,
        'total_2q_gates': total_2q,
        'extra_2q_vs_bare': total_2q - 16,
        'depth': flag_syndrome.depth(),
        'qubits': 13,
    },
    'weight2_errors': {
        f'g{s+1}': {
            'error': info['error_desc'],
            'syndrome': list(info['syndrome']),
            'std_correction': str(info['std_correction']),
        } for s, info in flag_weight2_errors.items()
    },
    'propagation_tests': propagation_tests,
    'flagged_correction': {
        'standard_correct': n_correct_standard,
        'flagged_correct': n_correct_flagged,
        'total': n_total_w2,
    },
    'single_qubit_tests': {
        'syndrome_matches': n_syn_ok,
        'no_false_flags': n_flag_ok,
        'total': len(single_qubit_tests),
    },
    'elapsed_seconds': round(elapsed, 1),
}

with open('results/sprint_027a_flag_circuit.json', 'w') as f:
    json.dump(output, f, indent=2)
print("Saved to results/sprint_027a_flag_circuit.json")
