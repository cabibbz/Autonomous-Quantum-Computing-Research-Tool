"""
Sprint 026b: Active syndrome correction vs passive encoding under gate-level noise.

Compare three strategies for [[5,1,3]]:
1. Uncoded (1 qubit, no protection)
2. Passive encoding (encode -> measure, no syndrome extraction)
3. Active correction (encode -> syndrome -> correct -> measure)

Uses Sprint 025a noise model (IBM hardware-calibrated).
Basis-averaged Holevo information as figure of merit.

Key question: Does syndrome extraction + correction improve over encode-only,
despite adding 16 more 2Q gates?
"""

import numpy as np
import json, time
from qiskit.circuit import QuantumCircuit
from qiskit.quantum_info import Clifford
from qiskit_aer import AerSimulator
from qiskit_aer.noise import NoiseModel, depolarizing_error, ReadoutError
from qiskit.transpiler import generate_preset_pass_manager

start = time.time()

# ============================================================
# Encoding circuit (from Sprint 025a)
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
# Stabilizer definitions and syndrome table (from 026a)
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
# Circuit builders
# ============================================================
def build_syndrome_extraction():
    """Syndrome extraction subcircuit for 5 data + 4 ancilla qubits."""
    qc = QuantumCircuit(9)
    for s_idx, stab in enumerate(stabilizers):
        anc = 5 + s_idx
        qc.h(anc)
        for data_q, pauli in stab:
            if pauli == 'X': qc.cx(anc, data_q)
            elif pauli == 'Z': qc.cz(anc, data_q)
        qc.h(anc)
    return qc

syndrome_circ = build_syndrome_extraction()

def build_uncoded(logical, basis):
    qc = QuantumCircuit(1, 1)
    if basis == 'Z':
        if logical == '1': qc.x(0)
    elif basis == 'X':
        if logical == '1': qc.x(0)
        qc.h(0)
    elif basis == 'Y':
        if logical == '1': qc.x(0)
        qc.h(0); qc.s(0)
    if basis == 'X': qc.h(0)
    elif basis == 'Y': qc.sdg(0); qc.h(0)
    qc.measure(0, 0)
    return qc

def build_passive(logical, basis):
    """Encode-only [[5,1,3]]: prepare, encode, measure in logical basis."""
    qc = QuantumCircuit(5, 5)
    if basis == 'Z':
        if logical == '1': qc.x(0)
    elif basis == 'X':
        if logical == '1': qc.x(0)
        qc.h(0)
    elif basis == 'Y':
        if logical == '1': qc.x(0)
        qc.h(0); qc.s(0)
    qc.compose(enc_513, inplace=True)
    if basis == 'X':
        for q in range(5): qc.h(q)
    elif basis == 'Y':
        for q in range(5): qc.sdg(q); qc.h(q)
    qc.measure(range(5), range(5))
    return qc

def build_active(logical, basis):
    """
    Encode + syndrome extraction + conditional correction + measure.

    Since Qiskit Aer supports mid-circuit measurement and classical conditioning,
    we use: encode -> syndrome measure -> classically-conditioned correction -> decode+measure.

    However, classical conditioning on 4-bit syndrome with lookup table is complex in Qiskit.
    Instead, we use a simpler approach:
    1. Run the circuit up to syndrome measurement
    2. Post-process: for each shot, read syndrome, determine correction,
       and adjust the final data measurement accordingly.

    Actually, for a proper comparison, we should do the full circuit.
    Let's use a hybrid approach: run syndrome measurement, get the syndrome,
    apply correction classically to the final bitstrings.

    The circuit: encode -> syndrome -> measure ancillas AND data qubits.
    Post-processing: for each shot, use syndrome to determine if a correction
    would have been applied, then adjust the decoded logical outcome.
    """
    qc = QuantumCircuit(9, 9)  # 5 data + 4 ancilla, 9 classical

    # Prepare logical state
    if basis == 'Z':
        if logical == '1': qc.x(0)
    elif basis == 'X':
        if logical == '1': qc.x(0)
        qc.h(0)
    elif basis == 'Y':
        if logical == '1': qc.x(0)
        qc.h(0); qc.s(0)

    # Encode
    qc.compose(enc_513, qubits=range(5), inplace=True)

    # Syndrome extraction
    qc.compose(syndrome_circ, qubits=range(9), inplace=True)

    # Measure ancillas (syndrome bits)
    qc.measure([5, 6, 7, 8], [5, 6, 7, 8])

    # Measure data qubits in logical basis
    if basis == 'X':
        for q in range(5): qc.h(q)
    elif basis == 'Y':
        for q in range(5): qc.sdg(q); qc.h(q)
    qc.measure(range(5), range(5))

    return qc

# ============================================================
# Noise model (IBM hardware-calibrated from Sprint 025)
# ============================================================
p1q = 0.0005
p2q = 0.008
p_ro = 0.015

noise = NoiseModel()
noise.add_all_qubit_quantum_error(depolarizing_error(p1q, 1), ['h', 'x', 's', 'sdg', 'z', 'rz', 'sx', 'id', 'y'])
noise.add_all_qubit_quantum_error(depolarizing_error(p2q, 2), ['cx', 'cz', 'ecr', 'swap'])
noise.add_all_qubit_readout_error(ReadoutError([[1-p_ro, p_ro], [p_ro, 1-p_ro]]))

sim = AerSimulator(noise_model=noise)
pm = generate_preset_pass_manager(optimization_level=2, basis_gates=['cx', 'id', 'rz', 'sx', 'x'])
shots = 50000

def decode_parity(bs):
    return sum(int(b) for b in bs) % 2

def holevo_from_err(p):
    p = max(1e-12, min(1 - 1e-12, p))
    return 1.0 + p * np.log2(p) + (1 - p) * np.log2(1 - p)

# ============================================================
# Run simulations
# ============================================================
print("Running noisy simulations...")
print(f"Noise: p1q={p1q}, p2q={p2q}, p_ro={p_ro}")
print(f"Shots: {shots}\n")

results = {}

# 1. UNCODED
print("--- Uncoded ---")
results['uncoded'] = {}
for basis in ['Z', 'X', 'Y']:
    qc0 = build_uncoded('0', basis)
    qc1 = build_uncoded('1', basis)
    qc0t = pm.run(qc0); qc1t = pm.run(qc1)
    counts0 = sim.run(qc0t, shots=shots).result().get_counts()
    counts1 = sim.run(qc1t, shots=shots).result().get_counts()
    err0 = sum(v for k, v in counts0.items() if int(k) == 1) / shots
    err1 = sum(v for k, v in counts1.items() if int(k) == 0) / shots
    avg_err = (err0 + err1) / 2
    h = holevo_from_err(avg_err)
    results['uncoded'][basis] = {'err': round(float(avg_err), 5), 'holevo': round(float(h), 4),
                                  'depth': qc0t.depth()}
    print(f"  {basis}: err={avg_err:.4f} Holevo={h:.4f} depth={qc0t.depth()}")

# 2. PASSIVE [[5,1,3]]
print("\n--- Passive [[5,1,3]] (encode-only) ---")
results['passive'] = {}
for basis in ['Z', 'X', 'Y']:
    qc0 = build_passive('0', basis)
    qc1 = build_passive('1', basis)
    qc0t = pm.run(qc0); qc1t = pm.run(qc1)
    counts0 = sim.run(qc0t, shots=shots).result().get_counts()
    counts1 = sim.run(qc1t, shots=shots).result().get_counts()
    # Decode by parity
    err0 = sum(v for k, v in counts0.items() if decode_parity(k) == 1) / shots
    err1 = sum(v for k, v in counts1.items() if decode_parity(k) == 0) / shots
    avg_err = (err0 + err1) / 2
    h = holevo_from_err(avg_err)
    results['passive'][basis] = {'err': round(float(avg_err), 5), 'holevo': round(float(h), 4),
                                  'depth': qc0t.depth()}
    print(f"  {basis}: err={avg_err:.4f} Holevo={h:.4f} depth={qc0t.depth()}")

# 3. ACTIVE [[5,1,3]] (encode + syndrome + correction)
print("\n--- Active [[5,1,3]] (encode + syndrome + correct) ---")
results['active'] = {}
for basis in ['Z', 'X', 'Y']:
    qc0 = build_active('0', basis)
    qc1 = build_active('1', basis)
    qc0t = pm.run(qc0); qc1t = pm.run(qc1)

    d0 = qc0t.depth()
    counts0 = sim.run(qc0t, shots=shots).result().get_counts()
    counts1 = sim.run(qc1t, shots=shots).result().get_counts()

    # Post-process: for each shot, extract syndrome and data bits
    # Apply correction in classical post-processing
    # Bitstring format: '876543210' where 0-4=data, 5-8=syndrome
    def decode_active_shot(bitstring, syndrome_table):
        """Given a 9-bit measurement, extract syndrome, determine correction, apply to data."""
        # Qiskit bitstring is reversed: bit 0 = rightmost
        bits = bitstring[::-1]  # now bits[i] = classical register i
        data_bits = [int(bits[i]) for i in range(5)]
        syn_bits = tuple(int(bits[i]) for i in range(5, 9))

        correction = syndrome_table.get(syn_bits)
        if correction is not None:
            c_q, c_p = correction
            # Correction in Z-basis measurement: X or Y flip bit, Z doesn't
            # But we're measuring in the logical basis already
            # The correction should be applied BEFORE measurement in the circuit
            # Since we can't do that classically, we need to account for how
            # the correction would affect the measurement outcome.
            #
            # For Z-basis measurement (measuring ZZZZZ):
            #   X correction on qubit q flips data_bits[q]
            #   Y correction on qubit q flips data_bits[q]
            #   Z correction on qubit q does NOT flip data_bits[q]
            #
            # For X-basis measurement (after H on all, measuring XXXXX->ZZZZZ):
            #   Z correction on qubit q flips data_bits[q] (Z before H = X after H)
            #   Y correction on qubit q flips data_bits[q]
            #   X correction on qubit q does NOT flip (X before H = Z after H)
            #
            # For Y-basis (after S†H, measuring YYYYY->ZZZZZ):
            #   All corrections flip (Y anticommutes with all single-qubit Paulis? No.)
            #   X on q: Y measurement = S†H X S†H†... this is complex
            #
            # Actually, simpler approach: the correction C acts on the data state.
            # We measure in basis B. The logical value = parity of all data bits.
            # If correction C anticommutes with the measurement operator on qubit q,
            # it flips that bit. The measurement operator per qubit is Z (after basis change).
            # So:
            #   Z-basis: meas_op = Z per qubit. X anticommutes -> flips. Z commutes -> no flip. Y anticommutes -> flips.
            #   X-basis: after H, meas_op = X per qubit in original frame. X commutes. Z anticommutes -> flips. Y anticommutes -> flips.
            #   Y-basis: after S†H, meas_op = Y per qubit in original frame. X anticommutes -> flips. Z anticommutes -> flips. Y commutes -> no flip.
            pass  # We'll handle this below

        return data_bits, syn_bits, correction

    # Process counts with basis-dependent correction
    def process_counts(counts, basis, logical_target):
        """Process measurement counts with syndrome-based correction."""
        n_total = 0
        n_error = 0
        syn_hist = {}

        for bitstring, count in counts.items():
            bits = bitstring[::-1]
            data_bits = [int(bits[i]) for i in range(5)]
            syn_bits = tuple(int(bits[i]) for i in range(5, 9))

            # Track syndrome distribution
            syn_key = str(syn_bits)
            syn_hist[syn_key] = syn_hist.get(syn_key, 0) + count

            correction = syndrome_table.get(syn_bits)

            # Apply correction classically: flip the appropriate data bit
            if correction is not None:
                c_q, c_p = correction
                # Does this correction flip the measurement outcome on qubit c_q?
                if basis == 'Z':
                    # Measuring Z. X and Y flip, Z doesn't.
                    if c_p in ('X', 'Y'):
                        data_bits[c_q] ^= 1
                elif basis == 'X':
                    # Measuring X (after H). Z and Y flip, X doesn't.
                    if c_p in ('Z', 'Y'):
                        data_bits[c_q] ^= 1
                elif basis == 'Y':
                    # Measuring Y (after S†H). X and Z flip, Y doesn't.
                    if c_p in ('X', 'Z'):
                        data_bits[c_q] ^= 1

            # Decode: parity of data bits
            parity = sum(data_bits) % 2
            n_total += count
            if parity != logical_target:
                n_error += count

        return n_error / n_total, syn_hist

    err0, syn0 = process_counts(counts0, basis, 0)
    err1, syn1 = process_counts(counts1, basis, 1)
    avg_err = (err0 + err1) / 2
    h = holevo_from_err(avg_err)

    # Also compute WITHOUT correction for comparison
    def process_no_correction(counts, logical_target):
        n_total = 0; n_error = 0
        for bitstring, count in counts.items():
            bits = bitstring[::-1]
            data_bits = [int(bits[i]) for i in range(5)]
            parity = sum(data_bits) % 2
            n_total += count
            if parity != logical_target:
                n_error += count
        return n_error / n_total

    err0_nc = process_no_correction(counts0, 0)
    err1_nc = process_no_correction(counts1, 1)
    avg_err_nc = (err0_nc + err1_nc) / 2
    h_nc = holevo_from_err(avg_err_nc)

    # Syndrome distribution for |0>_L
    top_syns = sorted(syn0.items(), key=lambda x: -x[1])[:5]

    results['active'][basis] = {
        'err_corrected': round(float(avg_err), 5),
        'holevo_corrected': round(float(h), 4),
        'err_uncorrected': round(float(avg_err_nc), 5),
        'holevo_uncorrected': round(float(h_nc), 4),
        'depth': d0,
        'top_syndromes': {k: v for k, v in top_syns},
    }
    print(f"  {basis}: corrected err={avg_err:.4f} Holevo={h:.4f} | uncorrected err={avg_err_nc:.4f} Holevo={h_nc:.4f} | depth={d0}")
    print(f"       Top syndromes: {top_syns[:3]}")

# ============================================================
# Summary
# ============================================================
print(f"\n{'='*70}")
print(f"SUMMARY: Active vs Passive QEC (p1q={p1q}, p2q={p2q}, ro={p_ro})")
print(f"{'='*70}")

for strategy in ['uncoded', 'passive', 'active']:
    r = results[strategy]
    if strategy == 'active':
        h_vals = [r[b]['holevo_corrected'] for b in ['Z', 'X', 'Y']]
    else:
        h_vals = [r[b]['holevo'] for b in ['Z', 'X', 'Y']]
    avg_h = sum(h_vals) / 3
    asym = max(h_vals) - min(h_vals)
    min_h = min(h_vals)
    results[strategy]['avg_holevo'] = round(float(avg_h), 4)
    results[strategy]['asymmetry'] = round(float(asym), 4)
    results[strategy]['min_holevo'] = round(float(min_h), 4)
    print(f"\n  {strategy:>10s}: Z={h_vals[0]:.4f}  X={h_vals[1]:.4f}  Y={h_vals[2]:.4f}")
    print(f"             avg={avg_h:.4f}  min={min_h:.4f}  asym={asym:.4f}")

# Correction benefit
if 'active' in results:
    h_active = results['active']['avg_holevo']
    h_passive = results['passive']['avg_holevo']
    h_uncoded = results['uncoded']['avg_holevo']
    print(f"\n  Correction benefit (active - passive): {h_active - h_passive:+.4f}")
    print(f"  Active vs uncoded: {h_active - h_uncoded:+.4f}")
    print(f"  Passive vs uncoded: {h_passive - h_uncoded:+.4f}")

    # Syndrome correction statistics
    for basis in ['Z', 'X', 'Y']:
        r = results['active'][basis]
        delta = r['holevo_corrected'] - r['holevo_uncorrected']
        print(f"  {basis}-basis correction delta: {delta:+.4f}")

elapsed = time.time() - start
print(f"\nTotal time: {elapsed:.1f}s")

output = {
    'sprint': '026b',
    'description': 'Active syndrome correction vs passive encoding under gate-level noise',
    'noise_model': {'p1q': p1q, 'p2q': p2q, 'p_readout': p_ro},
    'shots': shots,
    'results': results,
    'elapsed_seconds': round(elapsed, 1),
}

with open('results/sprint_026b_active_vs_passive.json', 'w') as f:
    json.dump(output, f, indent=2)
print("Saved to results/sprint_026b_active_vs_passive.json")
