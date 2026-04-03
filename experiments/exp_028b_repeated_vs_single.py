"""
Sprint 028b: Repeated syndrome (3 rounds, majority vote) vs single-round vs passive.

The critical test: does repeating syndrome measurement 3 times with majority vote
provide enough syndrome reliability for active correction to beat passive encoding?

Strategies compared:
1. Uncoded (1 qubit)
2. Passive [[5,1,3]] (encode-only, parity decode)
3. Active single-round (1 syndrome extraction + correction) — Sprint 026
4. Active 3-round (3 syndrome extractions + majority vote + correction) — NEW

Key parameters:
- Hardware noise: p1q=0.0005, p2q=0.008, p_readout=0.015
- 3-round circuit: 17 qubits, 48 2Q syndrome gates + ~10 encoding gates
- Basis-averaged Holevo as figure of merit
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
# Encoding circuit
# ============================================================
destab_x = np.array([[1,1,1,1,1],[1,0,1,1,1],[0,1,1,0,0],[1,0,0,0,1],[0,1,1,1,1]], dtype=bool)
destab_z = np.zeros((5,5), dtype=bool)
stab_x = np.array([[0,0,0,0,0],[1,0,0,1,0],[0,1,0,0,1],[1,0,1,0,0],[0,1,0,1,0]], dtype=bool)
stab_z = np.array([[1,1,1,1,1],[0,1,1,0,0],[0,0,1,1,0],[0,0,0,1,1],[1,0,0,0,1]], dtype=bool)
tableau = np.zeros((10, 11), dtype=bool)
tableau[:5, :5] = destab_x; tableau[:5, 5:10] = destab_z
tableau[5:, :5] = stab_x;  tableau[5:, 5:10] = stab_z
enc_513 = Clifford(tableau).to_circuit()

# ============================================================
# Stabilizer definitions and syndrome table
# ============================================================
stabilizers = [
    [(0,'X'), (1,'Z'), (2,'Z'), (3,'X')],
    [(1,'X'), (2,'Z'), (3,'Z'), (4,'X')],
    [(0,'X'), (2,'X'), (3,'Z'), (4,'Z')],
    [(0,'Z'), (1,'X'), (3,'X'), (4,'Z')],
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

syndrome_table = {(0,0,0,0): None}
for q in range(5):
    for p_label, p_mat in [('X', Xm), ('Y', Ym), ('Z', Zm)]:
        E = multi_kron({q: p_mat}, 5)
        syn = compute_syndrome(E)
        syndrome_table[syn] = (q, p_label)

# ============================================================
# Noise model
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

def holevo_from_err(p):
    p = max(1e-12, min(1 - 1e-12, p))
    return 1.0 + p * np.log2(p) + (1 - p) * np.log2(1 - p)

def decode_parity(bs):
    return sum(int(b) for b in bs) % 2

# ============================================================
# Circuit builders
# ============================================================
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

def build_active_single(logical, basis):
    """Single-round syndrome extraction (Sprint 026 approach)."""
    qc = QuantumCircuit(9, 9)
    if basis == 'Z':
        if logical == '1': qc.x(0)
    elif basis == 'X':
        if logical == '1': qc.x(0)
        qc.h(0)
    elif basis == 'Y':
        if logical == '1': qc.x(0)
        qc.h(0); qc.s(0)
    qc.compose(enc_513, qubits=range(5), inplace=True)

    # Single round syndrome
    for s_idx, stab in enumerate(stabilizers):
        anc = 5 + s_idx
        qc.h(anc)
        for data_q, pauli in stab:
            if pauli == 'X': qc.cx(anc, data_q)
            elif pauli == 'Z': qc.cz(anc, data_q)
        qc.h(anc)

    # Measure ancillas
    qc.measure([5,6,7,8], [5,6,7,8])

    # Measure data in logical basis
    if basis == 'X':
        for q in range(5): qc.h(q)
    elif basis == 'Y':
        for q in range(5): qc.sdg(q); qc.h(q)
    qc.measure(range(5), range(5))
    return qc

def build_active_repeated(logical, basis, n_rounds=3):
    """Repeated syndrome extraction with fresh ancillas per round."""
    n_total_qubits = 5 + 4 * n_rounds
    n_cl = 5 + 4 * n_rounds  # data + all ancillas
    qc = QuantumCircuit(n_total_qubits, n_cl)

    # Prepare
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

    # Repeated syndrome extraction
    for r in range(n_rounds):
        anc_start = 5 + 4 * r
        for s_idx, stab in enumerate(stabilizers):
            anc = anc_start + s_idx
            qc.h(anc)
            for data_q, pauli in stab:
                if pauli == 'X': qc.cx(anc, data_q)
                elif pauli == 'Z': qc.cz(anc, data_q)
            qc.h(anc)

    # Measure all ancillas
    for r in range(n_rounds):
        for s in range(4):
            qc.measure(5 + 4*r + s, 5 + 4*r + s)

    # Measure data in logical basis
    if basis == 'X':
        for q in range(5): qc.h(q)
    elif basis == 'Y':
        for q in range(5): qc.sdg(q); qc.h(q)
    qc.measure(range(5), range(5))
    return qc

# ============================================================
# Post-processing functions
# ============================================================
def process_single_round(counts, basis, logical_target):
    """Process single-round active correction results."""
    n_total = 0; n_err_corr = 0; n_err_nocorr = 0
    for bitstring, count in counts.items():
        bits = bitstring[::-1]
        data_bits = [int(bits[i]) for i in range(5)]
        syn_bits = tuple(int(bits[i]) for i in range(5, 9))

        # Without correction
        parity_raw = sum(data_bits) % 2
        if parity_raw != logical_target:
            n_err_nocorr += count

        # With correction
        correction = syndrome_table.get(syn_bits)
        corrected_bits = list(data_bits)
        if correction is not None:
            c_q, c_p = correction
            if basis == 'Z' and c_p in ('X', 'Y'): corrected_bits[c_q] ^= 1
            elif basis == 'X' and c_p in ('Z', 'Y'): corrected_bits[c_q] ^= 1
            elif basis == 'Y' and c_p in ('X', 'Z'): corrected_bits[c_q] ^= 1

        parity_corr = sum(corrected_bits) % 2
        if parity_corr != logical_target:
            n_err_corr += count

        n_total += count
    return n_err_corr / n_total, n_err_nocorr / n_total

def process_repeated(counts, basis, logical_target, n_rounds=3):
    """Process repeated syndrome results with majority vote."""
    n_total = 0; n_err_corr = 0; n_err_nocorr = 0
    n_err_last_round = 0
    syndrome_agreement = 0
    total_shots = 0

    for bitstring, count in counts.items():
        bits = bitstring[::-1]
        data_bits = [int(bits[i]) for i in range(5)]

        # Without correction
        parity_raw = sum(data_bits) % 2
        if parity_raw != logical_target:
            n_err_nocorr += count

        # Extract per-round syndromes
        round_syndromes = []
        for r in range(n_rounds):
            syn_r = tuple(int(bits[5 + r*4 + s]) for s in range(4))
            round_syndromes.append(syn_r)

        # Majority vote per stabilizer
        majority_syn = []
        for s in range(4):
            votes = [round_syndromes[r][s] for r in range(n_rounds)]
            majority_syn.append(1 if sum(votes) > n_rounds / 2 else 0)
        majority_syn = tuple(majority_syn)

        # Track agreement
        if all(rs == round_syndromes[0] for rs in round_syndromes):
            syndrome_agreement += count
        total_shots += count

        # With majority vote correction
        correction = syndrome_table.get(majority_syn)
        corrected_bits = list(data_bits)
        if correction is not None:
            c_q, c_p = correction
            if basis == 'Z' and c_p in ('X', 'Y'): corrected_bits[c_q] ^= 1
            elif basis == 'X' and c_p in ('Z', 'Y'): corrected_bits[c_q] ^= 1
            elif basis == 'Y' and c_p in ('X', 'Z'): corrected_bits[c_q] ^= 1
        parity_corr = sum(corrected_bits) % 2
        if parity_corr != logical_target:
            n_err_corr += count

        # Also try last-round-only correction (for comparison)
        last_syn = round_syndromes[-1]
        correction_last = syndrome_table.get(last_syn)
        corrected_last = list(data_bits)
        if correction_last is not None:
            c_q, c_p = correction_last
            if basis == 'Z' and c_p in ('X', 'Y'): corrected_last[c_q] ^= 1
            elif basis == 'X' and c_p in ('Z', 'Y'): corrected_last[c_q] ^= 1
            elif basis == 'Y' and c_p in ('X', 'Z'): corrected_last[c_q] ^= 1
        parity_last = sum(corrected_last) % 2
        if parity_last != logical_target:
            n_err_last_round += count

        n_total += count

    return {
        'err_majority': n_err_corr / n_total,
        'err_nocorr': n_err_nocorr / n_total,
        'err_last_round': n_err_last_round / n_total,
        'syndrome_agreement': syndrome_agreement / total_shots,
    }

# ============================================================
# Run simulations
# ============================================================
print(f"Noise: p1q={p1q}, p2q={p2q}, p_ro={p_ro}, shots={shots}")
results = {}

# 1. UNCODED
print("\n--- Uncoded ---")
results['uncoded'] = {}
for basis in ['Z', 'X', 'Y']:
    qc0 = pm.run(build_uncoded('0', basis))
    qc1 = pm.run(build_uncoded('1', basis))
    c0 = sim.run(qc0, shots=shots).result().get_counts()
    c1 = sim.run(qc1, shots=shots).result().get_counts()
    err0 = sum(v for k, v in c0.items() if int(k) == 1) / shots
    err1 = sum(v for k, v in c1.items() if int(k) == 0) / shots
    avg_err = (err0 + err1) / 2
    h = holevo_from_err(avg_err)
    results['uncoded'][basis] = {'err': round(float(avg_err), 5), 'holevo': round(float(h), 4)}
    print(f"  {basis}: err={avg_err:.4f} Holevo={h:.4f}")

# 2. PASSIVE
print("\n--- Passive [[5,1,3]] ---")
results['passive'] = {}
for basis in ['Z', 'X', 'Y']:
    qc0 = pm.run(build_passive('0', basis))
    qc1 = pm.run(build_passive('1', basis))
    c0 = sim.run(qc0, shots=shots).result().get_counts()
    c1 = sim.run(qc1, shots=shots).result().get_counts()
    err0 = sum(v for k, v in c0.items() if decode_parity(k) == 1) / shots
    err1 = sum(v for k, v in c1.items() if decode_parity(k) == 0) / shots
    avg_err = (err0 + err1) / 2
    h = holevo_from_err(avg_err)
    results['passive'][basis] = {'err': round(float(avg_err), 5), 'holevo': round(float(h), 4)}
    print(f"  {basis}: err={avg_err:.4f} Holevo={h:.4f}")

# 3. ACTIVE SINGLE-ROUND
print("\n--- Active single-round ---")
results['active_1round'] = {}
for basis in ['Z', 'X', 'Y']:
    qc0 = pm.run(build_active_single('0', basis))
    qc1 = pm.run(build_active_single('1', basis))
    c0 = sim.run(qc0, shots=shots).result().get_counts()
    c1 = sim.run(qc1, shots=shots).result().get_counts()
    err_corr_0, err_nocorr_0 = process_single_round(c0, basis, 0)
    err_corr_1, err_nocorr_1 = process_single_round(c1, basis, 1)
    avg_corr = (err_corr_0 + err_corr_1) / 2
    avg_nocorr = (err_nocorr_0 + err_nocorr_1) / 2
    h_corr = holevo_from_err(avg_corr)
    h_nocorr = holevo_from_err(avg_nocorr)
    results['active_1round'][basis] = {
        'err_corrected': round(float(avg_corr), 5),
        'holevo_corrected': round(float(h_corr), 4),
        'err_uncorrected': round(float(avg_nocorr), 5),
        'holevo_uncorrected': round(float(h_nocorr), 4),
    }
    print(f"  {basis}: corrected err={avg_corr:.4f} H={h_corr:.4f} | uncorrected err={avg_nocorr:.4f} H={h_nocorr:.4f}")

# 4. ACTIVE 3-ROUND (the new one)
print("\n--- Active 3-round (majority vote) ---")
results['active_3round'] = {}
for basis in ['Z', 'X', 'Y']:
    t0 = time.time()
    qc0 = pm.run(build_active_repeated('0', basis, n_rounds=3))
    qc1 = pm.run(build_active_repeated('1', basis, n_rounds=3))
    t_transpile = time.time() - t0

    t0 = time.time()
    c0 = sim.run(qc0, shots=shots).result().get_counts()
    c1 = sim.run(qc1, shots=shots).result().get_counts()
    t_sim = time.time() - t0

    res0 = process_repeated(c0, basis, 0, n_rounds=3)
    res1 = process_repeated(c1, basis, 1, n_rounds=3)

    avg_majority = (res0['err_majority'] + res1['err_majority']) / 2
    avg_nocorr = (res0['err_nocorr'] + res1['err_nocorr']) / 2
    avg_last = (res0['err_last_round'] + res1['err_last_round']) / 2
    avg_agreement = (res0['syndrome_agreement'] + res1['syndrome_agreement']) / 2

    h_majority = holevo_from_err(avg_majority)
    h_nocorr = holevo_from_err(avg_nocorr)
    h_last = holevo_from_err(avg_last)

    results['active_3round'][basis] = {
        'err_majority': round(float(avg_majority), 5),
        'holevo_majority': round(float(h_majority), 4),
        'err_nocorr': round(float(avg_nocorr), 5),
        'holevo_nocorr': round(float(h_nocorr), 4),
        'err_last_round': round(float(avg_last), 5),
        'holevo_last_round': round(float(h_last), 4),
        'syndrome_agreement': round(float(avg_agreement), 4),
        'depth': qc0.depth(),
        'transpile_time': round(t_transpile, 1),
        'sim_time': round(t_sim, 1),
    }
    print(f"  {basis}: majority err={avg_majority:.4f} H={h_majority:.4f} | "
          f"last-only err={avg_last:.4f} H={h_last:.4f} | "
          f"nocorr err={avg_nocorr:.4f} H={h_nocorr:.4f}")
    print(f"       syndrome agreement={avg_agreement:.3f} depth={qc0.depth()} "
          f"t_sim={t_sim:.1f}s")

# ============================================================
# Summary
# ============================================================
print(f"\n{'='*70}")
print(f"SUMMARY")
print(f"{'='*70}")

def summarize(label, results_dict, holevo_key='holevo'):
    h_vals = []
    for b in ['Z', 'X', 'Y']:
        h_vals.append(results_dict[b].get(holevo_key, results_dict[b].get('holevo', 0)))
    avg_h = sum(h_vals) / 3
    asym = max(h_vals) - min(h_vals)
    return avg_h, asym, h_vals

for label, rdict, hkey in [
    ('Uncoded', results['uncoded'], 'holevo'),
    ('Passive', results['passive'], 'holevo'),
    ('Active 1-round', results['active_1round'], 'holevo_corrected'),
    ('Active 3-round majority', results['active_3round'], 'holevo_majority'),
    ('Active 3-round last', results['active_3round'], 'holevo_last_round'),
    ('Active 3-round nocorr', results['active_3round'], 'holevo_nocorr'),
]:
    avg_h, asym, h_vals = summarize(label, rdict, hkey)
    print(f"  {label:>25s}: Z={h_vals[0]:.4f} X={h_vals[1]:.4f} Y={h_vals[2]:.4f} | avg={avg_h:.4f} asym={asym:.4f}")

# Deltas
avg_passive, _, _ = summarize('Passive', results['passive'], 'holevo')
avg_1round, _, _ = summarize('1-round', results['active_1round'], 'holevo_corrected')
avg_3round, _, _ = summarize('3-round', results['active_3round'], 'holevo_majority')
avg_uncoded, _, _ = summarize('Uncoded', results['uncoded'], 'holevo')

print(f"\n  Deltas vs passive:")
print(f"    1-round corrected: {avg_1round - avg_passive:+.4f}")
print(f"    3-round majority:  {avg_3round - avg_passive:+.4f}")
print(f"  Deltas vs uncoded:")
print(f"    Passive:           {avg_passive - avg_uncoded:+.4f}")
print(f"    3-round majority:  {avg_3round - avg_uncoded:+.4f}")

# Key question
print(f"\n  KEY QUESTION: Does 3-round majority beat passive?")
if avg_3round > avg_passive:
    print(f"    YES! 3-round majority ({avg_3round:.4f}) > passive ({avg_passive:.4f}) by {avg_3round - avg_passive:+.4f}")
else:
    print(f"    NO. 3-round majority ({avg_3round:.4f}) < passive ({avg_passive:.4f}) by {avg_3round - avg_passive:+.4f}")

elapsed = time.time() - start
print(f"\nTotal time: {elapsed:.1f}s")

# Save
output = {
    'sprint': '028b',
    'description': 'Repeated syndrome (3 rounds, majority vote) vs single-round vs passive',
    'noise_model': {'p1q': p1q, 'p2q': p2q, 'p_readout': p_ro},
    'shots': shots,
    'results': results,
    'summary': {
        'avg_holevo_uncoded': round(float(avg_uncoded), 4),
        'avg_holevo_passive': round(float(avg_passive), 4),
        'avg_holevo_1round': round(float(avg_1round), 4),
        'avg_holevo_3round': round(float(avg_3round), 4),
        'delta_3round_vs_passive': round(float(avg_3round - avg_passive), 4),
        'delta_3round_vs_uncoded': round(float(avg_3round - avg_uncoded), 4),
        'three_round_beats_passive': avg_3round > avg_passive,
    },
    'elapsed_seconds': round(elapsed, 1),
}

with open('results/sprint_028b_repeated_vs_single.json', 'w') as f:
    json.dump(output, f, indent=2)
print("Saved to results/sprint_028b_repeated_vs_single.json")
