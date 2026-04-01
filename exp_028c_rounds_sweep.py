"""
Sprint 028c: Error rate sweep for repeated syndrome measurement.

Key questions from 028b:
1. At what noise level does 3-round majority beat 1-round?
2. Is there ANY noise level where repeated syndrome beats passive?
3. How does the optimal number of rounds depend on noise?

We sweep p2q from 0.0001 to 0.02 and compare:
- Passive (encode-only)
- Active 1-round
- Active 3-round majority vote

Also sweep number of rounds (1, 3, 5) at select noise levels.
Use Z-basis only for speed (isotropy already established).
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
# Encoding + stabilizer setup (compact)
# ============================================================
destab_x = np.array([[1,1,1,1,1],[1,0,1,1,1],[0,1,1,0,0],[1,0,0,0,1],[0,1,1,1,1]], dtype=bool)
destab_z = np.zeros((5,5), dtype=bool)
stab_x = np.array([[0,0,0,0,0],[1,0,0,1,0],[0,1,0,0,1],[1,0,1,0,0],[0,1,0,1,0]], dtype=bool)
stab_z = np.array([[1,1,1,1,1],[0,1,1,0,0],[0,0,1,1,0],[0,0,0,1,1],[1,0,0,0,1]], dtype=bool)
tableau = np.zeros((10, 11), dtype=bool)
tableau[:5, :5] = destab_x; tableau[:5, 5:10] = destab_z
tableau[5:, :5] = stab_x;  tableau[5:, 5:10] = stab_z
enc_513 = Clifford(tableau).to_circuit()

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

pm = generate_preset_pass_manager(optimization_level=2, basis_gates=['cx', 'id', 'rz', 'sx', 'x'])

def holevo_from_err(p):
    p = max(1e-12, min(1 - 1e-12, p))
    return 1.0 + p * np.log2(p) + (1 - p) * np.log2(1 - p)

def decode_parity(bs):
    return sum(int(b) for b in bs) % 2

# ============================================================
# Circuit builders
# ============================================================
def build_passive(logical):
    qc = QuantumCircuit(5, 5)
    if logical == '1': qc.x(0)
    qc.compose(enc_513, inplace=True)
    qc.measure(range(5), range(5))
    return qc

def build_active(logical, n_rounds=1):
    n_qb = 5 + 4 * n_rounds
    n_cl = 5 + 4 * n_rounds
    qc = QuantumCircuit(n_qb, n_cl)
    if logical == '1': qc.x(0)
    qc.compose(enc_513, qubits=range(5), inplace=True)
    for r in range(n_rounds):
        anc_start = 5 + 4 * r
        for s_idx, stab in enumerate(stabilizers):
            anc = anc_start + s_idx
            qc.h(anc)
            for data_q, pauli in stab:
                if pauli == 'X': qc.cx(anc, data_q)
                elif pauli == 'Z': qc.cz(anc, data_q)
            qc.h(anc)
    for r in range(n_rounds):
        for s in range(4):
            qc.measure(5 + 4*r + s, 5 + 4*r + s)
    qc.measure(range(5), range(5))
    return qc

def process_active(counts, n_rounds):
    n_total = 0; n_err_corr = 0; n_err_nocorr = 0
    for bitstring, count in counts.items():
        bits = bitstring[::-1]
        data_bits = [int(bits[i]) for i in range(5)]
        parity_raw = sum(data_bits) % 2

        if n_rounds == 1:
            syn = tuple(int(bits[5 + s]) for s in range(4))
        else:
            round_syns = []
            for r in range(n_rounds):
                syn_r = tuple(int(bits[5 + r*4 + s]) for s in range(4))
                round_syns.append(syn_r)
            syn = tuple(
                1 if sum(round_syns[r][s] for r in range(n_rounds)) > n_rounds/2 else 0
                for s in range(4)
            )

        if parity_raw != 0:  # target is always 0 for |0>_L
            n_err_nocorr += count

        correction = syndrome_table.get(syn)
        corrected_bits = list(data_bits)
        if correction is not None:
            c_q, c_p = correction
            if c_p in ('X', 'Y'): corrected_bits[c_q] ^= 1
        parity_corr = sum(corrected_bits) % 2
        if parity_corr != 0:
            n_err_corr += count
        n_total += count
    return n_err_corr / n_total, n_err_nocorr / n_total

# ============================================================
# Part 1: Error rate sweep for passive, 1-round, 3-round
# ============================================================
print("Part 1: Error rate sweep (Z-basis, |0>_L)")
print("="*60)

p2q_values = [0.0001, 0.0005, 0.001, 0.002, 0.004, 0.008, 0.012, 0.016, 0.02]
shots = 30000

sweep_results = []

# Pre-transpile circuits (they don't depend on noise)
passive_circ = pm.run(build_passive('0'))
active_1r_circ = pm.run(build_active('0', n_rounds=1))
active_3r_circ = pm.run(build_active('0', n_rounds=3))

print(f"Circuit depths: passive={passive_circ.depth()}, 1-round={active_1r_circ.depth()}, 3-round={active_3r_circ.depth()}")
print(f"Circuit 2Q gates: passive={passive_circ.count_ops().get('cx',0)}, "
      f"1-round={active_1r_circ.count_ops().get('cx',0)}, "
      f"3-round={active_3r_circ.count_ops().get('cx',0)}")

for p2q in p2q_values:
    p1q = p2q / 16  # scale single-qubit noise proportionally
    p_ro = max(0.001, p2q * 1.5)  # scale readout

    noise = NoiseModel()
    noise.add_all_qubit_quantum_error(depolarizing_error(p1q, 1), ['h', 'x', 's', 'sdg', 'z', 'rz', 'sx', 'id', 'y'])
    noise.add_all_qubit_quantum_error(depolarizing_error(p2q, 2), ['cx', 'cz', 'ecr', 'swap'])
    noise.add_all_qubit_readout_error(ReadoutError([[1-p_ro, p_ro], [p_ro, 1-p_ro]]))
    sim = AerSimulator(noise_model=noise)

    # Passive
    c_passive = sim.run(passive_circ, shots=shots).result().get_counts()
    err_passive = sum(v for k, v in c_passive.items() if decode_parity(k) == 1) / shots
    h_passive = holevo_from_err(err_passive)

    # Active 1-round
    c_1r = sim.run(active_1r_circ, shots=shots).result().get_counts()
    err_1r_corr, err_1r_nocorr = process_active(c_1r, 1)
    h_1r = holevo_from_err(err_1r_corr)

    # Active 3-round
    c_3r = sim.run(active_3r_circ, shots=shots).result().get_counts()
    err_3r_corr, err_3r_nocorr = process_active(c_3r, 3)
    h_3r = holevo_from_err(err_3r_corr)

    result = {
        'p2q': p2q,
        'p1q': round(p1q, 6),
        'p_ro': round(p_ro, 4),
        'passive_holevo': round(float(h_passive), 4),
        'active_1r_holevo': round(float(h_1r), 4),
        'active_3r_holevo': round(float(h_3r), 4),
        'passive_err': round(float(err_passive), 5),
        'active_1r_err': round(float(err_1r_corr), 5),
        'active_3r_err': round(float(err_3r_corr), 5),
        'delta_1r_vs_passive': round(float(h_1r - h_passive), 4),
        'delta_3r_vs_passive': round(float(h_3r - h_passive), 4),
        'delta_3r_vs_1r': round(float(h_3r - h_1r), 4),
    }
    sweep_results.append(result)

    marker_1r = "BETTER" if h_1r > h_passive else "worse"
    marker_3r = "BETTER" if h_3r > h_passive else "worse"
    print(f"  p2q={p2q:.4f}: passive={h_passive:.4f} 1r={h_1r:.4f}({marker_1r}) "
          f"3r={h_3r:.4f}({marker_3r}) Δ(3r-1r)={h_3r-h_1r:+.4f}")

# ============================================================
# Part 2: Number of rounds sweep at select noise levels
# ============================================================
print(f"\nPart 2: Number of rounds sweep")
print("="*60)

round_values = [1, 3, 5]
test_p2q_values = [0.0001, 0.001, 0.008]

rounds_results = {}

for p2q in test_p2q_values:
    p1q = p2q / 16
    p_ro = max(0.001, p2q * 1.5)

    noise = NoiseModel()
    noise.add_all_qubit_quantum_error(depolarizing_error(p1q, 1), ['h', 'x', 's', 'sdg', 'z', 'rz', 'sx', 'id', 'y'])
    noise.add_all_qubit_quantum_error(depolarizing_error(p2q, 2), ['cx', 'cz', 'ecr', 'swap'])
    noise.add_all_qubit_readout_error(ReadoutError([[1-p_ro, p_ro], [p_ro, 1-p_ro]]))
    sim = AerSimulator(noise_model=noise)

    # Passive baseline
    c_passive = sim.run(passive_circ, shots=shots).result().get_counts()
    err_passive = sum(v for k, v in c_passive.items() if decode_parity(k) == 1) / shots
    h_passive = holevo_from_err(err_passive)

    rounds_data = {'passive': round(float(h_passive), 4)}
    print(f"\n  p2q={p2q}: passive={h_passive:.4f}")

    for n_rounds in round_values:
        n_qb = 5 + 4 * n_rounds
        if n_qb > 25:  # safety limit
            print(f"    {n_rounds} rounds: SKIPPED ({n_qb} qubits)")
            continue

        circ = pm.run(build_active('0', n_rounds=n_rounds))
        t0 = time.time()
        counts = sim.run(circ, shots=shots).result().get_counts()
        t_sim = time.time() - t0

        err_corr, err_nocorr = process_active(counts, n_rounds)
        h_corr = holevo_from_err(err_corr)
        h_nocorr = holevo_from_err(err_nocorr)

        delta = h_corr - h_passive
        marker = "BETTER" if delta > 0 else "worse"
        print(f"    {n_rounds} rounds: corrected={h_corr:.4f} nocorr={h_nocorr:.4f} "
              f"Δ={delta:+.4f} ({marker}) [{t_sim:.1f}s]")

        rounds_data[f'{n_rounds}_rounds_corrected'] = round(float(h_corr), 4)
        rounds_data[f'{n_rounds}_rounds_nocorr'] = round(float(h_nocorr), 4)
        rounds_data[f'{n_rounds}_rounds_delta'] = round(float(delta), 4)
        rounds_data[f'{n_rounds}_rounds_depth'] = circ.depth()
        rounds_data[f'{n_rounds}_rounds_cx'] = circ.count_ops().get('cx', 0)

    rounds_results[str(p2q)] = rounds_data

# ============================================================
# Part 3: Gate overhead analysis
# ============================================================
print(f"\nPart 3: Gate overhead analysis")
print("="*60)

# Count 2Q gates per strategy
gate_counts = {}
for n_rounds in [0, 1, 3, 5]:
    if n_rounds == 0:
        circ = pm.run(build_passive('0'))
        label = 'passive'
    else:
        circ = pm.run(build_active('0', n_rounds=n_rounds))
        label = f'{n_rounds}_round'
    cx = circ.count_ops().get('cx', 0)
    depth = circ.depth()
    gate_counts[label] = {'cx': cx, 'depth': depth}
    print(f"  {label:>10s}: {cx} CX gates, depth {depth}")

# At p2q=0.008, each CX has ~0.8% chance of error
# Expected errors from CX gates alone:
print(f"\n  Expected CX errors at p2q=0.008:")
for label, gc in gate_counts.items():
    expected = gc['cx'] * 0.008
    print(f"    {label:>10s}: {gc['cx']} × 0.008 = {expected:.2f} expected errors")

# ============================================================
# Summary
# ============================================================
print(f"\n{'='*70}")
print("FINAL SUMMARY")
print(f"{'='*70}")

any_3r_beats_passive = any(r['delta_3r_vs_passive'] > 0 for r in sweep_results)
any_1r_beats_passive = any(r['delta_1r_vs_passive'] > 0 for r in sweep_results)
any_3r_beats_1r = any(r['delta_3r_vs_1r'] > 0 for r in sweep_results)

print(f"  3-round ever beats passive? {'YES' if any_3r_beats_passive else 'NO'}")
print(f"  1-round ever beats passive? {'YES' if any_1r_beats_passive else 'NO'}")
print(f"  3-round ever beats 1-round? {'YES' if any_3r_beats_1r else 'NO'}")

# Find best noise level for each strategy
best_delta_3r = max(r['delta_3r_vs_passive'] for r in sweep_results)
best_p2q_3r = [r['p2q'] for r in sweep_results if r['delta_3r_vs_passive'] == best_delta_3r][0]
print(f"  Best 3-round delta vs passive: {best_delta_3r:+.4f} at p2q={best_p2q_3r}")

best_delta_1r = max(r['delta_1r_vs_passive'] for r in sweep_results)
best_p2q_1r = [r['p2q'] for r in sweep_results if r['delta_1r_vs_passive'] == best_delta_1r][0]
print(f"  Best 1-round delta vs passive: {best_delta_1r:+.4f} at p2q={best_p2q_1r}")

elapsed = time.time() - start
print(f"\nTotal time: {elapsed:.1f}s")

output = {
    'sprint': '028c',
    'description': 'Error rate and rounds sweep for repeated syndrome measurement',
    'error_rate_sweep': sweep_results,
    'rounds_sweep': rounds_results,
    'gate_counts': gate_counts,
    'summary': {
        'any_3r_beats_passive': any_3r_beats_passive,
        'any_1r_beats_passive': any_1r_beats_passive,
        'any_3r_beats_1r': any_3r_beats_1r,
        'best_3r_delta': round(float(best_delta_3r), 4),
        'best_3r_p2q': best_p2q_3r,
        'best_1r_delta': round(float(best_delta_1r), 4),
        'best_1r_p2q': best_p2q_1r,
    },
    'elapsed_seconds': round(elapsed, 1),
}

with open('results/sprint_028c_rounds_sweep.json', 'w') as f:
    json.dump(output, f, indent=2)
print("Saved to results/sprint_028c_rounds_sweep.json")
