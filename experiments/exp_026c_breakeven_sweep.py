"""
Sprint 026c: Break-even analysis — at what gate error rate does active correction help?

Sweep p2q from 0.0001 to 0.02. For each error rate, compare:
1. Passive [[5,1,3]] (encode-only)
2. Active [[5,1,3]] (encode + syndrome + correct)

Also test: what if we only correct when syndrome is "confident" (high-weight syndrome)?
And: what about repeated syndrome rounds (majority vote on syndrome)?

Key prediction: active correction needs p2q << 0.008 to help. The break-even
is where the syndrome reliability exceeds 50% for non-trivial syndromes.
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
n = 5
destab_x = np.array([[1,1,1,1,1],[1,0,1,1,1],[0,1,1,0,0],[1,0,0,0,1],[0,1,1,1,1]], dtype=bool)
destab_z = np.zeros((5,5), dtype=bool)
stab_x = np.array([[0,0,0,0,0],[1,0,0,1,0],[0,1,0,0,1],[1,0,1,0,0],[0,1,0,1,0]], dtype=bool)
stab_z = np.array([[1,1,1,1,1],[0,1,1,0,0],[0,0,1,1,0],[0,0,0,1,1],[1,0,0,0,1]], dtype=bool)
tableau = np.zeros((10, 11), dtype=bool)
tableau[:5, :5] = destab_x; tableau[:5, 5:10] = destab_z
tableau[5:, :5] = stab_x;  tableau[5:, 5:10] = stab_z
enc_513 = Clifford(tableau).to_circuit()

# Stabilizers and syndrome table
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

# Syndrome extraction subcircuit
stabilizers = [
    [(0,'X'), (1,'Z'), (2,'Z'), (3,'X')],
    [(1,'X'), (2,'Z'), (3,'Z'), (4,'X')],
    [(0,'X'), (2,'X'), (3,'Z'), (4,'Z')],
    [(0,'Z'), (1,'X'), (3,'X'), (4,'Z')],
]

def build_syndrome_extraction():
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

def holevo_from_err(p):
    p = max(1e-12, min(1 - 1e-12, p))
    return 1.0 + p * np.log2(p) + (1 - p) * np.log2(1 - p)

def decode_parity(bs):
    return sum(int(b) for b in bs) % 2

# ============================================================
# Sweep gate error rate
# ============================================================
p1q_base = 0.0005
p_ro = 0.005  # Use lower readout error to isolate gate effects

# Sweep: p2q from very low to current hardware
p2q_values = [0.0001, 0.0003, 0.0005, 0.001, 0.002, 0.003, 0.005, 0.008, 0.01, 0.015, 0.02]
shots = 30000

pm = generate_preset_pass_manager(optimization_level=2, basis_gates=['cx', 'id', 'rz', 'sx', 'x'])

results = {'sweep': [], 'p_ro': p_ro}

print(f"Sweeping p2q over {len(p2q_values)} values, p_ro={p_ro}")
print(f"{'p2q':>8s} | {'Uncoded':>8s} | {'Passive':>8s} | {'Active':>8s} | {'Act-Pas':>8s} | {'%trivSyn':>8s}")
print("-" * 60)

for p2q in p2q_values:
    t0 = time.time()
    p1q = p2q / 10  # Scale single-qubit error with 2Q error

    noise = NoiseModel()
    noise.add_all_qubit_quantum_error(depolarizing_error(p1q, 1), ['h', 'x', 's', 'sdg', 'z', 'rz', 'sx', 'id', 'y'])
    noise.add_all_qubit_quantum_error(depolarizing_error(p2q, 2), ['cx', 'cz', 'ecr', 'swap'])
    noise.add_all_qubit_readout_error(ReadoutError([[1-p_ro, p_ro], [p_ro, 1-p_ro]]))
    sim = AerSimulator(noise_model=noise)

    row = {'p2q': p2q, 'p1q': p1q}

    # Basis-average over Z, X, Y
    for strategy in ['uncoded', 'passive', 'active']:
        holevos = []
        trivial_frac = 0

        for basis in ['Z', 'X', 'Y']:
            # Build circuits
            if strategy == 'uncoded':
                qc0 = QuantumCircuit(1, 1)
                qc1 = QuantumCircuit(1, 1)
                qc1.x(0)
                if basis == 'Z': pass
                elif basis == 'X': qc0.h(0); qc1.h(0)
                elif basis == 'Y': qc0.h(0); qc0.s(0); qc1.h(0); qc1.s(0)
                if basis == 'X': qc0.h(0); qc1.h(0)
                elif basis == 'Y': qc0.sdg(0); qc0.h(0); qc1.sdg(0); qc1.h(0)
                qc0.measure(0, 0); qc1.measure(0, 0)
                qc0t = pm.run(qc0); qc1t = pm.run(qc1)
                c0 = sim.run(qc0t, shots=shots).result().get_counts()
                c1 = sim.run(qc1t, shots=shots).result().get_counts()
                e0 = sum(v for k, v in c0.items() if int(k) == 1) / shots
                e1 = sum(v for k, v in c1.items() if int(k) == 0) / shots

            elif strategy == 'passive':
                qc0 = QuantumCircuit(5, 5); qc1 = QuantumCircuit(5, 5)
                qc1.x(0)
                if basis == 'X': qc0.h(0); qc1.h(0)
                elif basis == 'Y': qc0.h(0); qc0.s(0); qc1.h(0); qc1.s(0)
                qc0.compose(enc_513, inplace=True); qc1.compose(enc_513, inplace=True)
                if basis == 'X':
                    for q in range(5): qc0.h(q); qc1.h(q)
                elif basis == 'Y':
                    for q in range(5): qc0.sdg(q); qc0.h(q); qc1.sdg(q); qc1.h(q)
                qc0.measure(range(5), range(5)); qc1.measure(range(5), range(5))
                qc0t = pm.run(qc0); qc1t = pm.run(qc1)
                c0 = sim.run(qc0t, shots=shots).result().get_counts()
                c1 = sim.run(qc1t, shots=shots).result().get_counts()
                e0 = sum(v for k, v in c0.items() if decode_parity(k) == 1) / shots
                e1 = sum(v for k, v in c1.items() if decode_parity(k) == 0) / shots

            elif strategy == 'active':
                qc0 = QuantumCircuit(9, 9); qc1 = QuantumCircuit(9, 9)
                qc1.x(0)
                if basis == 'X': qc0.h(0); qc1.h(0)
                elif basis == 'Y': qc0.h(0); qc0.s(0); qc1.h(0); qc1.s(0)
                qc0.compose(enc_513, qubits=range(5), inplace=True)
                qc1.compose(enc_513, qubits=range(5), inplace=True)
                qc0.compose(syndrome_circ, qubits=range(9), inplace=True)
                qc1.compose(syndrome_circ, qubits=range(9), inplace=True)
                qc0.measure([5,6,7,8], [5,6,7,8]); qc1.measure([5,6,7,8], [5,6,7,8])
                if basis == 'X':
                    for q in range(5): qc0.h(q); qc1.h(q)
                elif basis == 'Y':
                    for q in range(5): qc0.sdg(q); qc0.h(q); qc1.sdg(q); qc1.h(q)
                qc0.measure(range(5), range(5)); qc1.measure(range(5), range(5))
                qc0t = pm.run(qc0); qc1t = pm.run(qc1)
                c0 = sim.run(qc0t, shots=shots).result().get_counts()
                c1 = sim.run(qc1t, shots=shots).result().get_counts()

                # Post-process with correction
                def process_active(counts, target_parity, basis):
                    n_err = 0; n_tot = 0; n_trivial = 0
                    for bs, cnt in counts.items():
                        bits = bs[::-1]
                        data = [int(bits[i]) for i in range(5)]
                        syn = tuple(int(bits[i]) for i in range(5, 9))
                        if syn == (0,0,0,0):
                            n_trivial += cnt
                        corr = syndrome_table.get(syn)
                        if corr is not None:
                            cq, cp = corr
                            if basis == 'Z' and cp in ('X','Y'): data[cq] ^= 1
                            elif basis == 'X' and cp in ('Z','Y'): data[cq] ^= 1
                            elif basis == 'Y' and cp in ('X','Z'): data[cq] ^= 1
                        parity = sum(data) % 2
                        if parity != target_parity: n_err += cnt
                        n_tot += cnt
                    return n_err / n_tot, n_trivial / n_tot

                e0, tf0 = process_active(c0, 0, basis)
                e1, tf1 = process_active(c1, 1, basis)
                trivial_frac += (tf0 + tf1) / 2

            avg_err = (e0 + e1) / 2
            holevos.append(holevo_from_err(avg_err))

        avg_h = sum(holevos) / 3
        row[strategy] = round(float(avg_h), 4)
        if strategy == 'active':
            row['trivial_syndrome_frac'] = round(float(trivial_frac / 3), 4)

    delta = row['active'] - row['passive']
    row['active_minus_passive'] = round(float(delta), 4)
    results['sweep'].append(row)

    tf = row.get('trivial_syndrome_frac', 0)
    print(f"{p2q:8.4f} | {row['uncoded']:8.4f} | {row['passive']:8.4f} | {row['active']:8.4f} | {delta:+8.4f} | {tf:8.3f}")

# ============================================================
# Find break-even
# ============================================================
print("\n" + "="*60)
print("BREAK-EVEN ANALYSIS")
print("="*60)

# Find where active > passive
crossover = None
for i in range(len(results['sweep']) - 1):
    d1 = results['sweep'][i]['active_minus_passive']
    d2 = results['sweep'][i+1]['active_minus_passive']
    if d1 >= 0 and d2 < 0:
        # Linear interpolation
        p1 = results['sweep'][i]['p2q']
        p2 = results['sweep'][i+1]['p2q']
        crossover = p1 + (p2 - p1) * d1 / (d1 - d2)
        print(f"Active > Passive crossover at p2q ≈ {crossover:.5f}")
        break
    elif d1 > 0:
        print(f"Active still better at p2q={results['sweep'][i]['p2q']}")

if crossover is None:
    # Check if active is always worse or always better
    all_deltas = [r['active_minus_passive'] for r in results['sweep']]
    if all(d <= 0 for d in all_deltas):
        print("Active correction NEVER beats passive in this range!")
        print(f"Best delta: {max(all_deltas):.4f} at p2q={results['sweep'][np.argmax(all_deltas)]['p2q']}")
    elif all(d >= 0 for d in all_deltas):
        print("Active correction ALWAYS beats passive in this range!")

# Also find where passive > uncoded
crossover_passive = None
for i in range(len(results['sweep']) - 1):
    d1 = results['sweep'][i]['passive'] - results['sweep'][i]['uncoded']
    d2 = results['sweep'][i+1]['passive'] - results['sweep'][i+1]['uncoded']
    if d1 >= 0 and d2 < 0:
        p1 = results['sweep'][i]['p2q']
        p2 = results['sweep'][i+1]['p2q']
        crossover_passive = p1 + (p2 - p1) * d1 / (d1 - d2)
        print(f"Passive > Uncoded crossover at p2q ≈ {crossover_passive:.5f}")
        break

# Summary statistics
print("\nSummary:")
for r in results['sweep']:
    winner = 'uncoded'
    if r['passive'] > r['uncoded']: winner = 'passive'
    if r['active'] > max(r['uncoded'], r['passive']): winner = 'active'
    print(f"  p2q={r['p2q']:.4f}: winner={winner:>8s}  (active-passive={r['active_minus_passive']:+.4f})")

results['crossover_active_passive'] = crossover
results['crossover_passive_uncoded'] = crossover_passive

elapsed = time.time() - start
print(f"\nTotal time: {elapsed:.1f}s")
results['elapsed_seconds'] = round(elapsed, 1)

with open('results/sprint_026c_breakeven_sweep.json', 'w') as f:
    json.dump(results, f, indent=2)
print("Saved to results/sprint_026c_breakeven_sweep.json")
