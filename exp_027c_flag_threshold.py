"""
Sprint 027c: Error rate sweep — where does flag-FT beat passive/bare?

Sweep p2q from 0.0001 to 0.02. Compare:
1. Passive encoding (encode only)
2. Bare syndrome (non-FT correction, Sprint 026)
3. Flag-FT syndrome (flag-aware correction)

Key questions:
- Does flag-FT ever beat passive? (Need the correction to HELP, not just hurt less)
- Does flag-FT show O(p²) scaling vs bare's O(p)?
- What is the crossover noise level between flag-FT and bare?

Also compare bare-without-correction vs bare-with-correction to isolate
the effect of correction from the effect of extra syndrome gates.
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
# Encoding + stabilizers (compact, from 027b)
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

flag_weight2 = {
    0: {'syndrome': (0, 1, 0, 0), 'correction': [(2, 'Z'), (3, 'X')]},
    1: {'syndrome': (1, 0, 1, 0), 'correction': [(3, 'Z'), (4, 'X')]},
    2: {'syndrome': (1, 1, 0, 1), 'correction': [(3, 'Z'), (4, 'Z')]},
    3: {'syndrome': (0, 0, 1, 0), 'correction': [(3, 'X'), (4, 'Z')]},
}

# ============================================================
# Circuit builders (compact)
# ============================================================
def build_bare_syndrome():
    qc = QuantumCircuit(9)
    for s_idx, stab in enumerate(stabilizers):
        anc = 5 + s_idx
        qc.h(anc)
        for dq, pp in stab:
            if pp == 'X': qc.cx(anc, dq)
            elif pp == 'Z': qc.cz(anc, dq)
        qc.h(anc)
    return qc

def build_flag_syndrome():
    qc = QuantumCircuit(13)
    for s_idx, stab in enumerate(stabilizers):
        anc = 5 + s_idx; flg = 9 + s_idx
        qc.h(anc)
        dq, pp = stab[0]
        if pp == 'X': qc.cx(anc, dq)
        elif pp == 'Z': qc.cz(anc, dq)
        qc.cx(anc, flg)
        for gi in [1, 2]:
            dq, pp = stab[gi]
            if pp == 'X': qc.cx(anc, dq)
            elif pp == 'Z': qc.cz(anc, dq)
        qc.cx(anc, flg)
        dq, pp = stab[3]
        if pp == 'X': qc.cx(anc, dq)
        elif pp == 'Z': qc.cz(anc, dq)
        qc.h(anc)
    return qc

bare_synd = build_bare_syndrome()
flag_synd = build_flag_syndrome()

def holevo_from_err(p):
    p = max(1e-12, min(1-1e-12, p))
    return 1.0 + p * np.log2(p) + (1-p) * np.log2(1-p)

def decode_bare(counts, basis, logical_target):
    n_total = 0; n_error = 0
    for bitstring, count in counts.items():
        bits = bitstring[::-1]
        data_bits = [int(bits[i]) for i in range(5)]
        syn_bits = tuple(int(bits[i]) for i in range(5, 9))
        correction = syndrome_table.get(syn_bits)
        if correction is not None:
            cq, cp = correction
            if basis == 'Z' and cp in ('X','Y'): data_bits[cq] ^= 1
            elif basis == 'X' and cp in ('Z','Y'): data_bits[cq] ^= 1
            elif basis == 'Y' and cp in ('X','Z'): data_bits[cq] ^= 1
        parity = sum(data_bits) % 2
        n_total += count
        if parity != logical_target: n_error += count
    return n_error / n_total

def decode_bare_nocorrect(counts, logical_target):
    """Same circuit as bare, but don't apply correction — just decode parity."""
    n_total = 0; n_error = 0
    for bitstring, count in counts.items():
        bits = bitstring[::-1]
        data_bits = [int(bits[i]) for i in range(5)]
        parity = sum(data_bits) % 2
        n_total += count
        if parity != logical_target: n_error += count
    return n_error / n_total

def decode_flag(counts, basis, logical_target):
    n_total = 0; n_error = 0
    n_flagged = 0; n_flag_match = 0
    for bitstring, count in counts.items():
        bits = bitstring[::-1]
        data_bits = [int(bits[i]) for i in range(5)]
        syn_bits = tuple(int(bits[i]) for i in range(5, 9))
        flag_bits = tuple(int(bits[i]) for i in range(9, 13))
        if any(f == 1 for f in flag_bits):
            n_flagged += count

        correction = None; used_flag = False
        for f_idx in range(4):
            if flag_bits[f_idx] == 1 and syn_bits == flag_weight2[f_idx]['syndrome']:
                correction = flag_weight2[f_idx]['correction']
                used_flag = True; n_flag_match += count; break

        if not used_flag:
            std_corr = syndrome_table.get(syn_bits)
            if std_corr is not None: correction = [std_corr]

        if correction is not None:
            for cq, cp in correction:
                if basis == 'Z' and cp in ('X','Y'): data_bits[cq] ^= 1
                elif basis == 'X' and cp in ('Z','Y'): data_bits[cq] ^= 1
                elif basis == 'Y' and cp in ('X','Z'): data_bits[cq] ^= 1

        parity = sum(data_bits) % 2
        n_total += count
        if parity != logical_target: n_error += count
    return n_error / n_total, n_flagged / max(1,n_total), n_flag_match / max(1,n_total)

# ============================================================
# Sweep
# ============================================================
p2q_values = [0.0001, 0.0003, 0.0005, 0.001, 0.002, 0.004, 0.008, 0.012, 0.02]
shots = 30000
pm = generate_preset_pass_manager(optimization_level=2, basis_gates=['cx', 'id', 'rz', 'sx', 'x'])

print(f"Sweeping p2q over {len(p2q_values)} values, {shots} shots each")
print(f"{'p2q':>8s} | {'passive':>8s} | {'bare':>8s} | {'bare_nc':>8s} | {'flag_ft':>8s} | {'flag-bare':>10s} | {'flag-pass':>10s} | {'flag%':>6s}")
print("-" * 90)

sweep_results = []

for p2q in p2q_values:
    t0 = time.time()
    p1q = p2q / 16  # scale 1Q with 2Q
    p_ro = max(0.001, p2q * 2)

    noise = NoiseModel()
    noise.add_all_qubit_quantum_error(depolarizing_error(p1q, 1), ['h','x','s','sdg','z','rz','sx','id','y'])
    noise.add_all_qubit_quantum_error(depolarizing_error(p2q, 2), ['cx','cz','ecr','swap'])
    noise.add_all_qubit_readout_error(ReadoutError([[1-p_ro, p_ro],[p_ro, 1-p_ro]]))
    sim = AerSimulator(noise_model=noise)

    point = {'p2q': p2q, 'p1q': p1q, 'p_ro': p_ro}

    # Build all circuits for all bases and logical states
    for strategy in ['passive', 'bare', 'flag_ft']:
        basis_holevos = []
        total_flag_rate = 0; total_match_rate = 0

        for basis in ['Z', 'X', 'Y']:
            for logical in ['0', '1']:
                if strategy == 'passive':
                    qc = QuantumCircuit(5, 5)
                    if logical == '1': qc.x(0)
                    if basis == 'X': qc.h(0)
                    elif basis == 'Y': qc.h(0); qc.s(0)
                    qc.compose(enc_513, inplace=True)
                    if basis == 'X':
                        for q in range(5): qc.h(q)
                    elif basis == 'Y':
                        for q in range(5): qc.sdg(q); qc.h(q)
                    qc.measure(range(5), range(5))
                elif strategy == 'bare':
                    qc = QuantumCircuit(9, 9)
                    if logical == '1': qc.x(0)
                    if basis == 'X': qc.h(0)
                    elif basis == 'Y': qc.h(0); qc.s(0)
                    qc.compose(enc_513, qubits=range(5), inplace=True)
                    qc.compose(bare_synd, qubits=range(9), inplace=True)
                    qc.measure([5,6,7,8], [5,6,7,8])
                    if basis == 'X':
                        for q in range(5): qc.h(q)
                    elif basis == 'Y':
                        for q in range(5): qc.sdg(q); qc.h(q)
                    qc.measure(range(5), range(5))
                elif strategy == 'flag_ft':
                    qc = QuantumCircuit(13, 13)
                    if logical == '1': qc.x(0)
                    if basis == 'X': qc.h(0)
                    elif basis == 'Y': qc.h(0); qc.s(0)
                    qc.compose(enc_513, qubits=range(5), inplace=True)
                    qc.compose(flag_synd, qubits=range(13), inplace=True)
                    qc.measure([5,6,7,8], [5,6,7,8])
                    qc.measure([9,10,11,12], [9,10,11,12])
                    if basis == 'X':
                        for q in range(5): qc.h(q)
                    elif basis == 'Y':
                        for q in range(5): qc.sdg(q); qc.h(q)
                    qc.measure(range(5), range(5))

                qct = pm.run(qc)
                counts = sim.run(qct, shots=shots).result().get_counts()

                # Store counts for decoding
                if logical == '0':
                    counts0 = counts
                else:
                    counts1 = counts

            # Decode per basis
            if strategy == 'passive':
                err0 = sum(v for k,v in counts0.items() if sum(int(b) for b in k) % 2 == 1) / shots
                err1 = sum(v for k,v in counts1.items() if sum(int(b) for b in k) % 2 == 0) / shots
            elif strategy == 'bare':
                err0 = decode_bare(counts0, basis, 0)
                err1 = decode_bare(counts1, basis, 1)
            elif strategy == 'flag_ft':
                err0, fr0, fm0 = decode_flag(counts0, basis, 0)
                err1, fr1, fm1 = decode_flag(counts1, basis, 1)
                total_flag_rate += (fr0 + fr1) / 2
                total_match_rate += (fm0 + fm1) / 2

            avg_err = (err0 + err1) / 2
            h = holevo_from_err(avg_err)
            basis_holevos.append(h)

        avg_h = sum(basis_holevos) / 3
        point[strategy] = round(float(avg_h), 4)
        if strategy == 'flag_ft':
            point['flag_rate'] = round(total_flag_rate / 3, 4)
            point['match_rate'] = round(total_match_rate / 3, 5)

    # Also compute bare without correction (isolate syndrome circuit damage)
    bare_nc_holevos = []
    for basis in ['Z', 'X', 'Y']:
        for logical in ['0', '1']:
            qc = QuantumCircuit(9, 9)
            if logical == '1': qc.x(0)
            if basis == 'X': qc.h(0)
            elif basis == 'Y': qc.h(0); qc.s(0)
            qc.compose(enc_513, qubits=range(5), inplace=True)
            qc.compose(bare_synd, qubits=range(9), inplace=True)
            qc.measure([5,6,7,8], [5,6,7,8])
            if basis == 'X':
                for q in range(5): qc.h(q)
            elif basis == 'Y':
                for q in range(5): qc.sdg(q); qc.h(q)
            qc.measure(range(5), range(5))
            qct = pm.run(qc)
            counts = sim.run(qct, shots=shots).result().get_counts()
            if logical == '0': counts0 = counts
            else: counts1 = counts
        err0 = decode_bare_nocorrect(counts0, 0)
        err1 = decode_bare_nocorrect(counts1, 1)
        bare_nc_holevos.append(holevo_from_err((err0+err1)/2))
    point['bare_nocorrect'] = round(sum(bare_nc_holevos)/3, 4)

    dt = time.time() - t0
    delta_fb = point['flag_ft'] - point['bare']
    delta_fp = point['flag_ft'] - point['passive']
    fr = point.get('flag_rate', 0)
    print(f"{p2q:8.4f} | {point['passive']:8.4f} | {point['bare']:8.4f} | {point['bare_nocorrect']:8.4f} | {point['flag_ft']:8.4f} | {delta_fb:+10.4f} | {delta_fp:+10.4f} | {fr:5.1%} | {dt:.0f}s")

    sweep_results.append(point)

# ============================================================
# Analysis
# ============================================================
print(f"\n{'='*70}")
print("ANALYSIS")
print(f"{'='*70}")

# Find crossovers
for i in range(len(sweep_results)-1):
    r = sweep_results[i]; r2 = sweep_results[i+1]
    # Flag-FT vs bare crossover
    d1 = r['flag_ft'] - r['bare']
    d2 = r2['flag_ft'] - r2['bare']
    if d1 * d2 < 0:
        # Linear interpolation
        frac = d1 / (d1 - d2)
        p_cross = r['p2q'] + frac * (r2['p2q'] - r['p2q'])
        print(f"  Flag-FT = Bare crossover at p2q ≈ {p_cross:.4f}")

    # Flag-FT vs passive
    d1 = r['flag_ft'] - r['passive']
    d2 = r2['flag_ft'] - r2['passive']
    if d1 * d2 < 0:
        frac = d1 / (d1 - d2)
        p_cross = r['p2q'] + frac * (r2['p2q'] - r['p2q'])
        print(f"  Flag-FT = Passive crossover at p2q ≈ {p_cross:.4f}")

    # Bare vs passive
    d1 = r['bare'] - r['passive']
    d2 = r2['bare'] - r2['passive']
    if d1 * d2 < 0:
        frac = d1 / (d1 - d2)
        p_cross = r['p2q'] + frac * (r2['p2q'] - r['p2q'])
        print(f"  Bare = Passive crossover at p2q ≈ {p_cross:.4f}")

# Correction benefit: does correction ever help?
print("\n  Does correction help? (bare_corrected vs bare_nocorrect):")
for r in sweep_results:
    delta = r['bare'] - r['bare_nocorrect']
    print(f"    p2q={r['p2q']:.4f}: correction delta = {delta:+.4f}")

# Flag improvement over bare
print("\n  Flag-FT improvement over bare:")
for r in sweep_results:
    delta = r['flag_ft'] - r['bare']
    print(f"    p2q={r['p2q']:.4f}: {delta:+.4f}")

elapsed = time.time() - start
print(f"\nTotal time: {elapsed:.1f}s")

output = {
    'sprint': '027c',
    'description': 'Error rate sweep: flag-FT vs bare vs passive',
    'shots': shots,
    'sweep': sweep_results,
    'elapsed_seconds': round(elapsed, 1),
}

with open('results/sprint_027c_flag_threshold.json', 'w') as f:
    json.dump(output, f, indent=2)
print("Saved to results/sprint_027c_flag_threshold.json")
