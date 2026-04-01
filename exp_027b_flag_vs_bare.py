"""
Sprint 027b: Flag-FT vs bare syndrome vs passive encoding under noise.

Compare four strategies for [[5,1,3]]:
1. Uncoded (1 qubit)
2. Passive encoding (encode → measure, no syndrome)
3. Bare syndrome (encode → syndrome → correct, Sprint 026 approach)
4. Flag-FT syndrome (encode → flag syndrome → flag-aware correct)

Key question: Does flag information overcome the extra gate overhead
(24 vs 16 2Q gates in syndrome circuit)?

Uses IBM hardware-calibrated noise model from Sprint 025.
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
# [[5,1,3]] Encoding (from Sprint 025a)
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
# Stabilizer definitions and syndrome tables
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

# Standard syndrome -> correction lookup
syndrome_table = {(0,0,0,0): None}
for q in range(5):
    for p_label, p_mat in [('X', Xm), ('Y', Ym), ('Z', Zm)]:
        E = multi_kron({q: p_mat}, 5)
        syn = compute_syndrome(E)
        syndrome_table[syn] = (q, p_label)

# Flag-aware correction: weight-2 errors from propagation (from 027a)
# When flag s triggers AND syndrome matches, use weight-2 correction
flag_weight2 = {
    0: {'syndrome': (0, 1, 0, 0), 'correction': [(2, 'Z'), (3, 'X')]},  # g1: Z2X3
    1: {'syndrome': (1, 0, 1, 0), 'correction': [(3, 'Z'), (4, 'X')]},  # g2: Z3X4
    2: {'syndrome': (1, 1, 0, 1), 'correction': [(3, 'Z'), (4, 'Z')]},  # g3: Z3Z4
    3: {'syndrome': (0, 0, 1, 0), 'correction': [(3, 'X'), (4, 'Z')]},  # g4: X3Z4
}

# ============================================================
# Circuit builders
# ============================================================
def build_bare_syndrome():
    """Bare (non-FT) syndrome extraction: 5 data + 4 ancilla = 9 qubits."""
    qc = QuantumCircuit(9)
    for s_idx, stab in enumerate(stabilizers):
        anc = 5 + s_idx
        qc.h(anc)
        for data_q, pauli in stab:
            if pauli == 'X': qc.cx(anc, data_q)
            elif pauli == 'Z': qc.cz(anc, data_q)
        qc.h(anc)
    return qc

def build_flag_syndrome():
    """Flag-FT syndrome extraction: 5 data + 4 ancilla + 4 flag = 13 qubits."""
    qc = QuantumCircuit(13)
    for s_idx, stab in enumerate(stabilizers):
        anc = 5 + s_idx
        flg = 9 + s_idx
        qc.h(anc)
        # Gate 0
        dq, pp = stab[0]
        if pp == 'X': qc.cx(anc, dq)
        elif pp == 'Z': qc.cz(anc, dq)
        # Flag CNOT 1
        qc.cx(anc, flg)
        # Gates 1, 2 (inside bracket)
        for gi in [1, 2]:
            dq, pp = stab[gi]
            if pp == 'X': qc.cx(anc, dq)
            elif pp == 'Z': qc.cz(anc, dq)
        # Flag CNOT 2
        qc.cx(anc, flg)
        # Gate 3
        dq, pp = stab[3]
        if pp == 'X': qc.cx(anc, dq)
        elif pp == 'Z': qc.cz(anc, dq)
        qc.h(anc)
    return qc

bare_synd = build_bare_syndrome()
flag_synd = build_flag_syndrome()

def build_uncoded(logical, basis):
    qc = QuantumCircuit(1, 1)
    if logical == '1': qc.x(0)
    if basis == 'X': qc.h(0)
    elif basis == 'Y': qc.h(0); qc.s(0)
    # Measure in basis
    if basis == 'X': qc.h(0)
    elif basis == 'Y': qc.sdg(0); qc.h(0)
    qc.measure(0, 0)
    return qc

def build_passive(logical, basis):
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
    return qc

def build_bare_active(logical, basis):
    """Bare syndrome: 9 qubits, 9 classical bits."""
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
    return qc

def build_flag_active(logical, basis):
    """Flag-FT syndrome: 13 qubits, 13 classical bits."""
    qc = QuantumCircuit(13, 13)
    if logical == '1': qc.x(0)
    if basis == 'X': qc.h(0)
    elif basis == 'Y': qc.h(0); qc.s(0)
    qc.compose(enc_513, qubits=range(5), inplace=True)
    qc.compose(flag_synd, qubits=range(13), inplace=True)
    # Measure syndrome
    qc.measure([5,6,7,8], [5,6,7,8])
    # Measure flags
    qc.measure([9,10,11,12], [9,10,11,12])
    # Measure data in logical basis
    if basis == 'X':
        for q in range(5): qc.h(q)
    elif basis == 'Y':
        for q in range(5): qc.sdg(q); qc.h(q)
    qc.measure(range(5), range(5))
    return qc

# ============================================================
# Noise model (IBM hardware-calibrated)
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
    p = max(1e-12, min(1-1e-12, p))
    return 1.0 + p * np.log2(p) + (1-p) * np.log2(1-p)

# ============================================================
# Post-processing decoders
# ============================================================
def decode_bare(counts, basis, logical_target):
    """Decode bare syndrome: read syndrome, apply standard correction, decode parity."""
    n_total = 0; n_error = 0
    syn_hist = {}
    for bitstring, count in counts.items():
        bits = bitstring[::-1]
        data_bits = [int(bits[i]) for i in range(5)]
        syn_bits = tuple(int(bits[i]) for i in range(5, 9))
        syn_hist[str(syn_bits)] = syn_hist.get(str(syn_bits), 0) + count

        correction = syndrome_table.get(syn_bits)
        if correction is not None:
            cq, cp = correction
            if basis == 'Z' and cp in ('X', 'Y'): data_bits[cq] ^= 1
            elif basis == 'X' and cp in ('Z', 'Y'): data_bits[cq] ^= 1
            elif basis == 'Y' and cp in ('X', 'Z'): data_bits[cq] ^= 1

        parity = sum(data_bits) % 2
        n_total += count
        if parity != logical_target: n_error += count
    return n_error / n_total, syn_hist

def decode_flag(counts, basis, logical_target):
    """Decode flag-FT: read syndrome+flags, apply flag-aware correction, decode parity."""
    n_total = 0; n_error = 0
    syn_hist = {}
    flag_stats = {'flagged': 0, 'flag_match': 0, 'total': 0}

    for bitstring, count in counts.items():
        bits = bitstring[::-1]
        data_bits = [int(bits[i]) for i in range(5)]
        syn_bits = tuple(int(bits[i]) for i in range(5, 9))
        flag_bits = tuple(int(bits[i]) for i in range(9, 13))

        syn_hist[str(syn_bits)] = syn_hist.get(str(syn_bits), 0) + count
        flag_stats['total'] += count
        any_flag = any(f == 1 for f in flag_bits)
        if any_flag:
            flag_stats['flagged'] += count

        # Flag-aware correction
        correction = None
        used_flag = False

        # Check if any flag matches its weight-2 syndrome
        for f_idx in range(4):
            if flag_bits[f_idx] == 1:
                if syn_bits == flag_weight2[f_idx]['syndrome']:
                    # Apply weight-2 correction
                    correction = flag_weight2[f_idx]['correction']
                    used_flag = True
                    flag_stats['flag_match'] += count
                    break

        if not used_flag:
            # Standard correction
            std_corr = syndrome_table.get(syn_bits)
            if std_corr is not None:
                correction = [std_corr]  # single (qubit, pauli)

        # Apply correction to data bits
        if correction is not None:
            for cq, cp in correction:
                if basis == 'Z' and cp in ('X', 'Y'): data_bits[cq] ^= 1
                elif basis == 'X' and cp in ('Z', 'Y'): data_bits[cq] ^= 1
                elif basis == 'Y' and cp in ('X', 'Z'): data_bits[cq] ^= 1

        parity = sum(data_bits) % 2
        n_total += count
        if parity != logical_target: n_error += count
    return n_error / n_total, syn_hist, flag_stats

# ============================================================
# Run simulations
# ============================================================
print(f"Running noisy simulations: p1q={p1q}, p2q={p2q}, p_ro={p_ro}")
print(f"Shots: {shots}\n")

results = {}

# 1. UNCODED
print("--- Uncoded ---")
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
    err0 = sum(v for k,v in c0.items() if sum(int(b) for b in k) % 2 == 1) / shots
    err1 = sum(v for k,v in c1.items() if sum(int(b) for b in k) % 2 == 0) / shots
    avg_err = (err0 + err1) / 2
    h = holevo_from_err(avg_err)
    results['passive'][basis] = {'err': round(float(avg_err), 5), 'holevo': round(float(h), 4)}
    print(f"  {basis}: err={avg_err:.4f} Holevo={h:.4f}")

# 3. BARE SYNDROME (Sprint 026 approach)
print("\n--- Bare syndrome [[5,1,3]] ---")
results['bare'] = {}
for basis in ['Z', 'X', 'Y']:
    qc0 = pm.run(build_bare_active('0', basis))
    qc1 = pm.run(build_bare_active('1', basis))
    c0 = sim.run(qc0, shots=shots).result().get_counts()
    c1 = sim.run(qc1, shots=shots).result().get_counts()
    err0, sh0 = decode_bare(c0, basis, 0)
    err1, sh1 = decode_bare(c1, basis, 1)
    avg_err = (err0 + err1) / 2
    h = holevo_from_err(avg_err)
    results['bare'][basis] = {'err': round(float(avg_err), 5), 'holevo': round(float(h), 4)}
    print(f"  {basis}: err={avg_err:.4f} Holevo={h:.4f}")

# 4. FLAG-FT SYNDROME
print("\n--- Flag-FT syndrome [[5,1,3]] ---")
results['flag_ft'] = {}
total_flag_stats = {'flagged': 0, 'flag_match': 0, 'total': 0}
for basis in ['Z', 'X', 'Y']:
    qc0 = pm.run(build_flag_active('0', basis))
    qc1 = pm.run(build_flag_active('1', basis))
    c0 = sim.run(qc0, shots=shots).result().get_counts()
    c1 = sim.run(qc1, shots=shots).result().get_counts()
    err0, sh0, fs0 = decode_flag(c0, basis, 0)
    err1, sh1, fs1 = decode_flag(c1, basis, 1)
    avg_err = (err0 + err1) / 2
    h = holevo_from_err(avg_err)
    # Merge flag stats
    for k in total_flag_stats:
        total_flag_stats[k] += fs0[k] + fs1[k]
    results['flag_ft'][basis] = {
        'err': round(float(avg_err), 5),
        'holevo': round(float(h), 4),
        'flag_rate': round((fs0['flagged']+fs1['flagged'])/(fs0['total']+fs1['total']), 4),
        'flag_match_rate': round((fs0['flag_match']+fs1['flag_match'])/(fs0['total']+fs1['total']), 5),
    }
    print(f"  {basis}: err={avg_err:.4f} Holevo={h:.4f} flag_rate={results['flag_ft'][basis]['flag_rate']:.3f}")

# ============================================================
# Summary
# ============================================================
print(f"\n{'='*70}")
print(f"SUMMARY (p2q={p2q})")
print(f"{'='*70}")

for strategy in ['uncoded', 'passive', 'bare', 'flag_ft']:
    r = results[strategy]
    h_vals = [r[b]['holevo'] for b in ['Z', 'X', 'Y']]
    avg_h = sum(h_vals) / 3
    asym = max(h_vals) - min(h_vals)
    min_h = min(h_vals)
    results[strategy]['avg_holevo'] = round(float(avg_h), 4)
    results[strategy]['asymmetry'] = round(float(asym), 4)
    results[strategy]['min_holevo'] = round(float(min_h), 4)
    print(f"\n  {strategy:>10s}: Z={h_vals[0]:.4f}  X={h_vals[1]:.4f}  Y={h_vals[2]:.4f}")
    print(f"             avg={avg_h:.4f}  min={min_h:.4f}  asym={asym:.4f}")

# Deltas
print(f"\n  Bare vs passive:    {results['bare']['avg_holevo'] - results['passive']['avg_holevo']:+.4f}")
print(f"  Flag-FT vs passive: {results['flag_ft']['avg_holevo'] - results['passive']['avg_holevo']:+.4f}")
print(f"  Flag-FT vs bare:    {results['flag_ft']['avg_holevo'] - results['bare']['avg_holevo']:+.4f}")

# Flag statistics
frac_flagged = total_flag_stats['flagged'] / max(1, total_flag_stats['total'])
frac_match = total_flag_stats['flag_match'] / max(1, total_flag_stats['total'])
print(f"\n  Flag statistics:")
print(f"    Any flag triggered: {frac_flagged:.3%} of shots")
print(f"    Flag matched w2 syndrome: {frac_match:.3%} of shots")

elapsed = time.time() - start
print(f"\nTotal time: {elapsed:.1f}s")

output = {
    'sprint': '027b',
    'description': 'Flag-FT vs bare vs passive under hardware noise',
    'noise_model': {'p1q': p1q, 'p2q': p2q, 'p_readout': p_ro},
    'shots': shots,
    'results': results,
    'flag_statistics': {
        'frac_flagged': round(frac_flagged, 4),
        'frac_match': round(frac_match, 5),
        'raw': total_flag_stats,
    },
    'elapsed_seconds': round(elapsed, 1),
}

with open('results/sprint_027b_flag_vs_bare.json', 'w') as f:
    json.dump(output, f, indent=2)
print("Saved to results/sprint_027b_flag_vs_bare.json")
