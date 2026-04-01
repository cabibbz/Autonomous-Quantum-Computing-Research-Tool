"""
Sprint 025a: Simulator baseline for hardware comparison
Efficient Clifford encoding circuits for 3-qubit and [[5,1,3]] codes.
Gate-level noise simulation matching IBM hardware characteristics.
Measures basis-dependent Holevo information from bitstring decoding.

Key prediction: [[5,1,3]] should be basis-isotropic, 3-qubit highly anisotropic.
"""

import numpy as np
import json, time
from qiskit.quantum_info import Clifford, Statevector
from qiskit.circuit import QuantumCircuit
from qiskit_aer import AerSimulator
from qiskit_aer.noise import NoiseModel, depolarizing_error, ReadoutError
from qiskit.transpiler import generate_preset_pass_manager

start = time.time()

# ============================================================
# Build [[5,1,3]] encoding circuit from Clifford tableau
# ============================================================
# Verified destabilizers: D1=XIXXX, D2=IXXII, D3=XIIIX, D4=IXXXX (all pure-X)
# Stabilizers: g1=XZZXI, g2=IXZZX, g3=XIXZZ, g4=ZXIXZ
# Logical: X_L=XXXXX, Z_L=ZZZZZ

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
cx_count_513 = enc_513.count_ops().get('cx', 0) + enc_513.count_ops().get('swap', 0) * 3
print(f"[[5,1,3]] encoding: {enc_513.count_ops()}, effective CX~{cx_count_513}")

# Quick verification
I2 = np.eye(2, dtype=complex)
Xm = np.array([[0,1],[1,0]], dtype=complex)
Zm = np.array([[1,0],[0,-1]], dtype=complex)

def multi_kron(ops, nq):
    r = np.eye(1, dtype=complex)
    for q in range(nq):
        r = np.kron(r, ops.get(q, I2))
    return r

dim = 32
stab_ops = [{0:Xm,1:Zm,2:Zm,3:Xm},{1:Xm,2:Zm,3:Zm,4:Xm},{0:Xm,2:Xm,3:Zm,4:Zm},{0:Zm,1:Xm,3:Xm,4:Zm}]
proj = np.eye(dim, dtype=complex)
for s in stab_ops:
    S = multi_kron(s, 5)
    proj = proj @ (np.eye(dim) + S) / 2

s0 = np.zeros(dim, dtype=complex); s0[0] = 1
s0 = proj @ s0; s0 /= np.linalg.norm(s0)
X_all = multi_kron({q: Xm for q in range(5)}, 5)
s1 = X_all @ s0; s1 /= np.linalg.norm(s1)

# Verify encoding
for label, prep in [('|0>', None), ('|1>', 'x'), ('|+>', 'h')]:
    qc = QuantumCircuit(5)
    if prep == 'x': qc.x(0)
    elif prep == 'h': qc.h(0)
    qc.compose(enc_513, inplace=True)
    sv = Statevector.from_instruction(qc).data
    if prep is None: target = s0
    elif prep == 'x': target = s1
    else: target = (s0 + s1) / np.sqrt(2)
    fid = abs(np.dot(target.conj(), sv))**2
    print(f"  {label} fidelity: {fid:.6f}")

# ============================================================
# Circuit builders
# ============================================================
def build_circuit(code, logical, basis):
    """
    Build prepare-and-measure circuit.
    code: 'uncoded', '3q', '513'
    logical: '0' or '1'
    basis: 'Z', 'X', 'Y'
    """
    if code == 'uncoded':
        qc = QuantumCircuit(1, 1)
        # Prepare
        if basis == 'Z':
            if logical == '1': qc.x(0)
        elif basis == 'X':
            if logical == '1': qc.x(0)
            qc.h(0)
        elif basis == 'Y':
            if logical == '1': qc.x(0)
            qc.h(0); qc.s(0)
        # Measure in same basis
        if basis == 'X': qc.h(0)
        elif basis == 'Y': qc.sdg(0); qc.h(0)
        qc.measure(0, 0)
        return qc

    elif code == '3q':
        qc = QuantumCircuit(3, 3)
        # Prepare logical state on qubit 0
        if basis == 'Z':
            if logical == '1': qc.x(0)
        elif basis == 'X':
            if logical == '1': qc.x(0)
            qc.h(0)
        elif basis == 'Y':
            if logical == '1': qc.x(0)
            qc.h(0); qc.s(0)
        # Encode: CNOT ladder
        qc.cx(0, 1); qc.cx(0, 2)
        # Measure in matching basis
        if basis == 'X':
            qc.h(0); qc.h(1); qc.h(2)
        elif basis == 'Y':
            qc.sdg(0); qc.h(0); qc.sdg(1); qc.h(1); qc.sdg(2); qc.h(2)
        qc.measure([0, 1, 2], [0, 1, 2])
        return qc

    elif code == '513':
        qc = QuantumCircuit(5, 5)
        # Prepare logical state on qubit 0
        if basis == 'Z':
            if logical == '1': qc.x(0)
        elif basis == 'X':
            if logical == '1': qc.x(0)
            qc.h(0)
        elif basis == 'Y':
            if logical == '1': qc.x(0)
            qc.h(0); qc.s(0)
        # Encode
        qc.compose(enc_513, inplace=True)
        # Measure in matching basis
        # Z-basis: logical Z = ZZZZZ -> parity of all bits
        # X-basis: apply H to all -> XXXXX becomes ZZZZZ -> parity
        # Y-basis: apply S†H to all -> YYYYY becomes ZZZZZ -> parity
        if basis == 'X':
            for q in range(5): qc.h(q)
        elif basis == 'Y':
            for q in range(5): qc.sdg(q); qc.h(q)
        qc.measure(range(5), range(5))
        return qc

# ============================================================
# Noise model and simulation
# ============================================================
print("\nRunning gate-level noise simulations...")

# IBM Heron/Eagle typical error rates
p1q = 0.0005    # single-qubit gate error
p2q = 0.008     # two-qubit gate error
p_ro = 0.015    # readout error

noise = NoiseModel()
noise.add_all_qubit_quantum_error(depolarizing_error(p1q, 1), ['h', 'x', 's', 'sdg', 'z', 'rz', 'sx', 'id', 'y'])
noise.add_all_qubit_quantum_error(depolarizing_error(p2q, 2), ['cx', 'cz', 'ecr', 'swap'])
noise.add_all_qubit_readout_error(ReadoutError([[1-p_ro, p_ro], [p_ro, 1-p_ro]]))

sim = AerSimulator(noise_model=noise)
pm = generate_preset_pass_manager(optimization_level=2, basis_gates=['cx', 'id', 'rz', 'sx', 'x'])
shots = 50000

def decode_majority(bs):
    return 1 if sum(int(b) for b in bs) >= 2 else 0

def decode_parity(bs):
    return sum(int(b) for b in bs) % 2

def holevo_from_err(p):
    p = max(1e-12, min(1 - 1e-12, p))
    return 1.0 + p * np.log2(p) + (1 - p) * np.log2(1 - p)

results = {}
codes = [('uncoded', 'uncoded'), ('3-qubit', '3q'), ('[[5,1,3]]', '513')]

for label, code in codes:
    results[label] = {}
    print(f"\n--- {label} ---")

    for basis in ['Z', 'X', 'Y']:
        qc0 = build_circuit(code, '0', basis)
        qc1 = build_circuit(code, '1', basis)

        # Transpile
        qc0t = pm.run(qc0)
        qc1t = pm.run(qc1)
        d0 = qc0t.depth(); cx0 = qc0t.count_ops().get('cx', 0)
        d1 = qc1t.depth(); cx1 = qc1t.count_ops().get('cx', 0)

        # Simulate
        counts0 = sim.run(qc0t, shots=shots).result().get_counts()
        counts1 = sim.run(qc1t, shots=shots).result().get_counts()

        # Decode: 3-qubit uses majority for Z (corrects bit flips), parity for X/Y
        # Note: Y_L = -YYY for 3-qubit (n=3), so Y-basis parity must be inverted
        # This is because Y_L = i*X_L*Z_L = i*(XXX)*(ZZZ) = i*(-iY)^3 = -YYY
        # For [[5,1,3]] (n=5): Y_L = YYYYY (no sign), parity works directly
        if code == 'uncoded':
            err0 = sum(v for k, v in counts0.items() if int(k) == 1) / shots
            err1 = sum(v for k, v in counts1.items() if int(k) == 0) / shots
        elif code == '3q':
            if basis == 'Z':
                err0 = sum(v for k, v in counts0.items() if decode_majority(k) == 1) / shots
                err1 = sum(v for k, v in counts1.items() if decode_majority(k) == 0) / shots
            elif basis == 'X':
                err0 = sum(v for k, v in counts0.items() if decode_parity(k) == 1) / shots
                err1 = sum(v for k, v in counts1.items() if decode_parity(k) == 0) / shots
            else:  # Y: inverted parity due to Y_L = -YYY
                err0 = sum(v for k, v in counts0.items() if decode_parity(k) == 0) / shots
                err1 = sum(v for k, v in counts1.items() if decode_parity(k) == 1) / shots
        else:  # 513: parity works for all bases (X_L=XXXXX, Z_L=ZZZZZ, Y_L=YYYYY)
            err0 = sum(v for k, v in counts0.items() if decode_parity(k) == 1) / shots
            err1 = sum(v for k, v in counts1.items() if decode_parity(k) == 0) / shots

        avg_err = (err0 + err1) / 2
        holevo = holevo_from_err(avg_err)

        results[label][basis] = {
            'err_0': round(float(err0), 5),
            'err_1': round(float(err1), 5),
            'avg_err': round(float(avg_err), 5),
            'holevo': round(float(holevo), 4),
            'depth': [d0, d1],
            'cx': [cx0, cx1],
        }
        print(f"  {basis}: err={avg_err:.4f} Holevo={holevo:.4f} depth={d0}/{d1} CX={cx0}/{cx1}")

    h = [results[label][b]['holevo'] for b in ['Z', 'X', 'Y']]
    asym = max(h) - min(h)
    avg = sum(h) / 3
    results[label]['asymmetry'] = round(float(asym), 4)
    results[label]['avg_holevo'] = round(float(avg), 4)
    print(f"  => Asymmetry: {asym:.4f}  Avg Holevo: {avg:.4f}")

# ============================================================
# Also run ideal (no noise) for reference
# ============================================================
print("\n--- Ideal (no noise) ---")
sim_ideal = AerSimulator()
for label, code in codes:
    results[label]['ideal'] = {}
    for basis in ['Z', 'X', 'Y']:
        qc0 = build_circuit(code, '0', basis)
        qc1 = build_circuit(code, '1', basis)
        counts0 = sim_ideal.run(qc0, shots=shots).result().get_counts()
        counts1 = sim_ideal.run(qc1, shots=shots).result().get_counts()

        if code == 'uncoded':
            err0 = sum(v for k, v in counts0.items() if int(k) == 1) / shots
            err1 = sum(v for k, v in counts1.items() if int(k) == 0) / shots
        elif code == '3q':
            if basis == 'Z':
                err0 = sum(v for k, v in counts0.items() if decode_majority(k) == 1) / shots
                err1 = sum(v for k, v in counts1.items() if decode_majority(k) == 0) / shots
            elif basis == 'X':
                err0 = sum(v for k, v in counts0.items() if decode_parity(k) == 1) / shots
                err1 = sum(v for k, v in counts1.items() if decode_parity(k) == 0) / shots
            else:  # Y: inverted parity
                err0 = sum(v for k, v in counts0.items() if decode_parity(k) == 0) / shots
                err1 = sum(v for k, v in counts1.items() if decode_parity(k) == 1) / shots
        else:
            err0 = sum(v for k, v in counts0.items() if decode_parity(k) == 1) / shots
            err1 = sum(v for k, v in counts1.items() if decode_parity(k) == 0) / shots
        results[label]['ideal'][basis] = round(float((err0+err1)/2), 5)

    print(f"  {label} ideal errors: Z={results[label]['ideal']['Z']}, X={results[label]['ideal']['X']}, Y={results[label]['ideal']['Y']}")

# ============================================================
# Summary
# ============================================================
print(f"\n{'='*70}")
print(f"SUMMARY: Gate-Level Noise Baseline (p1q={p1q}, p2q={p2q}, ro={p_ro})")
print(f"{'='*70}")
print(f"{'Code':>12s} | {'Z-Hol':>7s} | {'X-Hol':>7s} | {'Y-Hol':>7s} | {'Avg':>7s} | {'Asym':>6s} | CX(Z)")
print(f"{'-'*65}")
for label, _ in codes:
    r = results[label]
    cx_z = r['Z']['cx'][0]
    print(f"{label:>12s} | {r['Z']['holevo']:7.4f} | {r['X']['holevo']:7.4f} | {r['Y']['holevo']:7.4f} | {r['avg_holevo']:7.4f} | {r['asymmetry']:6.4f} | {cx_z}")

print(f"\nKey predictions for hardware:")
print(f"  1. 3-qubit asymmetry >> [[5,1,3]] asymmetry (anisotropy prediction)")
print(f"  2. [[5,1,3]] avg Holevo vs 3-qubit avg Holevo (encoding overhead trade-off)")
print(f"  3. Uncoded should beat both codes (gate noise dominates at this scale)")

elapsed = time.time() - start
print(f"\nTotal time: {elapsed:.1f}s")

output = {
    'sprint': '025a',
    'description': 'Gate-level noise baseline for hardware comparison',
    'noise_model': {'p1q': p1q, 'p2q': p2q, 'p_readout': p_ro},
    'shots': shots,
    'encoding_circuit_513': {
        'cx_count': enc_513.count_ops().get('cx', 0),
        'depth': enc_513.depth(),
        'gates': dict(enc_513.count_ops()),
    },
    'results': results,
    'elapsed_seconds': round(elapsed, 1),
}

with open('results/sprint_025a_hardware_baseline.json', 'w') as f:
    json.dump(output, f, indent=2)
print("Saved to results/sprint_025a_hardware_baseline.json")
