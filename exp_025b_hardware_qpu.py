"""
Sprint 025b: Real IBM QPU test — basis-dependent QEC performance
Same circuits as 025a, submitted to real IBM quantum hardware.

QPU justification:
1. Testing whether [[5,1,3]]'s basis isotropy (asymmetry ~0.01 in simulation)
   survives on real hardware with correlated noise, crosstalk, and qubit-specific T1/T2.
2. Surprising outcome: if hardware asymmetry > 0.05 for [[5,1,3]], it means
   correlated/structured noise breaks the code's symmetry — unmapped territory.
3. If hardware matches simulation, we validate 24 sprints of noise models.
   If not, the deviation quantifies what our models miss.

QPU budget: 0s used of 600s. Expected usage: ~20-40s (18 circuits × 4096 shots).
"""

import numpy as np
import json, time
from qiskit.quantum_info import Clifford
from qiskit.circuit import QuantumCircuit
from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2 as Sampler
from qiskit.transpiler import generate_preset_pass_manager

start = time.time()

# ============================================================
# Build encoding circuits (same as 025a)
# ============================================================
n = 5
destab_x = np.array([[1,1,1,1,1],[1,0,1,1,1],[0,1,1,0,0],[1,0,0,0,1],[0,1,1,1,1]], dtype=bool)
destab_z = np.zeros((5,5), dtype=bool)
stab_x = np.array([[0,0,0,0,0],[1,0,0,1,0],[0,1,0,0,1],[1,0,1,0,0],[0,1,0,1,0]], dtype=bool)
stab_z = np.array([[1,1,1,1,1],[0,1,1,0,0],[0,0,1,1,0],[0,0,0,1,1],[1,0,0,0,1]], dtype=bool)

tableau = np.zeros((10, 11), dtype=bool)
tableau[:5, :5] = destab_x; tableau[:5, 5:10] = destab_z
tableau[5:, :5] = stab_x; tableau[5:, 5:10] = stab_z

cliff = Clifford(tableau)
enc_513 = cliff.to_circuit()
print(f"[[5,1,3]] encoding circuit: {enc_513.count_ops()}")

# ============================================================
# Circuit builders (identical to 025a)
# ============================================================
def build_circuit(code, logical, basis):
    if code == 'uncoded':
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

    elif code == '3q':
        qc = QuantumCircuit(3, 3)
        if basis == 'Z':
            if logical == '1': qc.x(0)
        elif basis == 'X':
            if logical == '1': qc.x(0)
            qc.h(0)
        elif basis == 'Y':
            if logical == '1': qc.x(0)
            qc.h(0); qc.s(0)
        qc.cx(0, 1); qc.cx(0, 2)
        if basis == 'X':
            qc.h(0); qc.h(1); qc.h(2)
        elif basis == 'Y':
            qc.sdg(0); qc.h(0); qc.sdg(1); qc.h(1); qc.sdg(2); qc.h(2)
        qc.measure([0, 1, 2], [0, 1, 2])
        return qc

    elif code == '513':
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

# ============================================================
# Connect to IBM and build all circuits
# ============================================================
print("\nConnecting to IBM Quantum...")
service = QiskitRuntimeService()
backend = service.least_busy(operational=True, simulator=False)
print(f"Backend: {backend.name}")
print(f"Qubits: {backend.num_qubits}")

# Get backend properties for context
try:
    props = backend.properties()
    if props:
        t1s = [props.qubit_property(i).get('T1', [None])[0] for i in range(min(5, backend.num_qubits))]
        t2s = [props.qubit_property(i).get('T2', [None])[0] for i in range(min(5, backend.num_qubits))]
        print(f"T1 (first 5 qubits): {t1s}")
        print(f"T2 (first 5 qubits): {t2s}")
except Exception as e:
    print(f"Could not get properties: {e}")

# Build and transpile all circuits
shots = 4096
pm = generate_preset_pass_manager(backend=backend, optimization_level=2)

circuits = {}
transpiled = {}
circuit_info = {}

for code_label, code_key in [('uncoded', 'uncoded'), ('3-qubit', '3q'), ('[[5,1,3]]', '513')]:
    for basis in ['Z', 'X', 'Y']:
        for logical in ['0', '1']:
            key = f"{code_label}_{basis}_{logical}"
            qc = build_circuit(code_key, logical, basis)
            qc.name = key
            circuits[key] = qc
            tc = pm.run(qc)
            transpiled[key] = tc
            circuit_info[key] = {
                'depth': tc.depth(),
                'cx': tc.count_ops().get('cx', 0) + tc.count_ops().get('ecr', 0),
                'gates': dict(tc.count_ops()),
            }
            print(f"  {key}: depth={tc.depth()}, 2Q gates={circuit_info[key]['cx']}")

# ============================================================
# Submit to QPU
# ============================================================
print(f"\nSubmitting {len(transpiled)} circuits to {backend.name}...")
print("QPU budget check: 0s used, ~30-40s expected")

# Submit all circuits in one batch
all_circuits = list(transpiled.values())
all_keys = list(transpiled.keys())

sampler = Sampler(mode=backend)
job = sampler.run(all_circuits, shots=shots)
print(f"Job ID: {job.job_id()}")
print("Waiting for results...")

result = job.result()

# Extract metrics
try:
    metrics = job.metrics()
    qpu_seconds = metrics.get('usage', {}).get('quantum_seconds', 'unknown')
    print(f"QPU time used: {qpu_seconds}s")
except Exception as e:
    qpu_seconds = 'unknown'
    print(f"Could not get metrics: {e}")

# ============================================================
# Decode results
# ============================================================
print("\nDecoding results...")

def decode_majority(bs):
    return 1 if sum(int(b) for b in bs) >= 2 else 0

def decode_parity(bs):
    return sum(int(b) for b in bs) % 2

def holevo_from_err(p):
    p = max(1e-12, min(1 - 1e-12, p))
    return 1.0 + p * np.log2(p) + (1 - p) * np.log2(1 - p)

results = {}
raw_counts = {}

for i, key in enumerate(all_keys):
    counts = result[i].data.meas.get_counts()
    raw_counts[key] = counts

for code_label, code_key in [('uncoded', 'uncoded'), ('3-qubit', '3q'), ('[[5,1,3]]', '513')]:
    results[code_label] = {}
    print(f"\n=== {code_label} ===")

    for basis in ['Z', 'X', 'Y']:
        key0 = f"{code_label}_{basis}_0"
        key1 = f"{code_label}_{basis}_1"
        counts0 = raw_counts[key0]
        counts1 = raw_counts[key1]

        if code_key == 'uncoded':
            err0 = sum(v for k, v in counts0.items() if int(k) == 1) / shots
            err1 = sum(v for k, v in counts1.items() if int(k) == 0) / shots
        elif code_key == '3q':
            if basis == 'Z':
                err0 = sum(v for k, v in counts0.items() if decode_majority(k) == 1) / shots
                err1 = sum(v for k, v in counts1.items() if decode_majority(k) == 0) / shots
            elif basis == 'X':
                err0 = sum(v for k, v in counts0.items() if decode_parity(k) == 1) / shots
                err1 = sum(v for k, v in counts1.items() if decode_parity(k) == 0) / shots
            else:  # Y: inverted parity (Y_L = -YYY for n=3)
                err0 = sum(v for k, v in counts0.items() if decode_parity(k) == 0) / shots
                err1 = sum(v for k, v in counts1.items() if decode_parity(k) == 1) / shots
        else:  # 513
            err0 = sum(v for k, v in counts0.items() if decode_parity(k) == 1) / shots
            err1 = sum(v for k, v in counts1.items() if decode_parity(k) == 0) / shots

        avg_err = (err0 + err1) / 2
        holevo = holevo_from_err(avg_err)

        results[code_label][basis] = {
            'err_0': round(float(err0), 5),
            'err_1': round(float(err1), 5),
            'avg_err': round(float(avg_err), 5),
            'holevo': round(float(holevo), 4),
            'depth': circuit_info[key0]['depth'],
            'cx': circuit_info[key0]['cx'],
        }
        print(f"  {basis}: err={avg_err:.4f} Holevo={holevo:.4f} depth={circuit_info[key0]['depth']} CX={circuit_info[key0]['cx']}")

    h = [results[code_label][b]['holevo'] for b in ['Z', 'X', 'Y']]
    asym = max(h) - min(h)
    avg = sum(h) / 3
    results[code_label]['asymmetry'] = round(float(asym), 4)
    results[code_label]['avg_holevo'] = round(float(avg), 4)
    print(f"  => Asymmetry: {asym:.4f}  Avg Holevo: {avg:.4f}")

# ============================================================
# Summary
# ============================================================
print(f"\n{'='*70}")
print(f"HARDWARE RESULTS: {backend.name}")
print(f"QPU time: {qpu_seconds}s, Shots: {shots}")
print(f"{'='*70}")
print(f"{'Code':>12s} | {'Z-Hol':>7s} | {'X-Hol':>7s} | {'Y-Hol':>7s} | {'Avg':>7s} | {'Asym':>6s} | CX")
print(f"{'-'*65}")
for cn in ['uncoded', '3-qubit', '[[5,1,3]]']:
    r = results[cn]
    cx = r['Z']['cx']
    print(f"{cn:>12s} | {r['Z']['holevo']:7.4f} | {r['X']['holevo']:7.4f} | {r['Y']['holevo']:7.4f} | {r['avg_holevo']:7.4f} | {r['asymmetry']:6.4f} | {cx}")

# Comparison with simulator predictions (from 025a)
print("\n--- Comparison: Simulator vs Hardware ---")
sim_results = {
    'uncoded': {'Z': 0.887, 'X': 0.887, 'Y': 0.885, 'asym': 0.002, 'avg': 0.886},
    '3-qubit': {'Z': 0.939, 'X': 0.694, 'Y': 0.701, 'asym': 0.245, 'avg': 0.778},
    '[[5,1,3]]': {'Z': 0.506, 'X': 0.516, 'Y': 0.508, 'asym': 0.010, 'avg': 0.510},
}
for cn in ['uncoded', '3-qubit', '[[5,1,3]]']:
    sr = sim_results[cn]
    hr = results[cn]
    print(f"\n{cn}:")
    print(f"  Sim avg Holevo: {sr['avg']:.3f}  Hardware: {hr['avg_holevo']:.3f}  Gap: {sr['avg'] - hr['avg_holevo']:+.3f}")
    print(f"  Sim asymmetry:  {sr['asym']:.3f}  Hardware: {hr['asymmetry']:.3f}  Gap: {sr['asym'] - hr['asymmetry']:+.3f}")

elapsed = time.time() - start
print(f"\nTotal time: {elapsed:.1f}s")

# Save results
output = {
    'sprint': '025b',
    'description': 'Real IBM QPU basis-dependent QEC performance',
    'backend': backend.name,
    'shots': shots,
    'qpu_seconds': qpu_seconds,
    'job_id': job.job_id(),
    'circuit_info': circuit_info,
    'results': results,
    'raw_counts': {k: dict(v) for k, v in raw_counts.items()},
    'elapsed_seconds': round(elapsed, 1),
}

with open('results/sprint_025b_hardware_qpu.json', 'w') as f:
    json.dump(output, f, indent=2)
print("Saved to results/sprint_025b_hardware_qpu.json")
