"""
Sprint 012c: Scrambling Efficiency — Structured vs Random Circuits

How do GHZ, W, and Cluster preparation circuits compare to random circuits
as scramblers? We measure:
1. Entanglement per gate (efficiency)
2. OTOC decay rate
3. Final scrambling quality (how close to Page entropy)

For each structured circuit, we also compare what happens when we continue
evolving with random gates after preparation — do structured states
scramble faster or slower as a starting point?
"""

import numpy as np
import json
from qiskit.circuit import QuantumCircuit
from qiskit.quantum_info import Statevector, partial_trace, entropy, random_unitary, Operator

n = 6

def von_neumann_entropy_sv(sv, qubits_to_keep, n_qubits):
    qubits_to_trace = [q for q in range(n_qubits) if q not in qubits_to_keep]
    if not qubits_to_trace:
        return 0.0
    rdm = partial_trace(sv, qubits_to_trace)
    return float(entropy(rdm, base=2))

def mutual_information_sv(sv, i, j, n_qubits):
    si = von_neumann_entropy_sv(sv, [i], n_qubits)
    sj = von_neumann_entropy_sv(sv, [j], n_qubits)
    sij = von_neumann_entropy_sv(sv, [i, j], n_qubits)
    return si + sj - sij

def compute_otoc_from_sv(psi0_vec, U_mat, n_qubits, w_qubit=0, v_qubit=None):
    """Compute average OTOC across all probe qubits."""
    Udag = U_mat.conj().T
    dim = 2**n_qubits

    # W = Z on w_qubit
    W = np.ones(dim)
    for idx in range(dim):
        if (idx >> (n_qubits - 1 - w_qubit)) & 1:
            W[idx] = -1.0

    f_vals = []
    probe_qubits = [v_qubit] if v_qubit is not None else range(n_qubits)
    for vq in probe_qubits:
        if vq == w_qubit:
            continue
        # V = Z on v_qubit
        V = np.ones(dim)
        for idx in range(dim):
            if (idx >> (n_qubits - 1 - vq)) & 1:
                V[idx] = -1.0

        # |e⟩ = V U† W U V |ψ0⟩
        e = V * psi0_vec
        e = U_mat @ e
        e = W * e
        e = Udag @ e
        e = V * e

        # |h⟩ = U† W U |ψ0⟩
        h = U_mat @ psi0_vec
        h = W * h
        h = Udag @ h

        F = np.abs(np.dot(h.conj(), e))**2
        f_vals.append(F)

    return float(np.mean(f_vals))

def make_ghz_circuit(n_qubits):
    qc = QuantumCircuit(n_qubits)
    qc.h(0)
    for i in range(1, n_qubits):
        qc.cx(0, i)
    return qc

def make_w_circuit(n_qubits):
    """Approximate W state preparation."""
    qc = QuantumCircuit(n_qubits)
    qc.x(0)
    for i in range(n_qubits - 1):
        angle = np.arccos(np.sqrt(1.0 / (n_qubits - i)))
        qc.ry(2 * angle, i)
        qc.cx(i, i + 1)
        qc.x(i)
    return qc

def make_cluster_1d_circuit(n_qubits):
    qc = QuantumCircuit(n_qubits)
    for i in range(n_qubits):
        qc.h(i)
    for i in range(n_qubits - 1):
        qc.cz(i, i + 1)
    return qc

def add_random_layers(qc, n_qubits, n_layers, seed):
    """Add random 1D brick-wall layers."""
    rng = np.random.default_rng(seed)
    for l in range(n_layers):
        offset = l % 2
        for i in range(offset, n_qubits - 1, 2):
            u = random_unitary(4, seed=int(rng.integers(0, 2**31)))
            qc.unitary(u, [i, i + 1])
    return qc

# Part 1: Scrambling quality of preparation circuits alone
print("=== Part 1: Scrambling quality of preparation circuits ===")
max_ent = n / 2  # max half-cut entropy

states = {
    'GHZ': make_ghz_circuit(n),
    'W': make_w_circuit(n),
    'Cluster_1D': make_cluster_1d_circuit(n),
}

state_metrics = {}
for name, qc in states.items():
    sv = Statevector.from_instruction(qc)
    half_ent = von_neumann_entropy_sv(sv, list(range(n // 2)), n)

    # Count 2-qubit gates
    n_2q_gates = sum(1 for inst in qc.data if inst.operation.num_qubits == 2)

    # Total MI
    total_mi = 0
    for i in range(n):
        for j in range(i + 1, n):
            total_mi += mutual_information_sv(sv, i, j, n)

    # OTOC from |0...0⟩ through this circuit
    U_mat = Operator(qc).data
    psi0 = np.zeros(2**n)
    psi0[0] = 1.0
    avg_otoc = compute_otoc_from_sv(psi0, U_mat, n)

    # Entanglement per gate
    ent_per_gate = half_ent / n_2q_gates if n_2q_gates > 0 else 0

    state_metrics[name] = {
        'n_2q_gates': n_2q_gates,
        'half_entropy': float(half_ent),
        'half_entropy_fraction': float(half_ent / max_ent),
        'total_mi': float(total_mi),
        'avg_otoc': avg_otoc,
        'entropy_per_gate': float(ent_per_gate),
    }

    print(f"  {name:12s}: gates={n_2q_gates}, S={half_ent:.3f} ({half_ent/max_ent:.1%} of max), "
          f"MI={total_mi:.3f}, OTOC={avg_otoc:.3f}, S/gate={ent_per_gate:.3f}")

# Part 2: Post-preparation scrambling — add random layers after preparation
print("\n=== Part 2: Post-preparation scrambling (adding random layers) ===")
n_extra_layers = 10
n_instances = 8

post_scramble = {}
for name, base_qc in states.items():
    layer_data = []
    for extra in range(n_extra_layers + 1):
        entropies = []
        otocs = []
        for inst in range(n_instances):
            qc = base_qc.copy()
            add_random_layers(qc, n, extra, seed=42 + inst * 1000)
            sv = Statevector.from_instruction(qc)
            half_ent = von_neumann_entropy_sv(sv, list(range(n // 2)), n)
            entropies.append(half_ent)

            U_mat = Operator(qc).data
            psi0 = np.zeros(2**n)
            psi0[0] = 1.0
            otocs.append(compute_otoc_from_sv(psi0, U_mat, n))

        layer_data.append({
            'extra_layers': extra,
            'half_entropy_mean': float(np.mean(entropies)),
            'half_entropy_std': float(np.std(entropies)),
            'otoc_mean': float(np.mean(otocs)),
            'otoc_std': float(np.std(otocs)),
        })

    post_scramble[name] = layer_data
    # Print first and last
    print(f"  {name:12s}: S(+0)={layer_data[0]['half_entropy_mean']:.3f} → "
          f"S(+{n_extra_layers})={layer_data[-1]['half_entropy_mean']:.3f}  "
          f"OTOC(+0)={layer_data[0]['otoc_mean']:.3f} → "
          f"OTOC(+{n_extra_layers})={layer_data[-1]['otoc_mean']:.3f}")

# Also do plain random circuit (no preparation) for comparison
print("\n  Random (no prep):")
random_layer_data = []
for layers in range(n_extra_layers + 1):
    entropies = []
    otocs = []
    for inst in range(n_instances):
        qc = QuantumCircuit(n)
        add_random_layers(qc, n, layers, seed=42 + inst * 1000)
        if layers > 0:
            sv = Statevector.from_instruction(qc)
        else:
            sv = Statevector.from_label('0' * n)
        half_ent = von_neumann_entropy_sv(sv, list(range(n // 2)), n)
        entropies.append(half_ent)

        if layers > 0:
            U_mat = Operator(qc).data
        else:
            U_mat = np.eye(2**n)
        psi0 = np.zeros(2**n)
        psi0[0] = 1.0
        otocs.append(compute_otoc_from_sv(psi0, U_mat, n))

    random_layer_data.append({
        'extra_layers': layers,
        'half_entropy_mean': float(np.mean(entropies)),
        'half_entropy_std': float(np.std(entropies)),
        'otoc_mean': float(np.mean(otocs)),
        'otoc_std': float(np.std(otocs)),
    })

post_scramble['Random'] = random_layer_data
print(f"  {'Random':12s}: S(+0)={random_layer_data[0]['half_entropy_mean']:.3f} → "
      f"S(+{n_extra_layers})={random_layer_data[-1]['half_entropy_mean']:.3f}  "
      f"OTOC(+0)={random_layer_data[0]['otoc_mean']:.3f} → "
      f"OTOC(+{n_extra_layers})={random_layer_data[-1]['otoc_mean']:.3f}")

# Part 3: Compare scrambling time (layers to reach 90% of Page entropy)
print("\n=== Part 3: Time to 90% Page entropy ===")
d_A = 2**(n // 2)
d_B = 2**(n - n // 2)
s_page = np.log2(d_A) - d_A / (2 * np.log(2) * d_A * d_B)
threshold = 0.9 * s_page

for name, data in post_scramble.items():
    t90 = None
    for d in data:
        if d['half_entropy_mean'] >= threshold:
            t90 = d['extra_layers']
            break
    # For structured states, add their prep depth
    n_prep_gates = states[name].count_ops().get('cx', 0) + states[name].count_ops().get('cz', 0) if name != 'Random' else 0
    print(f"  {name:12s}: t_90={t90 if t90 is not None else '>' + str(n_extra_layers)} extra layers "
          f"(prep has {n_prep_gates} 2q gates)")

results = {
    'n_qubits': n,
    'page_entropy': float(s_page),
    'state_metrics': state_metrics,
    'post_scramble': post_scramble,
}

with open('results/sprint_012c_structured_vs_random.json', 'w') as f:
    json.dump(results, f, indent=2)
print("\nResults saved to results/sprint_012c_structured_vs_random.json")
