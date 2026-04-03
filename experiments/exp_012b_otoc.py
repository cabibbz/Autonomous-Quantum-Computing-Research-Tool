"""
Sprint 012b: Out-of-Time-Ordered Correlators (OTOCs)

OTOC measures how a local perturbation spreads through a quantum system.
F(t) = |⟨ψ| W†(t) V† W(t) V |ψ⟩|²

where W(t) = U†WU is the Heisenberg-evolved operator.

When F(t) = 1: W and V commute → information hasn't spread from V to W
When F(t) → 0: complete scrambling → information from V has reached W

We compute F(t) for random circuits (time = circuit depth) with:
- W = Z on qubit 0 (local perturbation site)
- V = Z on qubit j (probe site at varying distance)

Compare 1D brick-wall vs all-to-all, tracking how distance affects scrambling time.
"""

import numpy as np
import json
from qiskit.circuit import QuantumCircuit
from qiskit.quantum_info import Statevector, Operator, random_unitary

def build_random_circuit_1d(n, n_layers, seed):
    """Build a 1D brick-wall random circuit."""
    qc = QuantumCircuit(n)
    rng = np.random.default_rng(seed)
    for l in range(n_layers):
        offset = l % 2
        for i in range(offset, n - 1, 2):
            u = random_unitary(4, seed=int(rng.integers(0, 2**31)))
            qc.unitary(u, [i, i + 1])
    return qc

def build_random_circuit_alltoall(n, n_layers, seed):
    """Build an all-to-all random circuit."""
    qc = QuantumCircuit(n)
    rng = np.random.default_rng(seed)
    for l in range(n_layers):
        available = list(range(n))
        rng.shuffle(available)
        for idx in range(0, n - 1, 2):
            u = random_unitary(4, seed=int(rng.integers(0, 2**31)))
            qc.unitary(u, [available[idx], available[idx + 1]])
    return qc

def compute_otoc(n, circuit, w_qubit, v_qubit):
    """
    Compute OTOC F(t) = |⟨0| W†(t) V† W(t) V |0⟩|²

    Using the identity: F(t) = |⟨0| U† W† U V† U† W U V |0⟩|²

    We compute this by evolving state vectors:
    |ψ1⟩ = V |0⟩
    |ψ2⟩ = U |ψ1⟩ = U V |0⟩
    |ψ3⟩ = W |ψ2⟩ = W U V |0⟩
    |ψ4⟩ = U† |ψ3⟩ = U† W U V |0⟩
    |ψ5⟩ = V† |ψ4⟩ = V† U† W U V |0⟩

    Similarly:
    |φ1⟩ = U |0⟩
    |φ2⟩ = W |φ1⟩ = W U |0⟩
    |φ3⟩ = U† |φ2⟩ = U† W U |0⟩

    F(t) = |⟨φ3|ψ5⟩|² ... actually let me just do it directly.

    Simpler: compute via operator if small enough, or via state overlaps.
    For n=6 (dim=64), operator approach is fine.
    """
    # Build U operator from circuit
    U = Operator(circuit)
    Udag = U.adjoint()

    # Build W = Z on w_qubit (n-qubit operator)
    W_list = [np.eye(2)] * n
    W_list[w_qubit] = np.array([[1, 0], [0, -1]])  # Z
    W = W_list[0]
    for i in range(1, n):
        W = np.kron(W, W_list[i])

    # Build V = Z on v_qubit
    V_list = [np.eye(2)] * n
    V_list[v_qubit] = np.array([[1, 0], [0, -1]])  # Z
    V = V_list[0]
    for i in range(1, n):
        V = np.kron(V, V_list[i])

    # State approach:
    # F(t) = |⟨0^n| U† W U V† U† W U V |0^n⟩|² ... no, that's the wrong order
    #
    # Actually, OTOC = ⟨W†(t) V† W(t) V⟩ where W(t) = U† W U
    # = ⟨0| (U†WU)† V† (U†WU) V |0⟩
    # = ⟨0| U†W†U V† U†WU V |0⟩
    #
    # |a⟩ = V |0⟩
    # |b⟩ = U |a⟩ = U V |0⟩
    # |c⟩ = W |b⟩ = W U V |0⟩
    # |d⟩ = U† |c⟩ = U† W U V |0⟩
    # |e⟩ = V† |d⟩ = V† U† W U V |0⟩
    #
    # |f⟩ = U |0⟩
    # |g⟩ = W† |f⟩ = W† U |0⟩  (W†=W for Z)
    # |h⟩ = U† |g⟩ = U† W† U |0⟩
    #
    # F = |⟨h|e⟩|² ... wait, let me be more careful.
    #
    # ⟨0| U†W†U V† U†WU V |0⟩ = ⟨h|e⟩ where:
    # |e⟩ = V† U† W U V |0⟩
    # ⟨h| = ⟨0| U† W† U ... so |h⟩ = U† W U |0⟩  (since W†=W for Pauli)
    #
    # Actually since Z†=Z, W†=W.
    # |h⟩ = U† W U |0⟩
    # But ⟨0| U†W†U = (U†W†U|0⟩)† so |h⟩ = U†WU|0⟩

    # Let me just do it with numpy arrays directly
    U_mat = U.data
    Udag_mat = Udag.data

    psi0 = np.zeros(2**n)
    psi0[0] = 1.0

    # |e⟩ = V† U† W U V |0⟩  (V†=V for Pauli Z)
    e = V @ psi0             # V|0⟩
    e = U_mat @ e            # U V|0⟩
    e = W @ e                # W U V|0⟩
    e = Udag_mat @ e         # U† W U V|0⟩
    e = V @ e                # V† U† W U V|0⟩ (V†=V)

    # |h⟩ = U† W U |0⟩
    h = U_mat @ psi0         # U|0⟩
    h = W @ h                # W U|0⟩
    h = Udag_mat @ h         # U† W U|0⟩

    F = np.abs(np.dot(h.conj(), e))**2
    return float(F)


# Parameters
n = 6
max_layers = 12
n_instances = 10
w_qubit = 0  # perturbation site

results = {
    'n_qubits': n,
    'max_layers': max_layers,
    'n_instances': n_instances,
    'w_qubit': w_qubit,
    'geometries': {}
}

for geom_name, builder in [('1D_brickwall', build_random_circuit_1d),
                             ('all_to_all', build_random_circuit_alltoall)]:
    print(f"\n=== {geom_name} ===")
    geom_results = {}

    for v_qubit in range(1, n):  # probe at each other qubit
        distance = v_qubit  # for 1D, this is the chain distance
        layer_otocs = []

        for layer in range(max_layers + 1):
            f_values = []
            for inst in range(n_instances):
                seed = 42 + inst * 1000
                qc = builder(n, layer, seed)
                f_val = compute_otoc(n, qc, w_qubit, v_qubit)
                f_values.append(f_val)

            avg_f = float(np.mean(f_values))
            std_f = float(np.std(f_values))
            layer_otocs.append({
                'layer': layer,
                'F_mean': avg_f,
                'F_std': std_f,
            })

        geom_results[f'v_qubit_{v_qubit}'] = {
            'v_qubit': v_qubit,
            'distance': distance,
            'layers': layer_otocs
        }

        # Print summary for this probe qubit
        f_values_by_layer = [d['F_mean'] for d in layer_otocs]
        # Find scrambling time: first layer where F < 0.5
        scramble_time = None
        for i, f in enumerate(f_values_by_layer):
            if f < 0.5:
                scramble_time = i
                break
        print(f"  V=Z_{v_qubit} (dist={distance}): "
              f"F(0)={f_values_by_layer[0]:.3f} → F({max_layers})={f_values_by_layer[-1]:.3f}  "
              f"t_scr={'>' + str(max_layers) if scramble_time is None else scramble_time}")

    results['geometries'][geom_name] = geom_results

# Save
with open('results/sprint_012b_otoc.json', 'w') as f:
    json.dump(results, f, indent=2)
print("\nResults saved to results/sprint_012b_otoc.json")
