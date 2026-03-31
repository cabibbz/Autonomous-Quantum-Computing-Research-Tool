"""
Quantum Explorer — Sprint 001
Bell States, CHSH Inequality Violation, and Entanglement Entropy Scaling

Questions:
- Do Bell states produce the expected perfect correlations?
- How close to the Tsirelson bound (2√2) can we get with CHSH?
- How does entanglement entropy scale with system size for GHZ states?
- What does the measurement distribution look like for entangled vs product states?
"""

from qiskit.circuit import QuantumCircuit
from qiskit_aer import AerSimulator
from qiskit.quantum_info import Statevector, partial_trace, entropy, DensityMatrix
import numpy as np
import json
import datetime
import os

RESULTS_DIR = "results"
os.makedirs(RESULTS_DIR, exist_ok=True)

def log_result(experiment_name, data):
    """Save experiment results to JSON."""
    timestamp = datetime.datetime.now().isoformat()
    entry = {
        "timestamp": timestamp,
        "experiment": experiment_name,
        "data": data
    }
    log_path = os.path.join(RESULTS_DIR, "experiment_log.jsonl")
    with open(log_path, "a") as f:
        f.write(json.dumps(entry, default=str) + "\n")
    print(f"Logged: {experiment_name} at {timestamp}")
    return entry


# ============================================================
# Experiment 1: Create and verify all 4 Bell states
# ============================================================
def experiment_bell_states():
    """Create Phi+, Phi-, Psi+, Psi- and verify via statevector and sampling."""
    sim = AerSimulator(method='statevector')
    results = {}

    bell_configs = {
        "Phi+": {"x_gate": False, "z_gate": False},  # (|00⟩ + |11⟩)/√2
        "Phi-": {"x_gate": False, "z_gate": True},    # (|00⟩ - |11⟩)/√2
        "Psi+": {"x_gate": True,  "z_gate": False},   # (|01⟩ + |10⟩)/√2
        "Psi-": {"x_gate": True,  "z_gate": True},    # (|01⟩ - |10⟩)/√2
    }

    for name, config in bell_configs.items():
        qc = QuantumCircuit(2)
        qc.h(0)
        qc.cx(0, 1)
        if config["x_gate"]:
            qc.x(0)
        if config["z_gate"]:
            qc.z(0)

        # Get statevector
        sv = Statevector.from_instruction(qc)
        amplitudes = {format(i, '02b'): complex(a) for i, a in enumerate(sv.data)}

        # Measure
        qc_meas = qc.copy()
        qc_meas.measure_all()
        job = sim.run(qc_meas, shots=10000)
        counts = job.result().get_counts()

        # Compute reduced density matrix and entanglement entropy
        dm = DensityMatrix(sv)
        rho_A = partial_trace(dm, [1])
        ent = entropy(rho_A, base=2)

        results[name] = {
            "amplitudes": {k: [v.real, v.imag] for k, v in amplitudes.items()},
            "counts": counts,
            "entanglement_entropy_bits": float(ent),
            "total_shots": sum(counts.values()),
        }

        # Check: should be maximally entangled (entropy = 1 bit)
        print(f"\n{name}:")
        print(f"  Amplitudes: {amplitudes}")
        print(f"  Counts: {counts}")
        print(f"  Entanglement entropy: {ent:.6f} bits (expect 1.0)")

        # Check correlations: in Phi+/Phi- the qubits always agree, in Psi+/Psi- they always disagree
        agreeing = sum(v for k, v in counts.items() if k[0] == k[1])
        disagreeing = sum(v for k, v in counts.items() if k[0] != k[1])
        correlation = (agreeing - disagreeing) / (agreeing + disagreeing)
        results[name]["correlation"] = float(correlation)
        print(f"  Correlation: {correlation:.4f}")

    return results


# ============================================================
# Experiment 2: CHSH inequality test
# ============================================================
def experiment_chsh():
    """
    Measure the CHSH parameter S using optimal measurement angles.
    Classical bound: |S| ≤ 2
    Quantum bound (Tsirelson): |S| ≤ 2√2 ≈ 2.828

    We prepare a Bell state and measure in 4 angle combinations:
    Alice: a=0, a'=π/2
    Bob: b=π/4, b'=-π/4
    """
    sim = AerSimulator(method='statevector')
    shots = 100000

    # Optimal CHSH angles
    alice_angles = [0, np.pi/2]          # a, a'
    bob_angles = [np.pi/4, 3*np.pi/4]    # b, b'

    correlations = {}

    for i, a in enumerate(alice_angles):
        for j, b in enumerate(bob_angles):
            qc = QuantumCircuit(2, 2)
            # Create Bell state Phi+
            qc.h(0)
            qc.cx(0, 1)

            # Alice measures in rotated basis
            qc.ry(-a, 0)
            # Bob measures in rotated basis
            qc.ry(-b, 1)

            qc.measure([0, 1], [0, 1])

            job = sim.run(qc, shots=shots)
            counts = job.result().get_counts()

            # Compute <AB> expectation value
            # outcome 0 -> +1, outcome 1 -> -1
            E = 0
            for bitstring, count in counts.items():
                bits = [int(x) for x in bitstring]
                a_val = 1 - 2*bits[1]  # qubit 0 is rightmost
                b_val = 1 - 2*bits[0]  # qubit 1 is leftmost
                E += a_val * b_val * count
            E /= shots

            label = f"a{i}_b{j}"
            correlations[label] = float(E)
            print(f"  E(a{i},b{j}) = {E:.4f}  [a={a:.3f}, b={b:.3f}]")

    # CHSH: S = E(a,b) - E(a,b') + E(a',b) + E(a',b')
    S = correlations["a0_b0"] - correlations["a0_b1"] + correlations["a1_b0"] + correlations["a1_b1"]

    print(f"\n  CHSH parameter S = {S:.4f}")
    print(f"  Classical bound: 2.0")
    print(f"  Tsirelson bound: {2*np.sqrt(2):.4f}")
    print(f"  Violation: {'YES' if abs(S) > 2 else 'NO'}")

    return {
        "correlations": correlations,
        "S": float(S),
        "classical_bound": 2.0,
        "tsirelson_bound": float(2*np.sqrt(2)),
        "violates_classical": bool(abs(S) > 2),
        "shots_per_setting": shots,
    }


# ============================================================
# Experiment 3: Entanglement entropy scaling for GHZ states
# ============================================================
def experiment_ghz_entropy_scaling():
    """
    Create GHZ states |GHZ_n⟩ = (|00...0⟩ + |11...1⟩)/√2 for n=2..25
    Measure entanglement entropy of subsystem A (first k qubits) for various bipartitions.

    GHZ states are interesting because:
    - They are maximally entangled across ANY bipartition (entropy = 1 bit always)
    - This is DIFFERENT from random states where entropy scales with min(|A|, |B|)
    - This flat entropy profile is a signature of "fragile" entanglement
    """
    results = {"by_n": {}, "scaling_summary": []}

    for n in range(2, 15):  # Cap at 14 qubits for reasonable runtime
        qc = QuantumCircuit(n)
        qc.h(0)
        for i in range(1, n):
            qc.cx(0, i)

        sv = Statevector.from_instruction(qc)

        entropies = {}
        # Check a representative set of bipartitions
        if n <= 10:
            ks = list(range(1, n))
        else:
            ks = [1, n//4, n//2, 3*n//4, n-1]
        for k in ks:
            subsystem_to_trace = list(range(k, n))
            rho_A = partial_trace(sv, subsystem_to_trace)
            ent = float(entropy(rho_A, base=2))
            entropies[f"A={k}_B={n-k}"] = ent

        results["by_n"][str(n)] = entropies

        # For the half-bipartition
        half = n // 2
        half_key = f"A={half}_B={n-half}"
        half_ent = entropies[half_key]
        results["scaling_summary"].append({
            "n": n,
            "half_entropy": half_ent,
            "all_bipartition_entropies": list(entropies.values()),
        })

        ent_vals = list(entropies.values())
        print(f"  n={n:2d}: entropy range [{min(ent_vals):.4f}, {max(ent_vals):.4f}], "
              f"half-cut={half_ent:.4f}")

    # Verify: all entropies should be exactly 1.0 for GHZ
    all_entropies = []
    for n_data in results["scaling_summary"]:
        all_entropies.extend(n_data["all_bipartition_entropies"])

    results["verification"] = {
        "all_entropies_are_1": all(abs(e - 1.0) < 1e-10 for e in all_entropies),
        "mean_entropy": float(np.mean(all_entropies)),
        "std_entropy": float(np.std(all_entropies)),
        "max_n_tested": 14,
    }
    print(f"\n  All GHZ entropies = 1.0? {results['verification']['all_entropies_are_1']}")

    return results


# ============================================================
# Experiment 4: Random states vs GHZ — entropy comparison
# ============================================================
def experiment_random_vs_ghz_entropy():
    """
    Compare entanglement entropy scaling between:
    1. GHZ states: entropy = 1 for all bipartitions (fragile entanglement)
    2. Random Haar states: entropy ≈ min(|A|, |B|) * log(2) (volume law)

    This shows fundamentally different entanglement structures.
    """
    results = {"ghz": [], "random_haar": [], "n_range": list(range(2, 13))}

    for n in range(2, 13):  # Cap at 12 qubits for random state partial traces
        # GHZ half-cut entropy
        qc = QuantumCircuit(n)
        qc.h(0)
        for i in range(1, n):
            qc.cx(0, i)
        sv_ghz = Statevector.from_instruction(qc)
        half = n // 2
        rho_ghz = partial_trace(sv_ghz, list(range(half, n)))
        ent_ghz = float(entropy(rho_ghz, base=2))
        results["ghz"].append(ent_ghz)

        # Average over random Haar states
        n_samples = 10 if n >= 10 else 20
        haar_ents = []
        for _ in range(n_samples):
            sv_rand = Statevector(np.random.randn(2**n) + 1j*np.random.randn(2**n))
            sv_rand = sv_rand / np.linalg.norm(sv_rand.data)
            rho_rand = partial_trace(sv_rand, list(range(half, n)))
            ent_rand = float(entropy(rho_rand, base=2))
            haar_ents.append(ent_rand)
        avg_haar = float(np.mean(haar_ents))
        results["random_haar"].append(avg_haar)

        if n <= 8 or n % 3 == 0:
            print(f"  n={n:2d}: GHZ={ent_ghz:.4f}, Random(avg)={avg_haar:.4f}, "
                  f"theoretical_max={min(half, n-half):.1f}")

    # Page curve: theoretical entropy for random states
    # S ≈ min(k, n-k) * log(2) - 1/(2*ln(2)*2^(|n-2k|+1))
    results["observation"] = (
        "GHZ states have constant 1-bit entropy regardless of system size. "
        "Random states follow volume law — entropy grows with subsystem size. "
        "This means GHZ entanglement is 'all-or-nothing': measuring ANY qubit "
        "collapses the entire state, while random states are robust to local measurement."
    )

    return results


# ============================================================
# Main
# ============================================================
def run_experiment():
    all_results = {}

    print("=" * 60)
    print("SPRINT 001: Bell States, CHSH, & Entanglement Scaling")
    print("=" * 60)

    print("\n--- Experiment 1: Bell States ---")
    all_results["bell_states"] = experiment_bell_states()

    print("\n--- Experiment 2: CHSH Inequality ---")
    all_results["chsh"] = experiment_chsh()

    print("\n--- Experiment 3: GHZ Entropy Scaling ---")
    all_results["ghz_entropy"] = experiment_ghz_entropy_scaling()

    print("\n--- Experiment 4: Random vs GHZ Entropy ---")
    all_results["random_vs_ghz"] = experiment_random_vs_ghz_entropy()

    print("\n" + "=" * 60)
    print("ALL EXPERIMENTS COMPLETE")
    print("=" * 60)

    return all_results


if __name__ == "__main__":
    results = run_experiment()
    if results:
        log_result("sprint_001_bell_chsh_entropy", results)

        # Also save detailed results
        with open(os.path.join(RESULTS_DIR, "sprint_001_results.json"), "w") as f:
            json.dump(results, f, indent=2, default=str)
        print(f"\nDetailed results saved to {RESULTS_DIR}/sprint_001_results.json")
