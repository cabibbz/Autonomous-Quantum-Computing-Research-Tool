"""
Sprint 019a: Coherent Information of Quantum Channels.

Coherent information I_coh(ρ, N) = S(N(ρ)) - S((I⊗N)(|ψ⟩⟨ψ|))
where |ψ⟩ is a purification of ρ.

For a qubit channel N acting on 1 qubit:
- Prepare maximally entangled state |Φ⟩ = (|00⟩+|11⟩)/√2 on AB (reference + system)
- Apply N to system qubit B only
- I_coh = S(B) - S(AB)  [equivalently S(B) - S(E) via complementary channel]

The quantum capacity Q(N) = lim (1/n) max I_coh(ρ^n, N^⊗n)
For degradable channels (amplitude damping): Q = max_ρ I_coh(ρ, N) (single-letter)
For depolarizing: Q ≥ max_ρ I_coh(ρ, N) (hashing bound, may not be tight)

We compute coherent information for:
1. Depolarizing channel: ρ → (1-p)ρ + p/3(XρX + YρY + ZρZ)
2. Amplitude damping: ρ → E0 ρ E0† + E1 ρ E1†
3. Phase damping (dephasing): ρ → (1-λ)ρ + λ ZρZ
"""

import numpy as np
import json

# Pauli matrices
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)

def von_neumann_entropy(rho):
    """Von Neumann entropy in bits."""
    evals = np.linalg.eigvalsh(rho)
    evals = evals[evals > 1e-15]
    return -np.sum(evals * np.log2(evals))

def partial_trace(rho, dims, keep):
    """Partial trace of a density matrix.
    dims: list of subsystem dimensions
    keep: list of indices to keep
    """
    n = len(dims)
    rho_r = rho.reshape(dims + dims)
    # Trace over all indices not in keep
    trace_over = sorted(set(range(n)) - set(keep))
    for i, idx in enumerate(trace_over):
        # Trace pairs: idx and idx+n (adjusted for previous traces)
        adj_idx = idx - i  # axes shift as we trace
        adj_n = n - i
        rho_r = np.trace(rho_r, axis1=adj_idx, axis2=adj_idx + adj_n)
    return rho_r.reshape(np.prod([dims[k] for k in keep]), -1)

def apply_channel_to_qubit_B(rho_AB, kraus_ops):
    """Apply channel to qubit B of a 2-qubit state rho_AB.
    kraus_ops: list of 2x2 Kraus operators for the channel.
    """
    d = 4  # 2 qubits
    result = np.zeros((d, d), dtype=complex)
    for K in kraus_ops:
        # Full operator: I_A ⊗ K_B
        full_K = np.kron(I2, K)
        result += full_K @ rho_AB @ full_K.conj().T
    return result

def depolarizing_kraus(p):
    """Kraus operators for depolarizing channel with error rate p.
    ρ → (1-p)ρ + (p/3)(XρX + YρY + ZρZ)
    """
    return [
        np.sqrt(1 - p) * I2,
        np.sqrt(p / 3) * X,
        np.sqrt(p / 3) * Y,
        np.sqrt(p / 3) * Z
    ]

def amplitude_damping_kraus(gamma):
    """Kraus operators for amplitude damping channel."""
    E0 = np.array([[1, 0], [0, np.sqrt(1 - gamma)]], dtype=complex)
    E1 = np.array([[0, np.sqrt(gamma)], [0, 0]], dtype=complex)
    return [E0, E1]

def phase_damping_kraus(lam):
    """Kraus operators for phase damping (dephasing) channel.
    ρ → (1-λ)ρ + λ ZρZ
    """
    return [
        np.sqrt(1 - lam) * I2,
        np.sqrt(lam) * Z
    ]

def coherent_info_max_entangled(kraus_ops):
    """Compute coherent information using maximally entangled input |Φ+⟩.
    I_coh = S(B) - S(AB)
    """
    # Maximally entangled state |Φ+⟩ = (|00⟩+|11⟩)/√2
    phi_plus = np.array([1, 0, 0, 1], dtype=complex) / np.sqrt(2)
    rho_AB = np.outer(phi_plus, phi_plus.conj())

    # Apply channel to qubit B
    rho_AB_out = apply_channel_to_qubit_B(rho_AB, kraus_ops)

    # S(B) = entropy of reduced state on B
    rho_B = partial_trace(rho_AB_out, [2, 2], [1])
    S_B = von_neumann_entropy(rho_B)

    # S(AB) = entropy of full output
    S_AB = von_neumann_entropy(rho_AB_out)

    return S_B - S_AB

def coherent_info_optimized(kraus_ops, n_samples=50):
    """Optimize coherent information over input states.
    For qubit channels, the optimal input for coherent info is diagonal
    in the eigenbasis of the channel's output, so we scan over |ψ⟩ = cos(θ)|0⟩ + sin(θ)|1⟩.
    """
    best_icoh = -np.inf
    best_theta = 0

    for theta in np.linspace(0, np.pi / 2, n_samples):
        # Input state
        psi = np.array([np.cos(theta), np.sin(theta)], dtype=complex)
        rho_in = np.outer(psi, psi.conj())

        # Purification: |ψ_AB⟩ = cos(θ)|00⟩ + sin(θ)|11⟩
        psi_AB = np.array([np.cos(theta), 0, 0, np.sin(theta)], dtype=complex)
        rho_AB = np.outer(psi_AB, psi_AB.conj())

        # Apply channel to B
        rho_AB_out = apply_channel_to_qubit_B(rho_AB, kraus_ops)

        rho_B = partial_trace(rho_AB_out, [2, 2], [1])
        S_B = von_neumann_entropy(rho_B)
        S_AB = von_neumann_entropy(rho_AB_out)

        icoh = S_B - S_AB
        if icoh > best_icoh:
            best_icoh = icoh
            best_theta = theta

    return best_icoh, best_theta

# ============================================================
# MAIN EXPERIMENT
# ============================================================

noise_params = np.linspace(0, 1.0, 51)  # 0 to 1 in steps of 0.02

results = {
    'noise_params': noise_params.tolist(),
    'depolarizing': {'icoh_max_ent': [], 'icoh_optimized': [], 'opt_theta': []},
    'amplitude_damping': {'icoh_max_ent': [], 'icoh_optimized': [], 'opt_theta': []},
    'phase_damping': {'icoh_max_ent': [], 'icoh_optimized': [], 'opt_theta': []},
}

print("Computing coherent information for quantum channels...")
print("=" * 60)

for channel_name, kraus_fn in [
    ('depolarizing', depolarizing_kraus),
    ('amplitude_damping', amplitude_damping_kraus),
    ('phase_damping', phase_damping_kraus),
]:
    print(f"\n{channel_name}:")
    for p in noise_params:
        kraus = kraus_fn(p)

        # Coherent info with maximally entangled input
        icoh_me = coherent_info_max_entangled(kraus)

        # Optimized coherent info
        icoh_opt, theta_opt = coherent_info_optimized(kraus)

        results[channel_name]['icoh_max_ent'].append(icoh_me)
        results[channel_name]['icoh_optimized'].append(icoh_opt)
        results[channel_name]['opt_theta'].append(theta_opt)

    # Find zero-crossing (where coherent info becomes negative)
    icoh_opt = results[channel_name]['icoh_optimized']
    zero_cross = None
    for i in range(len(icoh_opt) - 1):
        if icoh_opt[i] > 0 and icoh_opt[i + 1] <= 0:
            # Linear interpolation
            p1, p2 = noise_params[i], noise_params[i + 1]
            v1, v2 = icoh_opt[i], icoh_opt[i + 1]
            zero_cross = p1 + (p2 - p1) * v1 / (v1 - v2)
            break

    results[channel_name]['zero_crossing'] = zero_cross
    print(f"  Zero crossing (hashing bound): p = {zero_cross}")
    print(f"  I_coh at p=0: {icoh_opt[0]:.4f}")
    print(f"  I_coh at p=0.1: {icoh_opt[5]:.4f}")

    # Print a few key values
    for idx in [0, 5, 10, 15, 25]:
        p = noise_params[idx]
        print(f"  p={p:.2f}: I_coh(max_ent)={results[channel_name]['icoh_max_ent'][idx]:.4f}, "
              f"I_coh(opt)={icoh_opt[idx]:.4f}, θ_opt={results[channel_name]['opt_theta'][idx]:.3f}")

# Known theoretical values for comparison
print("\n" + "=" * 60)
print("THEORETICAL COMPARISONS:")
print("=" * 60)
print(f"Depolarizing hashing bound (theory): p ≈ 0.1893 (from 1 - H(3p/4) - (3p/4)log₂3)")
print(f"Depolarizing hashing bound (measured): p ≈ {results['depolarizing']['zero_crossing']}")
print(f"Amplitude damping capacity threshold: γ = 0.5 (theory)")
print(f"Amplitude damping zero crossing (measured): γ ≈ {results['amplitude_damping']['zero_crossing']}")
print(f"Phase damping zero crossing (measured): λ ≈ {results['phase_damping']['zero_crossing']}")

# Compute theoretical depolarizing hashing bound
# Q ≥ 1 - H(p) - p*log₂(3) where the channel is ρ→(1-4p/3)ρ + (4p/3)I/2
# For our parametrization: ρ→(1-p)ρ + p/3(XρX+YρY+ZρZ)
# Effective depolarizing: (1-4p/3) fidelity
# Hashing bound: 1 + (1-p)log₂(1-p) + p*log₂(p/3)
def depol_hashing_bound(p):
    """Hashing bound for depolarizing channel.
    Q ≥ 1 + (1-p)*log₂(1-p) + p*log₂(p/3)
    """
    if p <= 0 or p >= 1:
        return 0
    return 1 + (1-p)*np.log2(1-p) + p*np.log2(p/3)

theory_hashing = [max(0, depol_hashing_bound(p)) for p in noise_params]
results['depolarizing']['theory_hashing'] = theory_hashing

# Save results
with open('results/sprint_019a_coherent_information.json', 'w') as f:
    json.dump(results, f, indent=2)
print("\nResults saved to results/sprint_019a_coherent_information.json")

# Summary analysis
print("\n" + "=" * 60)
print("SUMMARY")
print("=" * 60)
print(f"\nChannel capacity thresholds (where Q → 0):")
print(f"  Depolarizing:     p ≈ {results['depolarizing']['zero_crossing']:.4f}")
print(f"  Amplitude damping: γ ≈ {results['amplitude_damping']['zero_crossing']:.4f}")
print(f"  Phase damping:    λ ≈ {results['phase_damping']['zero_crossing']:.4f}")
print(f"\nCapacity at low noise (p=0.02):")
for ch in ['depolarizing', 'amplitude_damping', 'phase_damping']:
    print(f"  {ch}: Q ≈ {results[ch]['icoh_optimized'][1]:.4f} bits/use")
print(f"\nMax-entangled vs optimized input:")
for ch in ['depolarizing', 'amplitude_damping', 'phase_damping']:
    me = results[ch]['icoh_max_ent'][5]  # p=0.1
    opt = results[ch]['icoh_optimized'][5]
    diff = opt - me
    print(f"  {ch} at p=0.1: opt-ME = {diff:.4f} ({'+' if diff > 0 else ''}{diff:.4f})")
