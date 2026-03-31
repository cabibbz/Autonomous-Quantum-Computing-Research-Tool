"""Sprint 014b: QEC error correction performance.

Apply single-qubit errors at varying rates, measure logical fidelity
with and without correction. Find break-even point.

Tests both the 3-qubit bit-flip code and [[5,1,3]] perfect code.
Uses density matrix simulation with noise channels.
"""

import numpy as np
import json
from itertools import combinations

# Pauli matrices
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)

def kron_list(mats):
    result = mats[0]
    for m in mats[1:]:
        result = np.kron(result, m)
    return result

def pauli_on_qubit(pauli, qubit, n):
    """Apply single Pauli to qubit in n-qubit system."""
    ops = [I2] * n
    ops[qubit] = pauli
    return kron_list(ops)

def depolarizing_channel(rho, p, n):
    """Apply depolarizing noise with probability p to each qubit independently."""
    for q in range(n):
        XI = pauli_on_qubit(X, q, n)
        YI = pauli_on_qubit(Y, q, n)
        ZI = pauli_on_qubit(Z, q, n)
        rho = (1 - p) * rho + (p / 3) * (XI @ rho @ XI + YI @ rho @ YI + ZI @ rho @ ZI)
    return rho

def bitflip_channel(rho, p, n):
    """Apply bit-flip noise with probability p to each qubit independently."""
    for q in range(n):
        XI = pauli_on_qubit(X, q, n)
        rho = (1 - p) * rho + p * (XI @ rho @ XI)
    return rho

def pauli_string(paulis, n):
    """Build n-qubit Pauli from string like 'XZZXI'."""
    pmap = {'I': I2, 'X': X, 'Y': Y, 'Z': Z}
    ops = [pmap[c] for c in paulis]
    return kron_list(ops)


# ============================================================
# 3-qubit bit-flip code
# ============================================================
def three_qubit_encode(logical_rho):
    """Encode 1-qubit logical state into 3-qubit code.
    |0>_L = |000>, |1>_L = |111>
    """
    n = 3
    dim = 2**n
    # Encoding isometry: |0> -> |000>, |1> -> |111>
    V = np.zeros((dim, 2), dtype=complex)
    V[0b000, 0] = 1.0  # |0> -> |000>
    V[0b111, 1] = 1.0  # |1> -> |111>
    return V @ logical_rho @ V.conj().T

def three_qubit_syndrome(rho):
    """Measure syndromes Z1Z2 and Z2Z3, apply correction."""
    n = 3
    # Syndrome operators
    Z0Z1 = pauli_on_qubit(Z, 0, n) @ pauli_on_qubit(Z, 1, n)
    Z1Z2 = pauli_on_qubit(Z, 1, n) @ pauli_on_qubit(Z, 2, n)

    # Projectors for syndrome outcomes
    P_pp = (np.eye(8) + Z0Z1) / 2 @ (np.eye(8) + Z1Z2) / 2  # no error
    P_pm = (np.eye(8) + Z0Z1) / 2 @ (np.eye(8) - Z1Z2) / 2  # error on q2
    P_mp = (np.eye(8) - Z0Z1) / 2 @ (np.eye(8) + Z1Z2) / 2  # error on q0
    P_mm = (np.eye(8) - Z0Z1) / 2 @ (np.eye(8) - Z1Z2) / 2  # error on q1

    # Apply correction based on syndrome
    X0 = pauli_on_qubit(X, 0, n)
    X1 = pauli_on_qubit(X, 1, n)
    X2 = pauli_on_qubit(X, 2, n)

    rho_corrected = (
        P_pp @ rho @ P_pp +            # no error -> no correction
        X2 @ P_pm @ rho @ P_pm @ X2 +  # syndrome 01 -> flip q2
        X0 @ P_mp @ rho @ P_mp @ X0 +  # syndrome 10 -> flip q0
        X1 @ P_mm @ rho @ P_mm @ X1    # syndrome 11 -> flip q1
    )
    return rho_corrected

def three_qubit_decode(rho):
    """Decode: project onto code space, extract logical qubit."""
    V = np.zeros((8, 2), dtype=complex)
    V[0b000, 0] = 1.0
    V[0b111, 1] = 1.0
    return V.conj().T @ rho @ V


# ============================================================
# [[5,1,3]] perfect code
# ============================================================
def build_five_qubit_code():
    """Build encoding isometry for [[5,1,3]] code."""
    n = 5
    dim = 2**n
    stabilizers = ['XZZXI', 'IXZZX', 'XIXZZ', 'ZXIXZ']

    # |0>_L: project |00000> onto +1 eigenspace
    state_0 = np.zeros(dim, dtype=complex)
    state_0[0] = 1.0
    for stab in stabilizers:
        S = pauli_string(stab, n)
        proj = (np.eye(dim) + S) / 2.0
        state_0 = proj @ state_0
        state_0 = state_0 / np.linalg.norm(state_0)

    # |1>_L = X_L |0>_L where X_L = XXXXX
    X_all = pauli_string('XXXXX', n)
    state_1 = X_all @ state_0

    # Encoding isometry
    V = np.zeros((dim, 2), dtype=complex)
    V[:, 0] = state_0
    V[:, 1] = state_1
    return V

def five_qubit_syndrome_correct(rho, V):
    """Syndrome measurement and correction for [[5,1,3]] code."""
    n = 5
    dim = 2**n
    stabilizers = ['XZZXI', 'IXZZX', 'XIXZZ', 'ZXIXZ']
    stab_ops = [pauli_string(s, n) for s in stabilizers]

    # There are 2^4 = 16 syndrome outcomes
    # Each corresponds to a specific error (or no error)
    # For single-qubit errors, syndrome uniquely identifies the error

    # Build syndrome projectors
    rho_corrected = np.zeros_like(rho)
    for syndrome in range(16):
        # Build projector for this syndrome
        proj = np.eye(dim, dtype=complex) / (2**4)
        for k, S in enumerate(stab_ops):
            sign = 1 if (syndrome >> k) & 1 == 0 else -1
            proj = proj @ (np.eye(dim) + sign * S) / 2.0
        # This is wrong — need to build projector differently
        # Use: P_s = (1/2^4) prod_k (I + s_k * S_k)
        proj = np.eye(dim, dtype=complex)
        for k, S in enumerate(stab_ops):
            sign = 1 if (syndrome >> k) & 1 == 0 else -1
            proj = (np.eye(dim) + sign * S) / 2.0 @ proj

        # Find correction operator
        # Syndrome 0000 = no error
        # Single X,Y,Z errors on each qubit give distinct syndromes
        if syndrome == 0:
            correction = np.eye(dim)
        else:
            # Find which single-qubit Pauli gives this syndrome
            found = False
            for q in range(n):
                for P, pname in [(X, 'X'), (Y, 'Y'), (Z, 'Z')]:
                    error = pauli_on_qubit(P, q, n)
                    # Check syndrome of this error
                    s = 0
                    for k, S in enumerate(stab_ops):
                        # Error E anticommutes with S_k if syndrome bit k = 1
                        comm = S @ error - error @ S
                        if np.linalg.norm(comm) > 1e-10:  # anticommutes
                            s |= (1 << k)
                    if s == syndrome:
                        correction = error  # Pauli errors are self-inverse
                        found = True
                        break
                if found:
                    break
            if not found:
                correction = np.eye(dim)  # multi-qubit error, can't correct

        # Project, correct
        projected = proj @ rho @ proj.conj().T
        rho_corrected += correction @ projected @ correction.conj().T

    return rho_corrected


# ============================================================
# Fidelity measurement
# ============================================================
def fidelity(rho, target_state):
    """Fidelity between density matrix and pure state."""
    return float(np.real(target_state.conj() @ rho @ target_state))


# ============================================================
# Run experiments
# ============================================================
results = {}

# Error rates to test
p_values = np.linspace(0, 0.5, 26)

# Logical states to test
logical_states = {
    'zero': np.array([1, 0], dtype=complex),
    'one': np.array([0, 1], dtype=complex),
    'plus': np.array([1, 1], dtype=complex) / np.sqrt(2),
}

# Build 5-qubit code
V5 = build_five_qubit_code()

print("=== 3-qubit bit-flip code (bit-flip noise) ===")
for lname, logical_sv in logical_states.items():
    logical_rho = np.outer(logical_sv, logical_sv.conj())

    fid_uncoded = []
    fid_coded_uncorrected = []
    fid_coded_corrected = []

    for p in p_values:
        # Uncoded: single physical qubit under bit-flip noise
        rho_phys = np.array(logical_rho, dtype=complex)
        rho_phys_noisy = (1 - p) * rho_phys + p * (X @ rho_phys @ X)
        f_uncoded = fidelity(rho_phys_noisy, logical_sv)

        # Coded: encode, apply noise, decode (no correction)
        rho_coded = three_qubit_encode(logical_rho)
        rho_noisy = bitflip_channel(rho_coded, p, 3)
        rho_decoded = three_qubit_decode(rho_noisy)
        f_uncorrected = fidelity(rho_decoded, logical_sv)

        # Coded + corrected
        rho_corrected = three_qubit_syndrome(rho_noisy)
        rho_decoded_c = three_qubit_decode(rho_corrected)
        f_corrected = fidelity(rho_decoded_c, logical_sv)

        fid_uncoded.append(round(float(f_uncoded), 6))
        fid_coded_uncorrected.append(round(float(f_uncorrected), 6))
        fid_coded_corrected.append(round(float(f_corrected), 6))

    results[f'3qubit_bitflip_{lname}'] = {
        'p_values': [round(float(p), 4) for p in p_values],
        'fid_uncoded': fid_uncoded,
        'fid_coded_uncorrected': fid_coded_uncorrected,
        'fid_coded_corrected': fid_coded_corrected
    }
    # Find break-even
    for i, p in enumerate(p_values):
        if fid_coded_corrected[i] < fid_uncoded[i] and p > 0:
            print(f"  |{lname}>: break-even at p≈{p:.3f} (coded={fid_coded_corrected[i]:.4f}, uncoded={fid_uncoded[i]:.4f})")
            break
    else:
        print(f"  |{lname}>: coded always better (or never better)")

    print(f"    p=0.1: uncoded={fid_uncoded[5]:.4f}, corrected={fid_coded_corrected[5]:.4f}")
    print(f"    p=0.3: uncoded={fid_uncoded[15]:.4f}, corrected={fid_coded_corrected[15]:.4f}")


print("\n=== [[5,1,3]] code (depolarizing noise) ===")
for lname, logical_sv in logical_states.items():
    logical_rho = np.outer(logical_sv, logical_sv.conj())

    fid_uncoded = []
    fid_coded_corrected = []

    for p in p_values:
        # Uncoded: single physical qubit under depolarizing noise
        rho_phys = np.array(logical_rho, dtype=complex)
        rho_phys_noisy = (1 - p) * rho_phys + (p / 3) * (X @ rho_phys @ X + Y @ rho_phys @ Y + Z @ rho_phys @ Z)
        f_uncoded = fidelity(rho_phys_noisy, logical_sv)

        # Coded + corrected
        rho_coded = V5 @ logical_rho @ V5.conj().T
        rho_noisy = depolarizing_channel(rho_coded, p, 5)
        rho_corrected = five_qubit_syndrome_correct(rho_noisy, V5)
        rho_decoded = V5.conj().T @ rho_corrected @ V5
        # Normalize
        tr = np.trace(rho_decoded).real
        if tr > 1e-12:
            rho_decoded = rho_decoded / tr
        f_corrected = fidelity(rho_decoded, logical_sv)

        fid_uncoded.append(round(float(f_uncoded), 6))
        fid_coded_corrected.append(round(float(f_corrected), 6))

    results[f'5qubit_depol_{lname}'] = {
        'p_values': [round(float(p), 4) for p in p_values],
        'fid_uncoded': fid_uncoded,
        'fid_coded_corrected': fid_coded_corrected
    }

    # Find break-even
    for i, p in enumerate(p_values):
        if fid_coded_corrected[i] < fid_uncoded[i] and p > 0:
            print(f"  |{lname}>: break-even at p≈{p:.3f} (coded={fid_coded_corrected[i]:.4f}, uncoded={fid_uncoded[i]:.4f})")
            break
    else:
        print(f"  |{lname}>: coded always better up to p=0.5")

    print(f"    p=0.1: uncoded={fid_uncoded[5]:.4f}, corrected={fid_coded_corrected[5]:.4f}")
    print(f"    p=0.3: uncoded={fid_uncoded[15]:.4f}, corrected={fid_coded_corrected[15]:.4f}")


# === 3-qubit code under depolarizing (should fail — can't correct phase errors) ===
print("\n=== 3-qubit bit-flip code under DEPOLARIZING noise (expected to fail) ===")
logical_sv = logical_states['plus']
logical_rho = np.outer(logical_sv, logical_sv.conj())

fid_3q_depol_uncoded = []
fid_3q_depol_corrected = []

for p in p_values:
    # Uncoded
    rho_phys = np.array(logical_rho, dtype=complex)
    rho_phys_noisy = (1 - p) * rho_phys + (p/3) * (X @ rho_phys @ X + Y @ rho_phys @ Y + Z @ rho_phys @ Z)
    f_uncoded = fidelity(rho_phys_noisy, logical_sv)

    # Coded + corrected
    rho_coded = three_qubit_encode(logical_rho)
    rho_noisy = depolarizing_channel(rho_coded, p, 3)
    rho_corrected = three_qubit_syndrome(rho_noisy)
    rho_decoded = three_qubit_decode(rho_corrected)
    tr = np.trace(rho_decoded).real
    if tr > 1e-12:
        rho_decoded = rho_decoded / tr
    f_corrected = fidelity(rho_decoded, logical_sv)

    fid_3q_depol_uncoded.append(round(float(f_uncoded), 6))
    fid_3q_depol_corrected.append(round(float(f_corrected), 6))

results['3qubit_depol_plus'] = {
    'p_values': [round(float(p), 4) for p in p_values],
    'fid_uncoded': fid_3q_depol_uncoded,
    'fid_coded_corrected': fid_3q_depol_corrected
}

print(f"  p=0.1: uncoded={fid_3q_depol_uncoded[5]:.4f}, corrected={fid_3q_depol_corrected[5]:.4f}")
print(f"  p=0.3: uncoded={fid_3q_depol_uncoded[15]:.4f}, corrected={fid_3q_depol_corrected[15]:.4f}")

# Save
with open('results/sprint_014b_qec_performance.json', 'w') as f:
    json.dump(results, f, indent=2)
print("\nResults saved to results/sprint_014b_qec_performance.json")
