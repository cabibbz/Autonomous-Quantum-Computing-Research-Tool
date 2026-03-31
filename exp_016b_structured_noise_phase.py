"""Sprint 016b: Phase damping (T2 / dephasing) on QEC code families.

[[5,1,3]] with full correction. Steane with projection fidelity.
Tests |0>_L and |+>_L — dephasing should be basis-dependent:
|0> immune (Z eigenstate), |+> vulnerable.
"""

import numpy as np
import json
import time

I2 = np.eye(2, dtype=complex)
X = np.array([[0,1],[1,0]], dtype=complex)
Y = np.array([[0,-1j],[1j,0]], dtype=complex)
Z = np.array([[1,0],[0,-1]], dtype=complex)
pauli_map = {'I': I2, 'X': X, 'Y': Y, 'Z': Z}

def kron_list(mats):
    r = mats[0]
    for m in mats[1:]: r = np.kron(r, m)
    return r

def pauli_on_qubit(p, q, n):
    ops = [I2]*n; ops[q] = p; return kron_list(ops)

def pauli_string(s):
    return kron_list([pauli_map[c] for c in s])

def fidelity(rho, sv):
    return float(np.real(sv.conj() @ rho @ sv))

def phase_damping_channel(rho, lam, n):
    """Phase damping: E0 = diag(1, sqrt(1-lam)), E1 = diag(0, sqrt(lam))"""
    for q in range(n):
        E0 = np.array([[1,0],[0,np.sqrt(1-lam)]], dtype=complex)
        E1 = np.array([[0,0],[0,np.sqrt(lam)]], dtype=complex)
        ops0 = [I2]*n; ops0[q] = E0
        ops1 = [I2]*n; ops1[q] = E1
        E0f = kron_list(ops0); E1f = kron_list(ops1)
        rho = E0f @ rho @ E0f.conj().T + E1f @ rho @ E1f.conj().T
    return rho

def encode_stabilizer(stabilizers, n):
    dim = 2**n; state = np.zeros(dim, dtype=complex); state[0] = 1.0
    for s in stabilizers:
        S = pauli_string(s)
        state = (np.eye(dim) + S) / 2.0 @ state
        nm = np.linalg.norm(state)
        if nm < 1e-12: raise ValueError(f"{s} annihilated!")
        state /= nm
    return state

# ============================================================
# Build codes
# ============================================================
print("=== Building codes ===")

stab5 = ['XZZXI','IXZZX','XIXZZ','ZXIXZ']
sv5_0 = encode_stabilizer(stab5, 5)
sv5_1 = pauli_string('XXXXX') @ sv5_0; sv5_1 /= np.linalg.norm(sv5_1)
V5 = np.zeros((32,2), dtype=complex); V5[:,0]=sv5_0; V5[:,1]=sv5_1

# Also build |+>_L
sv5_plus = (sv5_0 + sv5_1) / np.sqrt(2)

# Correction ops for [[5,1,3]]
stab_ops5 = [pauli_string(s) for s in stab5]
table5 = {0: np.eye(32, dtype=complex)}
for q in range(5):
    for P in [X,Y,Z]:
        error = pauli_on_qubit(P, q, 5)
        s = 0
        for k, S in enumerate(stab_ops5):
            if np.linalg.norm(S @ error - error @ S) > 1e-10: s |= (1<<k)
        if s not in table5: table5[s] = error
pc5 = []
for syn in range(16):
    proj = np.eye(32, dtype=complex)
    for k, S in enumerate(stab_ops5):
        sign = 1 if (syn>>k)&1==0 else -1
        proj = (np.eye(32)+sign*S)/2.0 @ proj
    if np.linalg.norm(proj) < 1e-12: continue
    pc5.append(table5.get(syn, np.eye(32, dtype=complex)) @ proj)
print(f"[[5,1,3]]: {len(pc5)} correction ops")

# Steane
stab7 = ['IIIXXXX','IXXIIXX','XIXIXIX','IIIZZZZ','IZZIIZZ','ZIZIZIZ']
sv7_0 = encode_stabilizer(stab7, 7)
sv7_1 = pauli_string('XXXXXXX') @ sv7_0; sv7_1 /= np.linalg.norm(sv7_1)
V7 = np.zeros((128,2), dtype=complex); V7[:,0]=sv7_0; V7[:,1]=sv7_1
sv7_plus = (sv7_0 + sv7_1) / np.sqrt(2)
print("Steane: projection")

# ============================================================
# Experiment
# ============================================================
lam_values = np.array([0.0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5])
results = {}

for code_name, V, n, pc, sv_plus in [('five_qubit',V5,5,pc5,sv5_plus),
                                       ('steane',V7,7,None,sv7_plus)]:
    has_corr = pc is not None
    print(f"\n=== {code_name} (n={n}) under phase damping ===")

    # Test |0>_L (Z eigenstate — should be immune)
    # and |+>_L (X eigenstate — should be vulnerable)
    states = [
        ('zero', np.array([1,0], dtype=complex)),
        ('plus', np.array([1,1], dtype=complex)/np.sqrt(2)),
    ]

    for lname, logical_sv in states:
        logical_rho = np.outer(logical_sv, logical_sv.conj())
        fid_unc, fid_cod, ret = [], [], []
        t0 = time.time()

        for lam in lam_values:
            # Uncoded single qubit
            E0 = np.array([[1,0],[0,np.sqrt(1-lam)]], dtype=complex)
            E1 = np.array([[0,0],[0,np.sqrt(lam)]], dtype=complex)
            rho_u = E0 @ logical_rho @ E0.conj().T + E1 @ logical_rho @ E1.conj().T
            fid_unc.append(round(fidelity(rho_u, logical_sv), 6))

            # Coded
            rho_coded = V @ logical_rho @ V.conj().T
            rho_noisy = phase_damping_channel(rho_coded, lam, n)

            if has_corr:
                rho_c = sum(cp @ rho_noisy @ cp.conj().T for cp in pc)
                rho_dec = V.conj().T @ rho_c @ V
            else:
                rho_dec = V.conj().T @ rho_noisy @ V

            tr = np.trace(rho_dec).real
            ret.append(round(float(tr), 6))
            if tr > 1e-12: rho_dec /= tr
            fid_cod.append(round(fidelity(rho_dec, logical_sv), 6))

        elapsed = time.time()-t0
        results[f'{code_name}_{lname}'] = {
            'n': n, 'noise': 'phase_damping', 'has_correction': has_corr,
            'lambda_values': [round(float(l),4) for l in lam_values],
            'fid_uncoded': fid_unc, 'fid_coded': fid_cod,
            'codespace_retention': ret,
        }
        be = next((round(float(l),4) for i,l in enumerate(lam_values) if l>0 and fid_cod[i]<fid_unc[i]), None)
        print(f"  |{lname}>: break-even ~{be}, {elapsed:.1f}s")

# ============================================================
# Summary
# ============================================================
print("\n\n=== PHASE DAMPING: |zero> ===")
print(f"{'lambda':<8} {'[[5,1,3]]*':<14} {'Steane':<14} {'uncoded':<12}")
print("-" * 48)
for i,l in enumerate(lam_values):
    print(f"{l:.3f}   {results['five_qubit_zero']['fid_coded'][i]:.4f}        "
          f"{results['steane_zero']['fid_coded'][i]:.4f}        "
          f"{results['five_qubit_zero']['fid_uncoded'][i]:.4f}")
print("* = with correction")

print("\n=== PHASE DAMPING: |plus> ===")
print(f"{'lambda':<8} {'[[5,1,3]]*':<14} {'Steane':<14} {'uncoded':<12}")
print("-" * 48)
for i,l in enumerate(lam_values):
    print(f"{l:.3f}   {results['five_qubit_plus']['fid_coded'][i]:.4f}        "
          f"{results['steane_plus']['fid_coded'][i]:.4f}        "
          f"{results['five_qubit_plus']['fid_uncoded'][i]:.4f}")

print("\n=== STATE ASYMMETRY (|0> - |+> fidelity) ===")
for cn in ['five_qubit','steane']:
    asym = [results[f'{cn}_zero']['fid_coded'][i] - results[f'{cn}_plus']['fid_coded'][i]
            for i in range(len(lam_values))]
    mi = int(np.argmax([abs(a) for a in asym]))
    print(f"  {cn}: max |asym| = {abs(asym[mi]):.4f} at lambda={lam_values[mi]:.3f}")

print("\n=== CODESPACE RETENTION ===")
for cn in ['five_qubit','steane']:
    for st in ['zero','plus']:
        r = results[f'{cn}_{st}']['codespace_retention']
        print(f"  {cn} |{st}>: l=0.05:{r[1]:.4f}, l=0.2:{r[4]:.4f}, l=0.5:{r[-1]:.4f}")

# Compare amp damping asymmetry vs phase damping asymmetry
print("\n=== NOISE ASYMMETRY COMPARISON ===")
print("  Amplitude damping: biased in logical basis (|0> vs |1>)")
print("  Phase damping: biased in conjugate basis (|0> vs |+>)")
print("  If code is symmetric under both: it's a true quantum code")
print("  If code is only symmetric under one: it has a preferred basis")

with open('results/sprint_016b_structured_noise_phase.json', 'w') as f:
    json.dump(results, f, indent=2)
print("\nResults saved.")
