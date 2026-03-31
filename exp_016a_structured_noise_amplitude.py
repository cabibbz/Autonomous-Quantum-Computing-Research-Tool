"""Sprint 016a: Amplitude damping (T1) on QEC code families.

[[5,1,3]]: full syndrome correction. Steane/Shor: projection fidelity.
Reduced parameter sweep (8 points) to fit time budget.
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

def precompute_kraus_ad(gamma, n):
    """Precompute amplitude damping Kraus operators for all qubits."""
    E0 = np.array([[1,0],[0,np.sqrt(1-gamma)]], dtype=complex)
    E1 = np.array([[0,np.sqrt(gamma)],[0,0]], dtype=complex)
    kraus = []
    for q in range(n):
        ops0 = [I2]*n; ops0[q] = E0
        ops1 = [I2]*n; ops1[q] = E1
        kraus.append((kron_list(ops0), kron_list(ops1)))
    return kraus

def apply_kraus(rho, kraus_list):
    for E0f, E1f in kraus_list:
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
print("Steane: projection")

# Shor
plus_b = np.zeros(8, dtype=complex); plus_b[0]=plus_b[7]=1/np.sqrt(2)
minus_b = np.zeros(8, dtype=complex); minus_b[0]=1/np.sqrt(2); minus_b[7]=-1/np.sqrt(2)
V9 = np.zeros((512,2), dtype=complex)
V9[:,0] = np.kron(np.kron(plus_b,plus_b),plus_b)
V9[:,1] = np.kron(np.kron(minus_b,minus_b),minus_b)
print("Shor: projection")

# ============================================================
# Experiment
# ============================================================
gamma_values = np.array([0.0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5])
results = {}

for code_name, V, n, pc in [('five_qubit',V5,5,pc5), ('steane',V7,7,None)]:
    has_corr = pc is not None
    print(f"\n=== {code_name} (n={n}) ===")
    for lname, logical_sv in [('zero', np.array([1,0], dtype=complex)),
                               ('one', np.array([0,1], dtype=complex))]:
        logical_rho = np.outer(logical_sv, logical_sv.conj())
        fid_unc, fid_cod, ret = [], [], []
        t0 = time.time()

        for gamma in gamma_values:
            # Uncoded
            E0s = np.array([[1,0],[0,np.sqrt(1-gamma)]], dtype=complex)
            E1s = np.array([[0,np.sqrt(gamma)],[0,0]], dtype=complex)
            rho_u = E0s @ logical_rho @ E0s.conj().T + E1s @ logical_rho @ E1s.conj().T
            fid_unc.append(round(fidelity(rho_u, logical_sv), 6))

            # Coded
            rho_coded = V @ logical_rho @ V.conj().T
            kraus = precompute_kraus_ad(gamma, n)
            rho_noisy = apply_kraus(rho_coded, kraus)

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
            'n': n, 'noise': 'amplitude_damping', 'has_correction': has_corr,
            'gamma_values': [round(float(g),4) for g in gamma_values],
            'fid_uncoded': fid_unc, 'fid_coded': fid_cod,
            'codespace_retention': ret,
        }
        be = next((round(float(g),4) for i,g in enumerate(gamma_values) if g>0 and fid_cod[i]<fid_unc[i]), None)
        print(f"  |{lname}>: break-even ~{be}, {elapsed:.1f}s")

# ============================================================
# Summary
# ============================================================
print("\n\n=== AMPLITUDE DAMPING: |zero> ===")
print(f"{'gamma':<8} {'[[5,1,3]]*':<14} {'Steane':<14} {'uncoded':<12}")
print("-" * 48)
for i,g in enumerate(gamma_values):
    print(f"{g:.3f}   {results['five_qubit_zero']['fid_coded'][i]:.4f}        "
          f"{results['steane_zero']['fid_coded'][i]:.4f}        "
          f"{results['five_qubit_zero']['fid_uncoded'][i]:.4f}")
print("* = with correction")

print("\n=== AMPLITUDE DAMPING: |one> ===")
print(f"{'gamma':<8} {'[[5,1,3]]*':<14} {'Steane':<14} {'uncoded':<12}")
print("-" * 48)
for i,g in enumerate(gamma_values):
    print(f"{g:.3f}   {results['five_qubit_one']['fid_coded'][i]:.4f}        "
          f"{results['steane_one']['fid_coded'][i]:.4f}        "
          f"{results['five_qubit_one']['fid_uncoded'][i]:.4f}")

print("\n=== STATE ASYMMETRY (|0> - |1> fidelity) ===")
for cn in ['five_qubit','steane']:
    asym = [results[f'{cn}_zero']['fid_coded'][i] - results[f'{cn}_one']['fid_coded'][i]
            for i in range(len(gamma_values))]
    mi = int(np.argmax([abs(a) for a in asym]))
    print(f"  {cn}: max |asym| = {abs(asym[mi]):.4f} at gamma={gamma_values[mi]:.3f}")

print("\n=== CODESPACE RETENTION (|zero>) ===")
for cn in ['five_qubit','steane']:
    r = results[f'{cn}_zero']['codespace_retention']
    print(f"  {cn}: g=0.05:{r[1]:.4f}, g=0.2:{r[4]:.4f}, g=0.5:{r[-1]:.4f}")

with open('results/sprint_016a_structured_noise_amplitude.json', 'w') as f:
    json.dump(results, f, indent=2)
print("\nResults saved.")
