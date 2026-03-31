"""Sprint 016c: Noise fingerprint for QEC codes.

Tests [[5,1,3]] (with correction) under 4 noise channels:
depolarizing, amplitude damping, phase damping, bit-flip.
For each: measures fidelity for |0>, |1>, |+> logical states.
This gives a 4×3 = 12-dimensional noise fingerprint.

Also compares Steane [[7,1,3]] (projection-only) for the same channels.
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

# ============================================================
# Noise channels
# ============================================================

def depolarizing_channel(rho, p, n):
    for q in range(n):
        XI = pauli_on_qubit(X,q,n); YI = pauli_on_qubit(Y,q,n); ZI = pauli_on_qubit(Z,q,n)
        rho = (1-p)*rho + (p/3)*(XI@rho@XI + YI@rho@YI + ZI@rho@ZI)
    return rho

def amplitude_damping_channel(rho, gamma, n):
    for q in range(n):
        E0 = np.array([[1,0],[0,np.sqrt(1-gamma)]],dtype=complex)
        E1 = np.array([[0,np.sqrt(gamma)],[0,0]],dtype=complex)
        ops0=[I2]*n; ops0[q]=E0; ops1=[I2]*n; ops1[q]=E1
        E0f=kron_list(ops0); E1f=kron_list(ops1)
        rho = E0f@rho@E0f.conj().T + E1f@rho@E1f.conj().T
    return rho

def phase_damping_channel(rho, lam, n):
    for q in range(n):
        E0 = np.array([[1,0],[0,np.sqrt(1-lam)]],dtype=complex)
        E1 = np.array([[0,0],[0,np.sqrt(lam)]],dtype=complex)
        ops0=[I2]*n; ops0[q]=E0; ops1=[I2]*n; ops1[q]=E1
        E0f=kron_list(ops0); E1f=kron_list(ops1)
        rho = E0f@rho@E0f.conj().T + E1f@rho@E1f.conj().T
    return rho

def bit_flip_channel(rho, p, n):
    for q in range(n):
        XI = pauli_on_qubit(X,q,n)
        rho = (1-p)*rho + p*(XI@rho@XI)
    return rho

def phase_flip_channel(rho, p, n):
    for q in range(n):
        ZI = pauli_on_qubit(Z,q,n)
        rho = (1-p)*rho + p*(ZI@rho@ZI)
    return rho

# ============================================================
# Build codes
# ============================================================
print("=== Building codes ===")

stab5 = ['XZZXI','IXZZX','XIXZZ','ZXIXZ']

def encode_stab(stabilizers, n):
    dim = 2**n; state = np.zeros(dim, dtype=complex); state[0] = 1.0
    for s in stabilizers:
        S = pauli_string(s); state = (np.eye(dim)+S)/2.0 @ state
        nm = np.linalg.norm(state)
        if nm < 1e-12: raise ValueError(f"{s} annihilated!")
        state /= nm
    return state

sv5_0 = encode_stab(stab5, 5)
sv5_1 = pauli_string('XXXXX') @ sv5_0; sv5_1 /= np.linalg.norm(sv5_1)
V5 = np.zeros((32,2), dtype=complex); V5[:,0]=sv5_0; V5[:,1]=sv5_1

stab_ops5 = [pauli_string(s) for s in stab5]
table5 = {0: np.eye(32, dtype=complex)}
for q in range(5):
    for P in [X,Y,Z]:
        error = pauli_on_qubit(P, q, 5); s = 0
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

stab7 = ['IIIXXXX','IXXIIXX','XIXIXIX','IIIZZZZ','IZZIIZZ','ZIZIZIZ']
sv7_0 = encode_stab(stab7, 7)
sv7_1 = pauli_string('XXXXXXX') @ sv7_0; sv7_1 /= np.linalg.norm(sv7_1)
V7 = np.zeros((128,2), dtype=complex); V7[:,0]=sv7_0; V7[:,1]=sv7_1
print("Steane: projection")

# ============================================================
# Noise fingerprint at fixed noise strength
# ============================================================
# Use p=0.1 for all channels (moderate noise)
p_test = 0.1

noise_channels = [
    ('depolarizing', lambda rho,n: depolarizing_channel(rho, p_test, n)),
    ('amplitude_damping', lambda rho,n: amplitude_damping_channel(rho, p_test, n)),
    ('phase_damping', lambda rho,n: phase_damping_channel(rho, p_test, n)),
    ('bit_flip', lambda rho,n: bit_flip_channel(rho, p_test, n)),
    ('phase_flip', lambda rho,n: phase_flip_channel(rho, p_test, n)),
]

logical_states = [
    ('zero', np.array([1,0], dtype=complex)),
    ('one', np.array([0,1], dtype=complex)),
    ('plus', np.array([1,1], dtype=complex)/np.sqrt(2)),
]

results = {'p_test': p_test, 'codes': {}}

for code_name, V, n, pc in [('five_qubit',V5,5,pc5), ('steane',V7,7,None)]:
    has_corr = pc is not None
    print(f"\n=== {code_name} (n={n}) ===")
    code_results = {}

    for noise_name, noise_fn in noise_channels:
        t0 = time.time()
        noise_data = {}
        for lname, logical_sv in logical_states:
            logical_rho = np.outer(logical_sv, logical_sv.conj())

            # Uncoded
            rho_unc = noise_fn(logical_rho, 1)  # n=1 for single qubit
            f_unc = fidelity(rho_unc, logical_sv)

            # Coded
            rho_coded = V @ logical_rho @ V.conj().T
            rho_noisy = noise_fn(rho_coded, n)

            if has_corr:
                rho_c = sum(cp @ rho_noisy @ cp.conj().T for cp in pc)
                rho_dec = V.conj().T @ rho_c @ V
            else:
                rho_dec = V.conj().T @ rho_noisy @ V

            tr = np.trace(rho_dec).real
            retention = tr
            if tr > 1e-12: rho_dec /= tr
            f_cod = fidelity(rho_dec, logical_sv)

            noise_data[lname] = {
                'fid_uncoded': round(f_unc, 6),
                'fid_coded': round(f_cod, 6),
                'retention': round(retention, 6),
                'improvement': round(f_cod - f_unc, 6),
            }

        code_results[noise_name] = noise_data
        elapsed = time.time()-t0
        fids = [noise_data[s]['fid_coded'] for s in ['zero','one','plus']]
        print(f"  {noise_name}: |0>={fids[0]:.4f} |1>={fids[1]:.4f} |+>={fids[2]:.4f} [{elapsed:.1f}s]")

    results['codes'][code_name] = {
        'n': n, 'has_correction': has_corr, 'noise_data': code_results
    }

# ============================================================
# Also do a sweep at multiple noise levels for the most interesting comparison
# ============================================================
print("\n\n=== [[5,1,3]] fidelity sweep (|plus>) ===")
p_values = np.array([0.0, 0.05, 0.1, 0.2, 0.3])
sweep_results = {}

for noise_name, noise_fn_maker in [
    ('depolarizing', lambda p: lambda rho,n: depolarizing_channel(rho, p, n)),
    ('phase_damping', lambda p: lambda rho,n: phase_damping_channel(rho, p, n)),
    ('bit_flip', lambda p: lambda rho,n: bit_flip_channel(rho, p, n)),
    ('phase_flip', lambda p: lambda rho,n: phase_flip_channel(rho, p, n)),
]:
    fids = []
    logical_sv = np.array([1,1], dtype=complex)/np.sqrt(2)
    logical_rho = np.outer(logical_sv, logical_sv.conj())
    for p in p_values:
        noise_fn = noise_fn_maker(p)
        rho_coded = V5 @ logical_rho @ V5.conj().T
        rho_noisy = noise_fn(rho_coded, 5)
        rho_c = sum(cp @ rho_noisy @ cp.conj().T for cp in pc5)
        rho_dec = V5.conj().T @ rho_c @ V5
        tr = np.trace(rho_dec).real
        if tr > 1e-12: rho_dec /= tr
        fids.append(round(fidelity(rho_dec, logical_sv), 6))
    sweep_results[noise_name] = fids
    print(f"  {noise_name}: {fids}")

results['sweep'] = {
    'p_values': [round(float(p),4) for p in p_values],
    'code': 'five_qubit', 'state': 'plus',
    'fidelities': sweep_results,
}

# ============================================================
# Analysis
# ============================================================
print("\n\n" + "="*70)
print("=== NOISE FINGERPRINT TABLE (fidelity at p=0.1) ===")
print("="*70)
print(f"\n{'Code':<15} {'Noise':<20} {'|0>':<8} {'|1>':<8} {'|+>':<8} {'asym(0-1)':<10} {'asym(0-+)':<10}")
print("-" * 79)

for cn in ['five_qubit','steane']:
    for nn in ['depolarizing','amplitude_damping','phase_damping','bit_flip','phase_flip']:
        d = results['codes'][cn]['noise_data'][nn]
        f0 = d['zero']['fid_coded']; f1 = d['one']['fid_coded']; fp = d['plus']['fid_coded']
        print(f"{cn:<15} {nn:<20} {f0:.4f}  {f1:.4f}  {fp:.4f}  {f0-f1:+.4f}    {f0-fp:+.4f}")
    print()

# Key metric: which noise is worst for each code?
print("=== WORST NOISE PER CODE (lowest avg fidelity) ===")
for cn in ['five_qubit','steane']:
    worst_noise = None; worst_avg = 1.0
    for nn in ['depolarizing','amplitude_damping','phase_damping','bit_flip','phase_flip']:
        d = results['codes'][cn]['noise_data'][nn]
        avg = (d['zero']['fid_coded'] + d['one']['fid_coded'] + d['plus']['fid_coded']) / 3
        if avg < worst_avg:
            worst_avg = avg; worst_noise = nn
    print(f"  {cn}: worst = {worst_noise} (avg fid = {worst_avg:.4f})")

# Symmetry metric: how symmetric is each code?
print("\n=== CODE SYMMETRY (state-independence) ===")
for cn in ['five_qubit','steane']:
    total_asym = 0
    for nn in ['depolarizing','amplitude_damping','phase_damping','bit_flip','phase_flip']:
        d = results['codes'][cn]['noise_data'][nn]
        fids = [d[s]['fid_coded'] for s in ['zero','one','plus']]
        total_asym += max(fids) - min(fids)
    print(f"  {cn}: total asymmetry = {total_asym:.4f} (lower = more symmetric)")

with open('results/sprint_016c_noise_fingerprint_codes.json', 'w') as f:
    json.dump(results, f, indent=2)
print("\nResults saved.")
