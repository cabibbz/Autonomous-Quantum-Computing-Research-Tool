#!/usr/bin/env python3
"""Sprint 052c: Refine q=10 crossing + correct FSS corrections + refit scaling law.

Key insight: FSS correction depends on size pair used.
From q=2,3 exact values:
  - n=4,6 pairs: correction ~4.8% (0.2386→0.250, 0.3178→0.333)
  - n=6,8 pairs: correction ~2.5% (0.2439→0.250, 0.3251→0.333)
  - n=8,10 pairs: correction ~1.5% (0.2462→0.250)

q=7 used n=4,6 → should be 4.8%, not 2.5%!
"""
import numpy as np, json, time, warnings
warnings.filterwarnings('ignore')
from scipy.sparse.linalg import eigsh
from scipy.sparse import kron, csr_matrix
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d

def potts_hamiltonian(q, n, g, J=1.0):
    I = csr_matrix(np.eye(q))
    X = np.zeros((q, q), dtype=complex)
    for a in range(q):
        X[(a+1) % q, a] = 1.0
    Xphc = csr_matrix(X + X.conj().T)
    projectors = []
    for a in range(q):
        P = np.zeros((q, q)); P[a, a] = 1.0
        projectors.append(csr_matrix(P))
    dim = q**n
    H = csr_matrix((dim, dim), dtype=complex)
    for i in range(n - 1):
        for a in range(q):
            op = csr_matrix(np.eye(1))
            for j in range(n):
                if j == i:
                    op = kron(op, projectors[a])
                elif j == i + 1:
                    op = kron(op, projectors[a])
                else:
                    op = kron(op, I)
            H += -J * op
    for i in range(n):
        op = csr_matrix(np.eye(1))
        for j in range(n):
            if j == i:
                op = kron(op, Xphc)
            else:
                op = kron(op, I)
        H += -g * op
    return H

def energy_gap(q, n, g, J=1.0):
    H = potts_hamiltonian(q, n, g, J)
    vals = eigsh(H, k=min(4, H.shape[0]-1), which='SA', return_eigenvectors=False)
    vals = np.sort(vals.real)
    return vals[1] - vals[0], vals[0]

results = {}

# === Step 1: Refine q=10 crossing with denser scan in [0.62, 0.70] ===
print("=== q=10 refined scan near crossing ===", flush=True)
g_refine = np.linspace(0.62, 0.70, 9)  # 0.01 spacing

# n=4 is cheap
data_n4 = []
for g in g_refine:
    gap, E0 = energy_gap(10, 4, g)
    data_n4.append({'g': float(g), 'gap_x_n': float(gap * 4)})
    print(f"  n=4 g={g:.3f}: Δ·4={gap*4:.5f}", flush=True)

# n=6 targeted (3 points to bracket)
data_n6 = []
g_n6_targets = [0.63, 0.65, 0.67, 0.69]
t0 = time.time()
for g in g_n6_targets:
    gap, E0 = energy_gap(10, 6, g)
    data_n6.append({'g': float(g), 'gap_x_n': float(gap * 6)})
    print(f"  n=6 g={g:.3f}: Δ·6={gap*6:.5f} ({time.time()-t0:.1f}s)", flush=True)

# Find refined crossing
f_n4 = interp1d([d['g'] for d in data_n4], [d['gap_x_n'] for d in data_n4], kind='cubic')
crossings = []
for i in range(len(data_n6) - 1):
    g1, g2 = data_n6[i]['g'], data_n6[i+1]['g']
    diff1 = float(f_n4(g1)) - data_n6[i]['gap_x_n']
    diff2 = float(f_n4(g2)) - data_n6[i+1]['gap_x_n']
    print(f"  g=[{g1:.2f},{g2:.2f}]: Δ·4-Δ·6 = [{diff1:.5f}, {diff2:.5f}]", flush=True)
    if diff1 * diff2 < 0:
        g_cross = g1 - diff1 * (g2 - g1) / (diff2 - diff1)
        crossings.append(float(g_cross))
        print(f"  *** Refined crossing: g = {g_cross:.5f} ***", flush=True)

q10_raw = crossings[0] if crossings else 0.652
print(f"\nq=10 raw crossing: {q10_raw:.4f}\n", flush=True)

results['q10_refined'] = {'n4_data': data_n4, 'n6_data': data_n6, 'raw_crossing': q10_raw}

# === Step 2: Apply correct FSS corrections ===
print("=== FSS corrections by size pair ===")
# Calibration from known exact values
# q=2 exact: 0.250
# q=3 exact: 0.333
fss_corrections = {
    '4,6': {
        'q2': {'raw': 0.2386, 'exact': 0.250, 'corr': 0.250/0.2386 - 1},
        'q3': {'raw': 0.3178, 'exact': 0.333, 'corr': 0.333/0.3178 - 1},
    },
    '6,8': {
        'q2': {'raw': 0.2439, 'exact': 0.250, 'corr': 0.250/0.2439 - 1},
        'q3': {'raw': 0.3251, 'exact': 0.333, 'corr': 0.333/0.3251 - 1},
    },
}

for pair, data in fss_corrections.items():
    corrs = [v['corr'] for v in data.values()]
    mean_corr = np.mean(corrs)
    print(f"  n={pair}: corrections = {[f'{c:.4f}' for c in corrs]}, mean = {mean_corr:.4f}")
    fss_corrections[pair]['mean'] = mean_corr

corr_46 = fss_corrections['4,6']['mean']
corr_68 = fss_corrections['6,8']['mean']
print(f"\n  n=4,6 correction: {corr_46:.4f} ({corr_46*100:.1f}%)")
print(f"  n=6,8 correction: {corr_68:.4f} ({corr_68*100:.1f}%)")

# Apply correct corrections to all q values
print("\n=== Corrected g_c values ===")
gc_corrected = {
    2: {'raw': 0.250, 'pair': 'exact', 'gc': 0.250},
    3: {'raw': 0.333, 'pair': 'exact', 'gc': 0.333},
    4: {'raw': 0.3823, 'pair': '6,8', 'gc': 0.3823 * (1 + corr_68)},
    5: {'raw': 0.4303, 'pair': '6,8', 'gc': 0.4303 * (1 + corr_68)},
    7: {'raw': 0.5106, 'pair': '4,6', 'gc': 0.5106 * (1 + corr_46)},  # WAS 2.5%, NOW 4.8%
    10: {'raw': q10_raw, 'pair': '4,6', 'gc': q10_raw * (1 + corr_46)},
}

for q in sorted(gc_corrected):
    d = gc_corrected[q]
    print(f"  q={q}: raw={d['raw']:.4f}, pair={d['pair']}, g_c={d['gc']:.4f}")

# Old vs new q=7
old_q7 = 0.5106 * 1.025
new_q7 = gc_corrected[7]['gc']
print(f"\n  q=7 OLD (2.5% corr): {old_q7:.4f}")
print(f"  q=7 NEW (4.8% corr): {new_q7:.4f}")
print(f"  Shift: {(new_q7-old_q7)/old_q7*100:.1f}%")

results['gc_corrected'] = {str(q): v for q, v in gc_corrected.items()}

# === Step 3: Refit scaling law with 6 points ===
print("\n=== Refitting g_c(q) with 6 points (including q=10) ===\n")

q_data = np.array(sorted(gc_corrected.keys()), dtype=float)
gc_data = np.array([gc_corrected[int(q)]['gc'] for q in q_data])
# Uncertainties: ~0 for exact, ~3% for gap-corrected
gc_err = np.array([0.001, 0.001, 0.012, 0.013, 0.026, 0.035])  # larger err for 4,6 pairs

def power_offset(q, a, b, c):
    return a * (q - 1)**b + c

def logarithmic(q, a, b):
    return a * np.log(q) + b

def power_law(q, a, b):
    return a * q**b

models = {
    'power_offset': (power_offset, [0.2, 0.5, 0.0], 3),
    'logarithmic': (logarithmic, [0.25, 0.05], 2),
    'power_law': (power_law, [0.15, 0.65], 2),
}

fits = {}
for name, (func, p0, nparams) in models.items():
    try:
        popt, pcov = curve_fit(func, q_data, gc_data, p0=p0, sigma=gc_err, absolute_sigma=True, maxfev=10000)
        perr = np.sqrt(np.diag(pcov))
        gc_fit = func(q_data, *popt)
        residuals = gc_data - gc_fit
        chi2 = np.sum((residuals / gc_err)**2)
        dof = len(q_data) - nparams
        chi2_red = chi2 / dof if dof > 0 else float('inf')
        rms = np.sqrt(np.mean(residuals**2))
        gc_20 = func(20, *popt)
        gc_50 = func(50, *popt)

        fits[name] = {'params': popt.tolist(), 'param_err': perr.tolist(),
                      'chi2_red': float(chi2_red), 'rms': float(rms),
                      'residuals': residuals.tolist(),
                      'gc_20': float(gc_20), 'gc_50': float(gc_50)}

        print(f"{name}:")
        print(f"  params = {popt}")
        print(f"  χ²/dof = {chi2_red:.4f}  (RMS = {rms:.5f})")
        print(f"  residuals = {[f'{r:.5f}' for r in residuals]}")
        print(f"  g_c(20) = {gc_20:.4f}, g_c(50) = {gc_50:.4f}\n")
    except Exception as e:
        print(f"{name}: FAILED - {e}\n")

# === Step 4: Blind prediction for q=20 ===
ranked = sorted([(k, v) for k, v in fits.items() if 'chi2_red' in v],
                key=lambda x: x[1]['chi2_red'])
print("=== Rankings ===")
for i, (name, v) in enumerate(ranked):
    print(f"  {i+1}. {name}: χ²/dof = {v['chi2_red']:.4f}, g_c(20) = {v['gc_20']:.4f}")

best_name, best = ranked[0]
print(f"\nBest fit: {best_name}")
print(f"  g_c(20) = {best['gc_20']:.4f}")
print(f"  g_c(50) = {best['gc_50']:.4f}")

# Check if simple formula fits
# Try: g_c = a * ln(q) + b
log_fit = fits.get('logarithmic', {})
if log_fit and 'params' in log_fit:
    a, b = log_fit['params']
    print(f"\n=== Simple formula: g_c ≈ {a:.4f}·ln(q) + {b:.4f} ===")
    for q in [2, 3, 4, 5, 7, 10, 20, 50]:
        pred = a * np.log(q) + b
        actual = gc_corrected.get(q, {}).get('gc', None)
        if actual:
            print(f"  q={q}: predicted={pred:.4f}, actual={actual:.4f}, error={abs(pred-actual)/actual*100:.1f}%")
        else:
            print(f"  q={q}: predicted={pred:.4f}")

results['fits'] = fits
results['ranking'] = [name for name, _ in ranked]
results['best_fit'] = best_name

with open('results/sprint_052c_refit.json', 'w') as f:
    json.dump(results, f, indent=2)
print("\nSaved to results/sprint_052c_refit.json")
