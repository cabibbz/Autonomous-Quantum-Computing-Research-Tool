"""Sprint 107: Compile all spectral decomposition data (106 + 107a + 107a2)."""
import numpy as np
import json

# Combine all data sources
all_data = {}

# q=2 from 106b (5 sizes, n=6-14)
all_data[2] = {
    'sizes': [6, 8, 10, 12, 14],
    'gap_m': [1.0352761804100794, 0.7803612880645261, 0.6257378601609034, 0.5221047688801672, 0.447857904413242],
    'me_sq': [3.7320508075688696, 3.847759065022569, 3.902113032590309, 3.931851652578146, 3.949855824363628],
    'chi_F': [0.5803418012614825, 0.7898169490339784, 0.9965864547265996, 1.2019891779372882, 1.4066068257315543],
    'dom_frac': [1.0]*5,
}

# q=3 from 106b (3 sizes, n=6,8,10)
all_data[3] = {
    'sizes': [6, 8, 10],
    'gap_m': [0.7683639309279018, 0.5710777189280911, 0.4538702293534147],
    'me_sq': [13.098871247927713, 14.644577515883052, 15.91068267480527],
    'chi_F': [3.69784483288012, 5.6130179090598995, 7.72370041396113],
    'dom_frac': [1.0]*3,
}

# q=5: 106b (n=6,8,9) + 107a (n=7) + 107a2 (n=10)
all_data[5] = {
    'sizes': [6, 7, 8, 9, 10],
    'gap_m': [0.5029468370081389, 0.4206341237968916, 0.36048291833108337, 0.31472893277771874, 0.2788274041018628],
    'me_sq': [54.14345168045072, 60.85020042290008, 67.44776847715522, 73.9837924524345, 80.48518464188606],
    'chi_F': [35.67389469442717, 49.13092824418727, 64.87960811398077, 82.98904025410144, 103.52495104316695],
    'dom_frac': [1.0]*5,
}

# q=7: 106b (n=6,7) + 107a2 (n=8)
all_data[7] = {
    'sizes': [6, 7, 8],
    'gap_m': [0.36530875473450575, 0.2984979792894853, 0.25045452523845313],
    'me_sq': [132.54522446535265, 154.87254671591157, 178.02974273772583],
    'chi_F': [165.53608277556873, 248.30965639706005, 354.76830360963896],
    'dom_frac': [1.0]*3,
}

print("=" * 80)
print("COMBINED SCALING FITS: Sprint 106 + 107 data")
print("=" * 80)

summary = {}
for q in [2, 3, 5, 7]:
    d = all_data[q]
    Ns = np.array(d['sizes'], dtype=float)
    log_N = np.log(Ns)
    log_gap = np.log(np.array(d['gap_m']))
    log_me = np.log(np.array(d['me_sq']))
    log_chi = np.log(np.array(d['chi_F']))

    z_m = -np.polyfit(log_N, log_gap, 1)[0]
    beta_me = np.polyfit(log_N, log_me, 1)[0]
    alpha = np.polyfit(log_N, log_chi, 1)[0]
    alpha_pred = beta_me + 2 * z_m - 1
    alpha_known = 0.315 * q + 0.469

    # Residuals for alpha fit
    p = np.polyfit(log_N, log_chi, 1)
    resid = log_chi - np.polyval(p, log_N)
    rms_resid = np.sqrt(np.mean(resid**2))

    print(f"\nq={q} ({len(Ns)} sizes, n={d['sizes']}):")
    print(f"  z_m      = {z_m:.5f}")
    print(f"  beta_me  = {beta_me:.5f}")
    print(f"  alpha    = {alpha:.5f}  (fit RMS residual: {rms_resid:.5f})")
    print(f"  alpha(β+2z-1) = {alpha_pred:.5f}")
    print(f"  alpha(known)  = {alpha_known:.3f}")
    print(f"  alpha - alpha(known) = {alpha - alpha_known:.4f}")

    # Pairwise
    print(f"  Pairwise z_m:     ", end="")
    for i in range(len(Ns)-1):
        pw = -(log_gap[i+1] - log_gap[i]) / (log_N[i+1] - log_N[i])
        print(f"({d['sizes'][i]},{d['sizes'][i+1]})→{pw:.4f}  ", end="")
    print()

    print(f"  Pairwise beta_me: ", end="")
    for i in range(len(Ns)-1):
        pw = (log_me[i+1] - log_me[i]) / (log_N[i+1] - log_N[i])
        print(f"({d['sizes'][i]},{d['sizes'][i+1]})→{pw:.4f}  ", end="")
    print()

    print(f"  Pairwise alpha:   ", end="")
    for i in range(len(Ns)-1):
        pw = (log_chi[i+1] - log_chi[i]) / (log_N[i+1] - log_N[i])
        print(f"({d['sizes'][i]},{d['sizes'][i+1]})→{pw:.4f}  ", end="")
    print()

    summary[q] = {'z_m': z_m, 'beta_me': beta_me, 'alpha': alpha, 'alpha_pred': alpha_pred}

# Linear fits of z_m(q), beta_me(q), alpha(q)
qs = np.array([2, 3, 5, 7], dtype=float)
z_ms = np.array([summary[q]['z_m'] for q in [2, 3, 5, 7]])
betas = np.array([summary[q]['beta_me'] for q in [2, 3, 5, 7]])
alphas = np.array([summary[q]['alpha'] for q in [2, 3, 5, 7]])

p_z = np.polyfit(qs, z_ms, 1)
p_b = np.polyfit(qs, betas, 1)
p_a = np.polyfit(qs, alphas, 1)

print(f"\n{'=' * 80}")
print("LINEAR FITS vs q:")
print(f"  z_m(q)    = {p_z[0]:.4f}q + {p_z[1]:.4f}  (Sprint 106: 0.065q + 0.843)")
print(f"  beta_me(q)= {p_b[0]:.4f}q + {p_b[1]:.4f}  (Sprint 106: 0.182q - 0.201)")
print(f"  alpha(q)  = {p_a[0]:.4f}q + {p_a[1]:.4f}  (Sprint 103: 0.315q + 0.469)")

# Reconstructed alpha from component fits
print(f"\n  Reconstructed alpha(q) = beta + 2z - 1 = {p_b[0]+2*p_z[0]:.4f}q + {p_b[1]+2*p_z[1]-1:.4f}")
print(f"  Direct alpha(q)                        = {p_a[0]:.4f}q + {p_a[1]:.4f}")

# Cross-check table
print(f"\n{'=' * 80}")
print("CROSS-CHECK TABLE:")
print(f"{'q':>3} {'z_m':>8} {'beta_me':>8} {'alpha':>8} {'pred':>8} {'diff':>8}")
print("-" * 50)
for q in [2, 3, 5, 7]:
    s = summary[q]
    diff = s['alpha'] - s['alpha_pred']
    print(f"{q:>3} {s['z_m']:>8.4f} {s['beta_me']:>8.4f} {s['alpha']:>8.4f} {s['alpha_pred']:>8.4f} {diff:>8.5f}")

# Save compiled results
compiled = {
    'experiment': '107_compiled_scaling',
    'data': {f'q{q}': {
        'sizes': all_data[q]['sizes'],
        'z_m': float(summary[q]['z_m']),
        'beta_me': float(summary[q]['beta_me']),
        'alpha_fit': float(summary[q]['alpha']),
        'alpha_predicted': float(summary[q]['alpha_pred']),
    } for q in [2, 3, 5, 7]},
    'linear_fits': {
        'z_m': [float(p_z[0]), float(p_z[1])],
        'beta_me': [float(p_b[0]), float(p_b[1])],
        'alpha': [float(p_a[0]), float(p_a[1])],
    }
}
with open('results/sprint_107_compiled.json', 'w') as f:
    json.dump(compiled, f, indent=2)
print("\nCompiled results saved.")
