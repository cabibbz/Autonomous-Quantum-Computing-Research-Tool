"""Sprint 121b: Log correction fits at q=4 for both hybrid and S_q models.

Fit chi_F = A * N^alpha * (ln N)^{-p} to extract log correction exponent p.

S_q q=4 literature (Salas-Sokal 1997, Balog et al. 2007): marginal operator
gives multiplicative log corrections. Asymptotic alpha=2.0, but finite-size
alpha~1.77. Expect p > 0 (log suppression).

Hybrid q=4 (Sprint 120): alpha drifts downward (1.64->1.49). Different
universality class. Compare p values.

Also extends S_q q=4 analysis to match hybrid sizes (n=4-11).

Data sources:
  S_q q=4: sprint_118a_q4_chif.json (n=4-11)
  Hybrid q=4: sprint_120b_hybrid_chif_q4.json (n=4-11)
"""
import numpy as np
import json, time, os
from lmfit import Model
from db_utils import record

results = {
    'experiment': '121b_log_corrections_q4',
    'sprint': 121,
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'data': {},
}

def save():
    outpath = os.path.join(os.path.dirname(__file__), '..', 'results', 'sprint_121b_log_corrections_q4.json')
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    with open(outpath, 'w') as f:
        json.dump(results, f, indent=2, default=str)

# ---- LOAD DATA ----
res_dir = os.path.join(os.path.dirname(__file__), '..', 'results')

with open(os.path.join(res_dir, 'sprint_118a_q4_chif.json')) as f:
    sq_data = json.load(f)

with open(os.path.join(res_dir, 'sprint_120b_hybrid_chif_q4.json')) as f:
    hyb_data = json.load(f)

# Extract S_q q=4 data (n=4-11)
sq_sizes = []
sq_chif = []
sq_gaps = []
sq_me_sq = []
for key in sorted(sq_data['data'].keys()):
    d = sq_data['data'][key]
    sq_sizes.append(d['n'])
    sq_chif.append(d['chi_F'])
    sq_gaps.append(d['gap_multiplet'])
    sq_me_sq.append(d['me_sq'])
sq_sizes = np.array(sq_sizes, dtype=float)
sq_chif = np.array(sq_chif)
sq_gaps = np.array(sq_gaps)
sq_me_sq = np.array(sq_me_sq)

# Extract hybrid q=4 data (n=4-11)
hyb_sizes = []
hyb_chif = []
hyb_gaps = []
hyb_me_sq = []
for key in sorted(hyb_data['data'].keys()):
    d = hyb_data['data'][key]
    if d.get('q', 4) == 4:
        hyb_sizes.append(d['n'])
        hyb_chif.append(d['chi_F'])
        hyb_gaps.append(d['gap_multiplet'])
        hyb_me_sq.append(d['me_sq'])
hyb_sizes = np.array(hyb_sizes, dtype=float)
hyb_chif = np.array(hyb_chif)
hyb_gaps = np.array(hyb_gaps)
hyb_me_sq = np.array(hyb_me_sq)

print("=" * 70)
print("EXPERIMENT 121b: LOG CORRECTION FITS AT q=4")
print("=" * 70)

print(f"\nS_q q=4: {len(sq_sizes)} sizes (n={int(sq_sizes[0])}-{int(sq_sizes[-1])})")
print(f"Hybrid q=4: {len(hyb_sizes)} sizes (n={int(hyb_sizes[0])}-{int(hyb_sizes[-1])})")

# ---- FIT MODELS ----
def pure_power(x, A, alpha):
    """chi_F = A * N^alpha"""
    return A * np.power(x, alpha)

def log_corrected(x, A, alpha, p):
    """chi_F = A * N^alpha * (ln N)^{-p}"""
    return A * np.power(x, alpha) * np.power(np.log(x), -p)

def log_corrected_fixed(x, A, p):
    """chi_F = A * N^2 * (ln N)^{-p}  (fixed alpha=2 for S_q q=4)"""
    return A * np.power(x, 2.0) * np.power(np.log(x), -p)

# ---- S_q q=4 FITS ----
print("\n" + "=" * 40)
print("S_q POTTS q=4")
print("=" * 40)

# Pure power law
m1 = Model(pure_power)
p1 = m1.make_params(A=1.0, alpha=1.8)
r1 = m1.fit(sq_chif, x=sq_sizes, params=p1)
ss_res = np.sum(r1.residual**2)
ss_tot = np.sum((sq_chif - sq_chif.mean())**2)
r2_1 = 1 - ss_res/ss_tot
print(f"\n1. Pure power: chi_F = {r1.params['A'].value:.4f} * N^{r1.params['alpha'].value:.4f}")
print(f"   R^2 = {r2_1:.7f}, AIC = {r1.aic:.2f}")

# Log-corrected (free alpha)
m2 = Model(log_corrected)
p2 = m2.make_params(A=1.0, alpha=2.0, p=0.5)
p2['p'].min = -5
p2['p'].max = 10
r2 = m2.fit(sq_chif, x=sq_sizes, params=p2)
ss_res = np.sum(r2.residual**2)
r2_2 = 1 - ss_res/ss_tot
print(f"\n2. Log-corrected: chi_F = {r2.params['A'].value:.4f} * N^{r2.params['alpha'].value:.4f} * (ln N)^{{-{r2.params['p'].value:.4f}}}")
print(f"   R^2 = {r2_2:.7f}, AIC = {r2.aic:.2f}")
print(f"   alpha_err = {r2.params['alpha'].stderr or 0:.4f}, p_err = {r2.params['p'].stderr or 0:.4f}")

# Log-corrected (fixed alpha=2, S_q q=4 theoretical)
m3 = Model(log_corrected_fixed)
p3 = m3.make_params(A=0.5, p=0.5)
p3['p'].min = -5
p3['p'].max = 10
r3 = m3.fit(sq_chif, x=sq_sizes, params=p3)
ss_res = np.sum(r3.residual**2)
r2_3 = 1 - ss_res/ss_tot
print(f"\n3. Fixed alpha=2: chi_F = {r3.params['A'].value:.4f} * N^2 * (ln N)^{{-{r3.params['p'].value:.4f}}}")
print(f"   R^2 = {r2_3:.7f}, AIC = {r3.aic:.2f}")
print(f"   p_err = {r3.params['p'].stderr or 0:.4f}")

# Pairwise effective alpha
print("\n   Pairwise alpha drift:")
log_n = np.log(sq_sizes)
log_c = np.log(sq_chif)
for i in range(len(sq_sizes) - 1):
    alpha_pair = (log_c[i+1] - log_c[i]) / (log_n[i+1] - log_n[i])
    print(f"   ({int(sq_sizes[i])},{int(sq_sizes[i+1])}): alpha = {alpha_pair:.4f}")

results['data']['sq_q4'] = {
    'sizes': sq_sizes.tolist(),
    'chi_F': sq_chif.tolist(),
    'pure_power': {'A': r1.params['A'].value, 'alpha': r1.params['alpha'].value,
                   'r2': r2_1, 'aic': r1.aic},
    'log_corrected': {'A': r2.params['A'].value, 'alpha': r2.params['alpha'].value,
                      'alpha_err': r2.params['alpha'].stderr or 0,
                      'p': r2.params['p'].value, 'p_err': r2.params['p'].stderr or 0,
                      'r2': r2_2, 'aic': r2.aic},
    'log_fixed_alpha2': {'A': r3.params['A'].value, 'p': r3.params['p'].value,
                         'p_err': r3.params['p'].stderr or 0,
                         'r2': r2_3, 'aic': r3.aic},
}

# ---- HYBRID q=4 FITS ----
print("\n" + "=" * 40)
print("HYBRID POTTS-CLOCK q=4")
print("=" * 40)

# Pure power law
r1h = m1.fit(hyb_chif, x=hyb_sizes, params=m1.make_params(A=1.0, alpha=1.5))
ss_res = np.sum(r1h.residual**2)
ss_tot_h = np.sum((hyb_chif - hyb_chif.mean())**2)
r2_1h = 1 - ss_res/ss_tot_h
print(f"\n1. Pure power: chi_F = {r1h.params['A'].value:.4f} * N^{r1h.params['alpha'].value:.4f}")
print(f"   R^2 = {r2_1h:.7f}, AIC = {r1h.aic:.2f}")

# Log-corrected (free alpha)
p2h = m2.make_params(A=0.5, alpha=1.5, p=0.5)
p2h['p'].min = -5
p2h['p'].max = 10
r2h = m2.fit(hyb_chif, x=hyb_sizes, params=p2h)
ss_res = np.sum(r2h.residual**2)
r2_2h = 1 - ss_res/ss_tot_h
print(f"\n2. Log-corrected: chi_F = {r2h.params['A'].value:.4f} * N^{r2h.params['alpha'].value:.4f} * (ln N)^{{-{r2h.params['p'].value:.4f}}}")
print(f"   R^2 = {r2_2h:.7f}, AIC = {r2h.aic:.2f}")
print(f"   alpha_err = {r2h.params['alpha'].stderr or 0:.4f}, p_err = {r2h.params['p'].stderr or 0:.4f}")

# Pairwise effective alpha
print("\n   Pairwise alpha drift:")
log_nh = np.log(hyb_sizes)
log_ch = np.log(hyb_chif)
for i in range(len(hyb_sizes) - 1):
    alpha_pair = (log_ch[i+1] - log_ch[i]) / (log_nh[i+1] - log_nh[i])
    print(f"   ({int(hyb_sizes[i])},{int(hyb_sizes[i+1])}): alpha = {alpha_pair:.4f}")

results['data']['hybrid_q4'] = {
    'sizes': hyb_sizes.tolist(),
    'chi_F': hyb_chif.tolist(),
    'pure_power': {'A': r1h.params['A'].value, 'alpha': r1h.params['alpha'].value,
                   'r2': r2_1h, 'aic': r1h.aic},
    'log_corrected': {'A': r2h.params['A'].value, 'alpha': r2h.params['alpha'].value,
                      'alpha_err': r2h.params['alpha'].stderr or 0,
                      'p': r2h.params['p'].value, 'p_err': r2h.params['p'].stderr or 0,
                      'r2': r2_2h, 'aic': r2h.aic},
}

# ---- COMPARISON ----
print("\n" + "=" * 40)
print("MODEL COMPARISON AT q=4")
print("=" * 40)

sq_alpha_eff = r1.params['alpha'].value
sq_alpha_log = r2.params['alpha'].value
sq_p = r2.params['p'].value
hyb_alpha_eff = r1h.params['alpha'].value
hyb_alpha_log = r2h.params['alpha'].value
hyb_p = r2h.params['p'].value

print(f"\n{'':>20} {'S_q Potts':>14} {'Hybrid':>14}")
print("-" * 50)
print(f"{'Pure power alpha':>20} {sq_alpha_eff:>14.4f} {hyb_alpha_eff:>14.4f}")
print(f"{'Log-corr alpha':>20} {sq_alpha_log:>14.4f} {hyb_alpha_log:>14.4f}")
print(f"{'Log exponent p':>20} {sq_p:>14.4f} {hyb_p:>14.4f}")
print(f"{'Pure power R^2':>20} {r2_1:>14.7f} {r2_1h:>14.7f}")
print(f"{'Log-corr R^2':>20} {r2_2:>14.7f} {r2_2h:>14.7f}")
dAIC_sq = r1.aic - r2.aic
dAIC_hyb = r1h.aic - r2h.aic
print(f"{'dAIC (power-log)':>20} {dAIC_sq:>14.2f} {dAIC_hyb:>14.2f}")
print(f"\nS_q: log corrections {'SIGNIFICANT' if dAIC_sq > 2 else 'marginal'} (dAIC={dAIC_sq:.1f})")
print(f"Hybrid: log corrections {'SIGNIFICANT' if dAIC_hyb > 2 else 'marginal'} (dAIC={dAIC_hyb:.1f})")

if sq_alpha_log > 1.8:
    print(f"\nS_q log-corrected alpha={sq_alpha_log:.3f} -> consistent with asymptotic 2.0")
if hyb_alpha_log < 1.6:
    print(f"Hybrid log-corrected alpha={hyb_alpha_log:.3f} -> NOT approaching 2.0, different class")

results['data']['comparison'] = {
    'sq_pure_alpha': sq_alpha_eff, 'sq_log_alpha': sq_alpha_log, 'sq_p': sq_p,
    'hyb_pure_alpha': hyb_alpha_eff, 'hyb_log_alpha': hyb_alpha_log, 'hyb_p': hyb_p,
    'dAIC_sq': dAIC_sq, 'dAIC_hyb': dAIC_hyb,
}

# ---- RECORD ----
record(sprint=121, model='sq', q=4, n=0, quantity='alpha_log_corrected',
       value=sq_alpha_log, error=r2.params['alpha'].stderr or 0,
       method='log_correction_fit', notes=f'p={sq_p:.3f}')
record(sprint=121, model='sq', q=4, n=0, quantity='log_exponent_p',
       value=sq_p, error=r2.params['p'].stderr or 0,
       method='log_correction_fit')
record(sprint=121, model='hybrid', q=4, n=0, quantity='alpha_log_corrected',
       value=hyb_alpha_log, error=r2h.params['alpha'].stderr or 0,
       method='log_correction_fit', notes=f'p={hyb_p:.3f}')
record(sprint=121, model='hybrid', q=4, n=0, quantity='log_exponent_p',
       value=hyb_p, error=r2h.params['p'].stderr or 0,
       method='log_correction_fit')

save()
print("\nResults saved to results/sprint_121b_log_corrections_q4.json")
print("\nDone.")
