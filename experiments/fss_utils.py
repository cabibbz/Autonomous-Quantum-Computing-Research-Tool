"""Finite-size scaling utilities with proper error bars.

Usage:
    from fss_utils import fit_power_law, pairwise_exponents

    # Power-law fit: y = A * N^alpha
    result = fit_power_law(sizes, chi_values)
    print(f"alpha = {result['alpha']:.4f} ± {result['alpha_err']:.4f}")

    # Pairwise slopes between consecutive sizes
    pairs = pairwise_exponents(sizes, chi_values)
    for p in pairs:
        print(f"({p['n1']},{p['n2']}): alpha = {p['alpha']:.4f}")

    # Decomposition: chi_F spectral exponents
    decomp = fit_spectral_exponents(sizes, chi_F, gaps, me_sq)
    print(f"alpha={decomp['alpha']:.3f}±{decomp['alpha_err']:.3f}, "
          f"z_m={decomp['z_m']:.3f}±{decomp['z_m_err']:.3f}")
"""
import numpy as np
from lmfit import Model


def _power_law(x, A, alpha):
    return A * np.power(x, alpha)


def fit_power_law(sizes, values, errors=None):
    """Fit y = A * N^alpha. Returns dict with alpha, alpha_err, A, A_err, r_squared."""
    sizes = np.asarray(sizes, dtype=float)
    values = np.asarray(values, dtype=float)
    # Initial guess from log-log polyfit
    p = np.polyfit(np.log(sizes), np.log(np.abs(values)), 1)
    model = Model(_power_law)
    params = model.make_params(A=np.exp(p[1]), alpha=p[0])
    weights = 1.0 / np.asarray(errors) if errors is not None else None
    result = model.fit(values, params, x=sizes, weights=weights)
    ss_res = np.sum(result.residual**2)
    ss_tot = np.sum((values - values.mean())**2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
    return {
        'alpha': result.params['alpha'].value,
        'alpha_err': result.params['alpha'].stderr or 0.0,
        'A': result.params['A'].value,
        'A_err': result.params['A'].stderr or 0.0,
        'r_squared': r2,
    }


def pairwise_exponents(sizes, values):
    """Pairwise log-log slopes between consecutive sizes.
    Returns list of dicts with n1, n2, alpha."""
    sizes = np.asarray(sizes, dtype=float)
    values = np.asarray(values, dtype=float)
    log_N = np.log(sizes)
    log_v = np.log(np.abs(values))
    pairs = []
    for i in range(len(sizes) - 1):
        dln = log_N[i + 1] - log_N[i]
        if abs(dln) < 1e-12:
            continue
        alpha = (log_v[i + 1] - log_v[i]) / dln
        pairs.append({'n1': int(sizes[i]), 'n2': int(sizes[i + 1]), 'alpha': alpha})
    return pairs


def fit_spectral_exponents(sizes, chi_F, gaps, me_sq):
    """Fit chi_F spectral decomposition exponents: alpha, z_m, beta_me.

    alpha = d(log chi_F)/d(log N)   — chi_F ~ N^alpha
    z_m   = -d(log gap)/d(log N)    — gap ~ N^{-z_m}
    beta_me = d(log |me|²)/d(log N) — |me|² ~ N^{beta_me}

    Checks: alpha ≈ beta_me + 2*z_m - 1 (exact decomposition).
    Returns dict with all exponents, errors, and reconstruction check.
    """
    r_chi = fit_power_law(sizes, chi_F)
    r_gap = fit_power_law(sizes, gaps)
    r_me = fit_power_law(sizes, me_sq)
    alpha = r_chi['alpha']
    z_m = -r_gap['alpha']  # gap decreases, so negate
    beta_me = r_me['alpha']
    recon = beta_me + 2 * z_m - 1
    return {
        'alpha': alpha, 'alpha_err': r_chi['alpha_err'],
        'z_m': z_m, 'z_m_err': r_gap['alpha_err'],
        'beta_me': beta_me, 'beta_me_err': r_me['alpha_err'],
        'reconstruction': recon,
        'recon_error': abs(recon - alpha),
        'nu_eff': 1.0 / (0.5 * (alpha + 1)) if alpha > -1 else None,
        'r_squared': {'chi_F': r_chi['r_squared'],
                      'gap': r_gap['r_squared'],
                      'me_sq': r_me['r_squared']},
    }
