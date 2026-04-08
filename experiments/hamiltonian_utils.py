"""Shared Hamiltonian builders for Potts/clock models.

Usage:
    from hamiltonian_utils import build_hybrid_parts, build_sq_potts_parts

    # Get coupling + field separately (for chi_F spectral decomposition):
    H_coup, H_field = build_hybrid_parts(n, q)
    H = H_coup + g * H_field

    # Or get full Hamiltonian directly:
    from hamiltonian_utils import build_hybrid_H, build_sq_potts_H
    H = build_hybrid_H(n, q, g)
    H = build_sq_potts_H(n, q)       # defaults to g=1/q (self-dual)

Two models:
    hybrid:  H = -Σ δ(s_i,s_j) - g(X + X†)           Z_q symmetry
    sq:      H = -Σ δ(s_i,s_j) - g Σ_{k=1}^{q-1} X^k  S_q symmetry
"""
import numpy as np
from scipy.sparse import coo_matrix, csr_matrix


def _decode_states(n, q):
    """Decode all q^n basis states into digit arrays."""
    dim = q**n
    all_idx = np.arange(dim, dtype=np.int64)
    digits = np.zeros((dim, n), dtype=np.int64)
    tmp = all_idx.copy()
    for site in range(n):
        digits[:, site] = tmp % q
        tmp //= q
    powers = q ** np.arange(n, dtype=np.int64)
    return dim, all_idx, digits, powers


def _build_coupling(all_idx, digits, n, dim):
    """Potts delta coupling: -Σ δ(s_i, s_{i+1}), periodic BC."""
    diag_vals = np.zeros(dim, dtype=np.float64)
    for site in range(n):
        nxt = (site + 1) % n
        diag_vals -= (digits[:, site] == digits[:, nxt]).astype(np.float64)
    return csr_matrix((diag_vals, (all_idx, all_idx)), shape=(dim, dim))


def _build_field(all_idx, digits, powers, n, q, dim, shifts):
    """Off-diagonal field: -Σ_site Σ_{k in shifts} |shifted><original|."""
    rows_list, cols_list, vals_list = [], [], []
    for site in range(n):
        pw = powers[site]
        old_digit = digits[:, site]
        for k in shifts:
            new_digit = (old_digit + k) % q
            delta = (new_digit.astype(np.int64) - old_digit.astype(np.int64)) * pw
            rows_list.append(all_idx + delta)
            cols_list.append(all_idx)
            vals_list.append(np.full(dim, -1.0, dtype=np.float64))
    rows = np.concatenate(rows_list)
    cols = np.concatenate(cols_list)
    vals = np.concatenate(vals_list)
    return coo_matrix((vals, (rows, cols)), shape=(dim, dim)).tocsr()


def build_hybrid_parts(n, q):
    """Hybrid Potts-clock: delta coupling + (X + X†) clock field.
    Returns (H_coup, H_field). Use: H = H_coup + g * H_field."""
    dim, all_idx, digits, powers = _decode_states(n, q)
    H_coup = _build_coupling(all_idx, digits, n, dim)
    H_field = _build_field(all_idx, digits, powers, n, q, dim, [1, q - 1])
    return H_coup, H_field


def build_sq_potts_parts(n, q):
    """Standard S_q Potts: delta coupling + Σ_{k=1}^{q-1} X^k field.
    Returns (H_coup, H_field). Use: H = H_coup + g * H_field."""
    dim, all_idx, digits, powers = _decode_states(n, q)
    H_coup = _build_coupling(all_idx, digits, n, dim)
    H_field = _build_field(all_idx, digits, powers, n, q, dim, range(1, q))
    return H_coup, H_field


def build_hybrid_H(n, q, g):
    """Full hybrid Hamiltonian at coupling g."""
    H_coup, H_field = build_hybrid_parts(n, q)
    return H_coup + g * H_field


def build_sq_potts_H(n, q, g=None):
    """Full S_q Potts Hamiltonian. g defaults to 1/q (self-dual point)."""
    if g is None:
        g = 1.0 / q
    H_coup, H_field = build_sq_potts_parts(n, q)
    return H_coup + g * H_field
