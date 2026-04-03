"""GPU-accelerated eigensolver — drop-in replacement for scipy.sparse.linalg.eigsh.

Usage: Replace `from scipy.sparse.linalg import eigsh` with `from gpu_utils import eigsh`
That's it. Same API. Automatically uses GPU when beneficial, falls back to CPU.

    from gpu_utils import eigsh
    evals, evecs = eigsh(H, k=4, which='SA')
"""
import numpy as np
from scipy.sparse.linalg import eigsh as cpu_eigsh

_has_cupy = False
try:
    import cupy as cp
    from cupyx.scipy.sparse import csr_matrix as cp_csr
    from cupyx.scipy.sparse.linalg import eigsh as _cp_eigsh
    _has_cupy = True
except ImportError:
    pass

# Use GPU when matrix dimension exceeds this threshold
GPU_THRESHOLD = 50000

def eigsh(A, k=6, which='SA', return_eigenvectors=True, **kwargs):
    """Drop-in replacement for scipy.sparse.linalg.eigsh.
    Uses GPU (CuPy) when matrix dimension > 50k and CuPy is available.
    Falls back to CPU otherwise or on GPU error."""
    n = A.shape[0]

    if _has_cupy and n > GPU_THRESHOLD:
        try:
            A_gpu = cp_csr(A)
            if return_eigenvectors:
                evals_gpu, evecs_gpu = _cp_eigsh(A_gpu, k=k, which=which, **kwargs)
                return cp.asnumpy(evals_gpu), cp.asnumpy(evecs_gpu)
            else:
                evals_gpu = _cp_eigsh(A_gpu, k=k, which=which,
                                       return_eigenvectors=False, **kwargs)
                return cp.asnumpy(evals_gpu)
        except Exception:
            # Fall back to CPU on any GPU error
            pass

    return cpu_eigsh(A, k=k, which=which, return_eigenvectors=return_eigenvectors, **kwargs)
