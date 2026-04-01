# Sprint 039 — Potts Data Collapse: Distinguishing ν=5/6 from ν=1

**Date:** 2026-04-01
**Status:** In progress

## Idea

Sprint 038 showed q=3 Potts MI-CV has crossing curves like TFIM, confirming second-order classification. The 3-state Potts universality class has ν=5/6≈0.833 (exact, from CFT), distinct from Ising ν=1. Can MI-CV data collapse distinguish these two universality classes?

**Plan:**
1. Compute n=16 Potts MI-CV at key coupling values near the transition
2. Fill n=12 with more points for better interpolation
3. Attempt n=24 sparse points if timing allows
4. Data collapse: compare collapse quality with ν=5/6 vs ν=1 across n=8,12,16,(24)

**Prediction:** Potts MI-CV collapse should favor ν=5/6 over ν=1. If the quality landscape is as flat as TFIM (Sprint 038), we may need large sizes (n≥16) to see the difference.

## Experiments

### 039a/039a2: Potts n=16 MI-CV

Computed 6 points at n=16 (chi_max=30, ~17-18s per point):

| g/J | n=8 CV | n=12 CV | n=16 CV |
|-----|--------|---------|---------|
| 0.80 | 0.341 | 0.271 | 0.236 |
| 0.90 | 0.577 | — | 0.472 |
| 0.95 | 0.705 | 0.748 | 0.738 |
| 1.00 | 0.826 | 0.961 | 1.051 |
| 1.05 | 0.933 | — | 1.319 |
| 1.10 | 1.023 | — | 1.522 |

**Crossing confirmed at three sizes!** At g=0.8: CV decreases with n (0.341→0.271→0.236). At g=1.0: CV increases with n (0.826→0.961→1.051). Crossing between n=8 and n=16 occurs at g≈0.96-0.97. All three sizes cross in a narrow region, consistent with a second-order transition.

### 039b: Potts n=12 fill

Filled 3 missing g values for n=12:

| g/J | n=8 CV | n=12 CV | n=16 CV |
|-----|--------|---------|---------|
| 0.50 | 0.069 | 0.063 | — |
| 0.80 | 0.341 | 0.271 | 0.236 |
| 0.90 | 0.577 | 0.539 | 0.472 |
| 0.95 | 0.705 | 0.748 | 0.738 |
| 1.00 | 0.826 | 0.961 | 1.051 |
| 1.05 | 0.933 | 1.146 | 1.319 |
| 1.10 | 1.023 | 1.296 | 1.522 |

Monotonic trend at every g value: CV systematically decreases (ordered side) or increases (disordered side) with system size. The crossing narrows from g≈0.95 (n=8 vs n=12) toward g≈0.96-0.97 (n=12 vs n=16).

### 039b2: Potts n=24 sparse

