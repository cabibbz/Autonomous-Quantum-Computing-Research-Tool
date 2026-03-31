# Sprint 002 — CHSH Violation Landscape & Noise Sensitivity

**Date:** 2026-03-31
**Status:** In progress

## Idea

Sprint 001 showed S=2.834 at the optimal CHSH angles. Two questions follow:
1. What does the full angle landscape look like? Where exactly does violation occur?
2. How fragile is the violation? How much noise kills it?

Both experiments use only 2 qubits — well within resource limits.

## Experiments

### Experiment 2a: CHSH Angle Sweep

Computed exact CHSH S over a 50x50 grid of (a1, b0) angles with a0=0, b1=b0+pi/2.

Used the analytic formula E(a,b)=cos(a-b) for Phi+ state — verified against statevector at 3 test points (all match to 10 decimal places).

**Results:**
- S ranges from **-2.825 to +2.828** across the landscape
- Maximum at a1=1.539 (pi/2), b0=0.769 (pi/4) — the known optimal
- **Only 19.9% of the angle space violates the classical bound |S|>2**
- The violation region forms distinct "islands" in angle space
- At a1=pi, b0=pi/2: S=2.000 exactly — right at the classical boundary
- The landscape has 4-fold symmetry (S and -S regions are mirror images)

**Insight:** Quantum advantage in CHSH is not generic — it requires specific angle choices. ~80% of measurement settings don't violate the classical bound at all. This makes the violation more remarkable: nature supports correlations stronger than classical physics allows, but only if you know where to look.

### Experiment 2b: CHSH Under Depolarizing Noise

Swept depolarizing noise rate p from 0% to 30% on all gates (1q and 2q), measuring CHSH S at optimal angles with 20k shots per setting.

**Results:**
| Noise rate | S | Violates? |
|-----------|------|-----------|
| 0.00 | 2.826 | Yes |
| 0.05 | 2.374 | Yes |
| 0.09 | 2.040 | Barely |
| **0.095** | **2.000** | **Critical point** |
| 0.10 | 1.963 | No |
| 0.15 | 1.619 | No |
| 0.20 | 1.303 | No |
| 0.30 | 0.811 | No |

**Critical noise rate: ~9.5% depolarizing error per gate.**

At this point, quantum correlations become indistinguishable from classical ones.

**Insights:**
- Quantum advantage in CHSH is surprisingly fragile — only ~10% noise kills it
- The degradation is roughly exponential: S decays as ~(1-p)^k where k relates to circuit depth
- At p=0.30, S=0.81 — well into classical territory, approaching the uncorrelated limit of 0
- This has implications for real hardware: IBM's current error rates (~0.1-1% per gate) should allow CHSH violation, but more complex protocols may not survive
- The critical threshold is a hard physical boundary — no error mitigation strategy changes when S *actually* violates the bound, only when we can *detect* the violation through noise
