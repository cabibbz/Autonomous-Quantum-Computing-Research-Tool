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

### Experiment 2b: CHSH Under Noise
(pending)
