# Sprint 018 — Syndrome Information: The Measurement Side of Error Correction

**Date:** 2026-03-31
**Goal:** Investigate the information content of error syndromes. The syndrome is a classical measurement that QEC uses to diagnose errors — how much information does it carry about the error, and how does this change at the threshold phase transition?

**Hypothesis:**
- The syndrome should carry exactly enough MI about the error to enable correction (for correctable errors)
- At the threshold, syndrome-error MI should undergo its own transition — above threshold, syndrome becomes ambiguous
- The [[5,1,3]] code should show a sharper phase transition in Holevo information than the 3-qubit code (corrects more error types)

**Connection:** Sprint 017 (threshold as phase transition) + Sprint 014 (QEC basics) + Sprint 015 (code families) + Sprint 016 (structured noise)

---

## Experiment 18a: Syndrome-Error Mutual Information

**Setup:** For 3-qubit bit-flip code (under bit-flip noise) and [[5,1,3]] code (under depolarizing noise), enumerate all possible error patterns at each noise rate. Compute: syndrome entropy H(S), total error entropy H(E), residual ambiguity H(E|S), correction success probability, and MI between syndrome and correction outcome.

**Result:** The syndrome is a lossy compression of the error — and its efficiency degrades with noise.

| p     | 3q H(S)/H(E) | 5q H(S)/H(E) | 3q P(corr) | 5q P(corr) | 3q MI(S:L) | 5q MI(S:L) |
|-------|---------------|---------------|------------|------------|------------|------------|
| 0.01  | 99.0%         | 98.0%         | 0.9997     | 0.9993     | 0.0015     | 0.0028     |
| 0.05  | 95.0%         | 90.5%         | 0.9928     | 0.9851     | 0.0194     | 0.0316     |
| 0.10  | 90.2%         | 82.0%         | 0.9720     | 0.9470     | 0.0467     | 0.0664     |
| 0.20  | 81.2%         | 67.8%         | 0.8960     | 0.8339     | 0.0754     | 0.0850     |
| 0.30  | 73.7%         | 57.1%         | 0.7840     | 0.7117     | 0.0581     | 0.0553     |
| 0.50  | 66.7%         | 44.6%         | 0.5000     | 0.5432     | 0.0000     | 0.0036     |

**Key findings:**

1. **Syndrome efficiency = H(S)/H(E)** — at low noise (~99%), syndrome captures nearly all error information (almost all errors are weight-1 and uniquely identified). At high noise, efficiency drops to 67% (3-qubit) and 45% ([[5,1,3]]) as multi-qubit errors create ambiguity that the syndrome can't resolve.

2. **The efficiency gap is structural** — 3-qubit code has 2 syndrome bits for 3 qubits (max efficiency = 2/3 = 66.7%). [[5,1,3]] has 4 syndrome bits for 10 error-space bits (4^5 patterns → 10 bits). At high noise, the syndrome is fundamentally unable to capture most of the error information.

3. **MI(S:L) peaks at intermediate noise** (~p=0.15-0.20) — at low noise, almost everything succeeds so there's nothing to predict. At high noise, outcomes become random. The syndrome is most informative in the transition region.

4. **Syndrome predictive power decreases monotonically** — even at best (low p), the syndrome only predicts ~35-38% of correction outcome uncertainty. The syndrome tells you WHICH error happened, not WHETHER correction will succeed.

5. **[[5,1,3]] P(correct) > 0.5 at p=0.5** — degeneracy! Some weight-2 errors are effectively corrected because error × correction = stabilizer element. The 3-qubit code hits exactly 0.5 (its threshold).

**Surprise:** The syndrome is a surprisingly poor predictor of its own success. Even at low noise where correction almost always works, knowing the exact syndrome only resolves ~35% of the remaining uncertainty about whether it worked. This is because most syndromes always lead to correct outcomes (weight-1 errors) — the rare failures (weight-2+) are spread across many syndromes. The syndrome tells you what went wrong, but not how badly.

---

## Experiment 18b: [[5,1,3]] Holevo Information — Phase Transition for a General Code

**Setup:** Compare Holevo information (logical distinguishability of |0_L⟩ vs |1_L⟩) for uncoded, 3-qubit bit-flip code, and [[5,1,3]] code, all under depolarizing noise. Also compare fidelity for |0⟩ and |+⟩ logical states.

**Result:** The 3-qubit code dominates [[5,1,3]] at ALL depolarizing noise levels for Holevo information.

| p     | Uncoded | 3-qubit | [[5,1,3]] |
|-------|---------|---------|-----------|
| 0.01  | 0.942   | 0.998   | 0.992     |
| 0.05  | 0.789   | 0.968   | 0.888     |
| 0.10  | 0.647   | 0.902   | 0.701     |
| 0.15  | 0.531   | 0.816   | 0.513     |
| 0.20  | 0.434   | 0.720   | 0.351     |
| 0.30  | 0.278   | 0.519   | 0.134     |
| 0.50  | 0.082   | 0.174   | 0.005     |

**Key findings:**

1. **3-qubit code ALWAYS beats [[5,1,3]] under depolarizing** — extends Sprint 017b's fidelity finding to the information-theoretic domain. The advantage is even larger for Holevo: at p=0.15, 3-qubit retains 82% vs [[5,1,3]]'s 51%.

2. **[[5,1,3]] Holevo crosses below uncoded at p≈0.13** — encoding actively destroys logical information at moderate noise. The 5 qubits are 5 noise targets; the ability to correct all error types can't compensate for the increased attack surface.

3. **Massive state asymmetry in 3-qubit code:** |0⟩ fidelity excellent (0.987 at p=0.10) but |+⟩ fidelity WORSE than uncoded at all p (0.826 vs 0.933 at p=0.10). The code can't correct phase flips. Yet Holevo (which measures |0⟩ vs |1⟩ distinguishability) stays high because it only needs Z-basis discrimination.

4. **[[5,1,3]] is state-symmetric** — |0⟩ and |+⟩ fidelities are identical (both 0.947 at p=0.10). The code treats all errors equally. But this generality costs 5 qubits = 5 noise targets.

5. **[[5,1,3]] Holevo curve is steeper** — drops from 0.99 to 0.005 between p=0.01 and p=0.50 (200x drop). Uncoded: 0.94 to 0.08 (12x). The code amplifies the transition but around the wrong "threshold" for depolarizing noise.

**Surprise:** Holevo information reveals a fundamental tension: specialized codes (3-qubit) preserve more information under symmetric noise because they concentrate protection on one error type. General codes ([[5,1,3]]) spread protection evenly but pay an overhead tax. The "best" code depends entirely on what information you're trying to protect, not just the noise rate. This is the information-theoretic version of Sprint 017b's fidelity result — and the effect is even more dramatic.

---

## Experiment 18c: Syndrome MI at the Phase Transition

*Status: pending*
