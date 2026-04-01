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

**Setup:** Fine-grained sweep of all syndrome measures across the phase transition (p=0 to 1.0) for both 3-qubit and 9-qubit (2-level concatenated) repetition codes. Track syndrome entropy, MI(S:logical_outcome), P(correct), Holevo information, and syndrome hierarchy (inner vs outer) in the concatenated code.

**Result:** The syndrome and logical information have fundamentally different transition behaviors.

Three transition points in the 3-qubit code:
- H(S) = 75% of max at p ≈ 0.14 (syndrome nearly saturated)
- Holevo = 0.5 at p ≈ 0.24 (half logical info retained)
- P(correct) = 0.75 at p ≈ 0.34 (correction still mostly works)

| p     | H(S)/Hmax | Holevo | MI(S:L)/H(L) | P(correct) |
|-------|-----------|--------|---------------|------------|
| 0.05  | 0.408     | 0.957  | 0.313         | 0.993      |
| 0.10  | 0.635     | 0.862  | 0.253         | 0.972      |
| 0.20  | 0.880     | 0.594  | 0.157         | 0.896      |
| 0.30  | 0.975     | 0.305  | 0.077         | 0.784      |
| 0.40  | 0.998     | 0.084  | 0.021         | 0.648      |
| 0.50  | 1.000     | 0.000  | 0.000         | 0.500      |

**Key findings:**

1. **Syndrome entropy saturates BEFORE logical info dies** — syndrome reaches 90% of max at p≈0.10, but Holevo is still 0.86. There's a wide gap (p=0.15 to 0.50) where the syndrome is nearly maximally noisy but the code still works! Error correction doesn't need a clean syndrome — it works through collective majority vote even when individual syndrome bits are highly uncertain.

2. **Only Holevo shows a phase transition** — syndrome entropy saturates smoothly without any sharp transition. MI(S:L) rises and falls smoothly. Only the Holevo information (via concatenation) sharpens into a step function. The syndrome is a "thermometer" that measures noise intensity — it doesn't itself undergo a transition.

3. **Concatenation sharpens Holevo but NOT syndrome** — 9q Holevo/3q Holevo ratio reaches 3.0 near threshold (dramatic sharpening). But syndrome entropy just scales linearly (8.0 = 4 × 2.0). Concatenation amplifies the logical transition without changing the syndrome's behavior.

4. **Syndrome hierarchy** — at low noise, the outer syndrome carries only 1.6% of total syndrome entropy (inner correction handles everything). At threshold, outer carries 25% (= 2/8, the theoretical maximum for uniform syndromes). The outer syndrome "turns on" as inner correction fails.

5. **MI(S:L) peaks at p≈0.20** — the syndrome is most informative about its own success in the transition region, exactly where the outcome is most uncertain. At low and high noise, the outcome is predictable without the syndrome.

**Surprise:** The syndrome is most useful precisely when it's most noisy. At p=0.20, syndrome entropy is 88% of max (very noisy), but MI(S:L) is at its peak. The syndrome's value comes not from being clean but from being differentially noisy across outcomes — some syndromes still have high correction probability while others don't. This differential is what MI(S:L) captures.

---

## Sprint 018 Summary

**Three experiments, one unified picture:**

1. **18a (syndrome-error MI):** The syndrome is a lossy compression of the error. Efficiency drops from 99% at low noise to 45-67% at high noise. Syndrome predictive power (MI with correction outcome) peaks at ~35% and decreases with noise. The syndrome tells you WHICH error happened but poorly predicts WHETHER correction succeeds.

2. **18b ([[5,1,3]] Holevo):** Under depolarizing noise, the 3-qubit bit-flip code ALWAYS beats [[5,1,3]] for Holevo information. The general code's ability to correct all error types can't compensate for its 5-qubit attack surface. Encoding with [[5,1,3]] actively destroys information at p > 0.13.

3. **18c (syndrome at transition):** Syndrome and logical information have completely different transition behaviors. Syndrome entropy saturates smoothly; only Holevo shows a phase transition. Concatenation sharpens the logical transition but leaves the syndrome unchanged. The syndrome is most useful (highest MI with outcome) precisely when it's most noisy.

**Unified insight:** Error correction works through collective majority vote, not through clean syndrome measurements. The syndrome is a noisy, lossy compression of the error that becomes maximally uncertain long before error correction fails. This reveals a fundamental asymmetry: the syndrome's job is not to identify the exact error (it can't, at practical noise levels) but to guide a correction that works on average. The phase transition in logical information is invisible in the syndrome — it emerges only at the collective, decoded level. This is the information-theoretic essence of fault tolerance: the macro (logical) behavior transitions sharply while the micro (syndrome) behavior degrades gradually.

**Next sprint ideas:**
- Toric/surface code entanglement structure — topological codes
- Combined T1+T2 noise model — realistic hardware noise
- Real hardware QEC — simulator vs device
- Syndrome decoding strategies — does ML decoding change the transition?
- Quantum channel capacity — what's the theoretical limit?
