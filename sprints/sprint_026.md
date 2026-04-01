# Sprint 026 — Active Syndrome Extraction: Can Correction Beat Encoding Overhead?

**Date:** 2026-03-31
**Status:** Complete (3/3 experiments)

## Motivation

Sprint 025 confirmed [[5,1,3]] basis isotropy on real IBM hardware and showed that hardware 2Q error rate (0.93%) is below the active-correction break-even threshold (1.86%). However, encoding-only (passive) QEC was insufficient — the encoding circuit's 10 CX gates add too much noise. The critical question: does active syndrome measurement + correction close the gap?

The practical challenge is that syndrome extraction requires 4 ancilla qubits and ~8 additional CX gates per round, adding depth and noise. The hypothesis: at current error rates, one round of syndrome extraction + correction should improve logical fidelity over passive encoding, because the error rate is below threshold.

## Literature Context

- Flag qubits on IBM hardware shown effective for syndrome extraction (Quantum journal, 2025)
- Circuit depth during syndrome measurement identified as key noise source (arXiv:2504.07258)
- Syndrome circuit scheduling optimization is active research (AlphaSyndrome, 2026)
- Gap: whether active correction helps at smallest possible scale (5+4 qubits, d=3) under realistic noise

## Experiments

### 26a: Syndrome Extraction Circuit — Noiseless Verification
*Build [[5,1,3]] syndrome extraction circuit with 4 ancilla qubits. Verify correct syndrome identification for all single-qubit Pauli errors.*

**Results:**
- Syndrome extraction circuit: 8 CX + 8 CZ = 16 two-qubit gates, depth 10
- All 15 single-qubit Pauli errors produce unique syndromes (distance-3 confirmed)
- 30/30 circuit-level syndrome measurements correct
- 30/30 full encode-error-correct cycles recover perfectly
- Combined circuit (encode + syndrome): 26+ two-qubit gates total — significant depth

### 26b: Noisy Active Correction vs Passive Encoding
*Gate-level noise from Sprint 025. Compare: (1) uncoded, (2) encode-only, (3) encode + 1 round syndrome + correction. Basis-averaged Holevo.*

**Results:**
- **Active correction is WORSE than passive encoding!**
- Uncoded (depth 1): avg Holevo 0.883, asymmetry 0.005
- Passive [[5,1,3]] (depth ~20): avg Holevo 0.511, asymmetry 0.012
- Active [[5,1,3]] (depth ~37): avg Holevo 0.331, asymmetry 0.079
- Correction delta: -0.006 (Z), -0.039 (X), -0.077 (Y) — correction HURTS in all bases
- ~79% of shots have trivial syndrome (0,0,0,0), but the remaining 21% often trigger WRONG corrections
- The syndrome extraction adds 16 more 2Q gates (total ~26), making the circuit so deep that noise during syndrome measurement corrupts the syndrome itself
- Correction backfire: applying the wrong correction based on a noisy syndrome INTRODUCES new errors
- Y-basis hit hardest (correction delta -0.077) — syndrome errors interact worst with Y measurements
- **Key insight:** At current noise levels (p2q=0.8%), syndrome extraction generates more errors than it corrects. The syndrome is too noisy to be useful.

### 26c: Depth-Noise Tradeoff and Break-Even Analysis
*Sweep p2q from 0.0001 to 0.02. For each, compare passive vs active [[5,1,3]].*

**Results:**
- **Active correction NEVER beats passive encoding — not even at p2q=0.0001 (100x below current hardware)**
- The gap widens with increasing noise: delta ranges from -0.06 (p2q=0.0001) to -0.24 (p2q=0.02)
- Uncoded always wins at every error rate tested (encoding overhead always exceeds benefit)
- Trivial syndrome fraction: 97.8% at p2q=0.0001, drops to 61.7% at p2q=0.02
- Even when syndromes are 97.8% trivial, the 2.2% non-trivial corrections still hurt

**Root cause: Error propagation through the syndrome circuit.**

The standard syndrome extraction circuit is NOT fault-tolerant. Each stabilizer measurement uses CX/CZ gates from the ancilla to 4 data qubits. A single error on the ancilla qubit *propagates* through these gates to create correlated multi-qubit errors on the data:

- Example: X error on g1 ancilla before CX gates → X errors on data qubits 0 AND 3 → weight-2 error
- Weight-2 errors are UNCORRECTABLE by a distance-3 code (can only correct weight-1)
- The syndrome extraction creates the exact class of errors the code cannot handle

This is the fundamental reason fault-tolerant syndrome extraction was invented. Without it, the syndrome circuit is a net error amplifier.

## Key Insights

1. **The "below threshold" claim from Sprint 025 was misleading.** The 0.93% < 1.86% comparison assumed ideal syndrome extraction. Real syndrome circuits have their own error budget, and non-fault-tolerant extraction amplifies errors beyond what the code can correct.

2. **Error propagation is the killer.** The syndrome circuit converts ancilla errors into multi-qubit data errors. A distance-3 code corrects weight-1 errors, but the syndrome circuit generates weight-2 errors from single faults. This is fundamentally worse than no correction at all.

3. **Fault-tolerant syndrome extraction is not optional — it's the entire game.** Flag qubits, verified extraction, or repeated measurement with majority vote are required to prevent error propagation. The "bare" syndrome circuit (what we implemented) will always hurt.

4. **The correction backfire asymmetry.** Y-basis correction hurts most (-0.077) because Y errors are the worst-case scenario for the correction logic — they flip the data bit AND are the hardest to classify correctly from a noisy syndrome.

5. **Scale matters doubly.** Not only does [[5,1,3]] have a small distance (d=3), but the syndrome circuit is a large fraction of the total circuit. At larger code distances, syndrome circuits are proportionally smaller relative to the error correction capacity.

## Connection to Prior Sprints

- Sprint 025: Hardware is "below threshold" only for ideal syndrome extraction. For realistic circuits, the effective threshold is much lower.
- Sprint 017-018: The threshold theorem assumes fault-tolerant gadgets at every level. Our non-FT syndrome extraction violates this assumption.
- Sprint 016: The "correction can backfire" finding for [[5,1,3]] under structured noise was a precursor to this much larger backfire effect.

## Next Steps

1. **Fault-tolerant syndrome extraction** — implement flag-qubit protocol for [[5,1,3]] (adds 1-2 flag qubits per stabilizer, catches dangerous error propagation)
2. **Repeated syndrome measurement** — majority-vote over 3 syndrome rounds to suppress measurement errors (triples circuit depth but may cross the correction threshold)
3. **Steane-style extraction** — use encoded ancilla for syndrome measurement (transversal = fault-tolerant, but requires preparing encoded ancilla states)
4. **Shift to larger codes** — at d=5 or d=7, syndrome circuit overhead is a smaller fraction of total error budget
