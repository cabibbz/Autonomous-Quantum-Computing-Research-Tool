# Sprint 028 — Repeated Syndrome Measurement: Gate Overhead Kills

**Date:** 2026-03-31
**Goal:** Test whether repeated syndrome extraction (3+ rounds) with majority vote provides enough syndrome reliability for active correction to beat passive encoding on [[5,1,3]].

**Literature check:**
- Searched: "repeated syndrome measurement majority vote quantum error correction [[5,1,3]] fault tolerant threshold 2025 2026"
- Key findings: Standard FT protocol requires d rounds for distance-d. Google's below-threshold result uses repeated measurement on surface codes. No specific study on whether this works at [[5,1,3]] scale.
- Gap: Whether repeated measurement can overcome the single-round limitation proven in Sprints 026-027.

**Hypothesis:** 3 rounds of syndrome measurement with majority vote should improve syndrome reliability enough for correction to beat passive.

**Result:** WRONG. Repeated syndrome measurement makes things **worse**, not better. Gate overhead always exceeds correction benefit.

## Experiments

### 28a: Repeated syndrome circuit — build and verify
**Status:** Complete

Built 3-round repeated syndrome extraction circuit:
- 17 qubits total (5 data + 12 ancilla, fresh ancillas each round)
- 48 two-qubit gates (16 per round × 3)
- Depth 30

Noiseless verification:
- All 30 single-qubit data errors correctly identified by majority vote (30/30)
- All 12 ancilla readout errors correctly overridden by majority vote (12/12)
- Circuit works perfectly in the noiseless limit

### 28b: Repeated syndrome vs single-round vs passive
**Status:** Complete

At hardware noise (p1q=0.0005, p2q=0.008, p_readout=0.015):

| Strategy | Avg Holevo | Z | X | Y | Asymmetry |
|---|---|---|---|---|---|
| Uncoded | 0.886 | 0.887 | 0.884 | 0.887 | 0.003 |
| Passive [[5,1,3]] | 0.510 | 0.509 | 0.511 | 0.510 | 0.002 |
| Active 1-round | 0.329 | 0.365 | 0.334 | 0.288 | 0.077 |
| Active 3-round majority | 0.273 | 0.288 | 0.270 | 0.263 | 0.026 |
| Active 3-round nocorr | 0.200 | 0.202 | 0.199 | 0.198 | 0.004 |

**Key findings:**
- 3-round majority (0.273) is WORSE than 1-round (0.329), which is WORSE than passive (0.510)
- Each additional round degrades performance — the extra 32 2Q gates per 2 rounds introduce more errors than the improved syndrome can fix
- Syndrome agreement across rounds is only 59.4% — rounds disagree 40% of the time
- Without correction, 3-round circuit has Holevo 0.200 vs passive 0.510 — the syndrome extraction gates alone destroy information

**Surprise:** 3-round majority has BETTER isotropy (0.026) than 1-round (0.077) — majority vote does improve syndrome quality, but the gate overhead overwhelms this benefit.

### 28c: Error rate and rounds sweep
**Status:** Complete

Swept p2q from 0.0001 to 0.02 and round counts 1, 3, 5:

| p2q | Passive | 1-round | 3-round | Δ(3r-passive) |
|---|---|---|---|---|
| 0.0001 | 0.953 | 0.929 | 0.939 | -0.014 |
| 0.0005 | 0.941 | 0.910 | 0.884 | -0.057 |
| 0.001 | 0.898 | 0.855 | 0.822 | -0.076 |
| 0.004 | 0.714 | 0.599 | 0.527 | -0.187 |
| 0.008 | 0.540 | 0.404 | 0.307 | -0.233 |
| 0.020 | 0.276 | 0.135 | 0.068 | -0.208 |

**Key findings:**
- Active correction NEVER beats passive at ANY noise level (p2q=0.0001 to 0.02)
- More rounds ALWAYS worse (except 3-round marginally beats 1-round at p2q=0.0001)
- Gate overhead analysis: passive=10 CX, 1-round=26 CX, 3-round=58 CX, 5-round=90 CX
- At p2q=0.008: expected CX errors are 0.08 (passive), 0.21 (1-round), 0.46 (3-round), 0.72 (5-round)
- The code can correct at most 1 error; the syndrome circuit reliably produces more than 1

## Key Insight

**Gate overhead is the fundamental barrier to active QEC at [[5,1,3]] scale.**

The [[5,1,3]] syndrome extraction requires 16 two-qubit gates per round (4 stabilizers × 4 data interactions). Each round introduces ~16p errors in expectation. For p2q=0.008, one round adds 0.13 expected errors — and the code can only correct 1 total error. Even with perfect syndrome reliability (infinite rounds, perfect majority vote), the accumulated gate damage from syndrome extraction produces uncorrectable multi-qubit errors faster than the syndrome can identify single-qubit errors.

This resolves the complete Sprint 026-028 arc:
- **Sprint 026:** Non-FT syndrome fails (error propagation)
- **Sprint 027:** Flag-FT syndrome fails (single round insufficient)
- **Sprint 028:** Repeated syndrome fails (gate overhead exceeds correction capacity)

The conclusion: **active error correction on [[5,1,3]] requires the physical error rate to be much lower than the syndrome circuit's per-gate error rate divided by the number of syndrome gates.** For 16 gates/round, the break-even requires p2q << 1/16 ≈ 0.06, but even at p2q=0.0001 (600× below this), correction still doesn't help because of readout noise and the accumulated effect of encoding gates.

**The real path to QEC advantage requires:**
1. **Larger codes** where distance grows faster than gate count (d=5,7,... surface codes)
2. **Lower per-gate error rates** — current ~0.1% is insufficient for 5-qubit codes
3. **Both together** — the threshold theorem works in the asymptotic limit, not at d=3

This is why Google's below-threshold result uses distance 5 and 7 surface codes with 0.14% per-cycle error rates — much more favorable gate-to-distance ratio than [[5,1,3]].

## Next directions
- **Memory experiment**: multiple correction cycles to measure error accumulation rate (even if each cycle hurts, the rate tells us about the threshold)
- **Steane-style extraction**: encoded ancillas that avoid the gate overhead problem
- **Shift to new territory**: variational circuits, QRNG, or other non-QEC exploration (the QEC small-scale story is now complete)
- **Analytical bound**: compute the exact per-gate threshold for [[5,1,3]] active correction
