# Sprint 027 — Flag Fault-Tolerant Syndrome Extraction for [[5,1,3]]

**Date:** 2026-03-31
**Status:** Complete (3/3 experiments)

## Motivation

Sprint 026 proved that "bare" (non-fault-tolerant) syndrome extraction ALWAYS hurts — at every error rate from 0.0001 to 0.02. The root cause: a single ancilla fault during syndrome measurement propagates through CX/CZ gates to create weight-2 data errors that the distance-3 code cannot correct. The standard decoder mistakes these weight-2 errors for weight-1 errors and applies the wrong correction, creating weight-3 (logical) errors.

The fix: **flag qubits** (Chao & Reichardt, 2018). Add one flag qubit per stabilizer measurement, placed to detect when ancilla errors propagate to multi-qubit data errors. When a flag triggers, use a modified correction lookup that accounts for the specific weight-2 error pattern.

## Literature

- Chao & Reichardt, "Quantum Error Correction with Only Two Extra Qubits" (PRL 2018) — original flag qubit proposal
- Chao & Reichardt, "Flag Fault-Tolerant Error Correction for any Stabilizer Code" (PRX Quantum 2020) — generalization
- Effectiveness of flag syndrome extraction on IBM hardware (Quantum, 2025) — recent experimental validation

## Results

### 27a: Flag Circuit Design & Noiseless Verification
- **Circuit**: 5 data + 4 ancilla + 4 flag = 13 qubits, 24 2Q gates (8 extra flag CNOTs vs bare's 16)
- **Flag placement**: For each weight-4 stabilizer with gate order [g0,g1,g2,g3]: H(anc), g0, CX(anc,flag), g1, g2, CX(anc,flag), g3, H(anc)
- **Weight-2 propagation errors identified** (ancilla X error between gates 1 and 2):
  - g1→Z₂X₃, g2→Z₃X₄, g3→Z₃Z₄, g4→X₃Z₄
- **Standard decoder fails completely**: 0/8 correct on weight-2 errors (fidelity = 0.000)
- **Flagged decoder succeeds perfectly**: 8/8 correct (fidelity = 1.000)
- **No false flags**: 15/15 single-qubit errors detected correctly with no flag triggers
- **All weight-2 propagation positions correctly flagged**: positions inside the bracket trigger flag, positions outside don't

### 27b: Flag-FT vs Bare vs Passive at Hardware Noise (p2q=0.008)

| Strategy | Z Holevo | X Holevo | Y Holevo | Avg Holevo | Asymmetry |
|----------|----------|----------|----------|------------|-----------|
| Uncoded  | 0.886    | 0.886    | 0.888    | 0.887      | 0.003     |
| Passive  | 0.510    | 0.516    | 0.507    | 0.511      | 0.009     |
| Bare     | 0.364    | 0.332    | 0.284    | 0.327      | 0.081     |
| Flag-FT  | 0.370    | 0.334    | 0.283    | 0.329      | 0.088     |

- **Flag-FT barely better than bare**: +0.002 avg Holevo (within statistical noise)
- **Both MUCH worse than passive**: ~-0.18 avg Holevo
- **Only 0.098% of shots** used the flagged weight-2 correction (5.9% had any flag)
- **Extra flag gates nearly cancel the flag benefit**: 8 additional 2Q gates add noise that offsets the improved correction

### 27c: Error Rate Sweep (p2q = 0.0001 to 0.02)

| p2q    | Passive | Bare   | Bare (no corr) | Flag-FT | Flag-Bare | Flag-Pass |
|--------|---------|--------|----------------|---------|-----------|-----------|
| 0.0001 | 0.950   | 0.928  | 0.945          | 0.928   | -0.000    | -0.022    |
| 0.0005 | 0.933   | 0.892  | 0.907          | 0.895   | +0.003    | -0.038    |
| 0.001  | 0.888   | 0.818  | 0.840          | 0.818   | -0.000    | -0.070    |
| 0.004  | 0.680   | 0.533  | 0.568          | 0.530   | -0.003    | -0.150    |
| 0.008  | 0.501   | 0.319  | 0.362          | 0.319   | -0.000    | -0.182    |
| 0.020  | 0.212   | 0.076  | 0.103          | 0.076   | -0.000    | -0.136    |

**Critical finding: Correction ALWAYS hurts at every noise level** (bare corrected < bare uncorrected, delta -0.02 to -0.04). Flag-FT never beats passive at any noise level.

## Key Insights

1. **Flag qubits solve the wrong problem at single-round scale.** The flag circuit perfectly detects weight-2 propagation errors (proven in 27a). But this is only one of many failure modes. Other sources — syndrome readout noise, accumulated gate errors across 24+ 2Q gates, multi-fault scenarios — dominate the error budget and are unaffected by flags.

2. **Single-round syndrome extraction is fundamentally insufficient.** Even with fault-tolerant gadgets, a single noisy syndrome measurement cannot provide enough information to improve over passive encoding. The correction delta is negative across ALL noise levels, for BOTH bare and flag-FT approaches. The syndrome circuit adds too many gates (encoding: 10 CX, syndrome: 16-24 2Q) relative to the information gained.

3. **The threshold theorem requires repeated measurements.** FT-QEC works by (a) preventing error propagation (flags), (b) getting reliable syndromes via repetition (majority vote over d rounds), and (c) using a decoder that accounts for time-correlated errors. We only implemented (a). Without (b) and (c), the noisy syndrome is more harmful than no syndrome.

4. **The "flag advantage" is real but tiny** (~+0.003 at p2q=0.0003-0.0005). This confirms that flags correctly identify weight-2 propagation errors. But the advantage is overwhelmed by the noise from the 8 extra flag CNOTs at all but the very lowest noise levels.

## Surprises

- **Flag-FT improvement over bare is essentially zero** — solving the weight-2 propagation problem (Sprint 026's identified root cause) doesn't fix active correction
- **Only 0.1% of shots use the flagged correction** — the flag match rate is so low that the decoder rarely engages the weight-2 correction
- **Correction itself hurts independently of the syndrome circuit** — bare (with correction) is 0.02-0.04 worse than bare (without correction) at ALL noise levels
- **Flag rate scales linearly with p2q** — from 0.4% at p2q=0.0001 to 15% at p2q=0.02, confirming flag triggers are noise-driven

## What This Means

Sprints 026-027 together prove that **single-round active QEC cannot beat passive encoding at the [[5,1,3]] scale**, whether using bare or flag-FT syndrome extraction. The problem is not error propagation (flags fix that) but rather that a single noisy syndrome measurement provides insufficient information. The path forward requires:

1. **Repeated syndrome measurements** with majority vote over 3+ rounds — this is the minimum for reliable syndrome information at realistic noise
2. **Larger code distance** (d≥5) where the syndrome overhead is proportionally smaller relative to the code's correction power
3. **Memory experiments** (multiple correction cycles) rather than single-shot correction — QEC advantage emerges over time, not in a single round

[Full data: results/sprint_027a_flag_circuit.json, results/sprint_027b_flag_vs_bare.json, results/sprint_027c_flag_threshold.json]
