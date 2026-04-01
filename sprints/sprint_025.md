# Sprint 025 — Real Hardware QPU Test: Does [[5,1,3]] Basis Isotropy Survive?

**Date:** 2026-03-31
**Status:** In Progress

## Motivation

24 sprints of simulator work established that [[5,1,3]]'s basis isotropy makes it the universally best small QEC code under combined T1+T2 noise. But all simulator models assume independent, identical noise on every qubit. Real hardware has:
- Spatially correlated noise
- Qubit-specific T1/T2 values
- Crosstalk between gates
- Readout errors
- Non-Markovian effects

The question: does isotropy survive these real-world complications?

## Plan

- **25a:** Simulator baseline — run exact encoding+measurement circuits under gate-level noise model matching IBM backend specs. Establish predictions.
- **25b:** Real hardware — same circuits on IBM QPU. Measure logical error rates per basis (Z, X, Y) for [[5,1,3]], 3-qubit repetition, and uncoded.
- **25c:** Gap analysis — compare simulator vs hardware. Where does the model break? What does the gap reveal about correlated noise?

## Literature Check

**Searched:** "[[5,1,3]] QEC real IBM hardware basis-dependent performance 2025-2026" on arXiv, Nature, APS

**Key findings:**
- Real IBM hardware shows basis-dependent QEC: Z-basis logical error ~0.040 vs X-basis ~0.088 (~2x gap) [Nature Comms 2023, subsystem codes]
- Noise-adapted 3-qubit code achieved break-even on IBM hardware under amplitude damping [arXiv:2603.04564]
- Comparative study of QEC on heavy-hexagonal lattice shows architecture-dependent performance [arXiv:2402.02185]
- No paper systematically compares basis-averaged Holevo for [[5,1,3]] vs repetition codes on real hardware

**Gap we're filling:** First basis-averaged Holevo comparison of [[5,1,3]] vs 3-qubit repetition on real IBM hardware. Our simulator predicts [[5,1,3]] isotropy (asymmetry ~0.006) vs 3-qubit anisotropy (~0.9). Real hardware's correlated noise may break this.

## Experiments

### 25a: Simulator Baseline with Gate-Level Noise

**Setup:** Built efficient [[5,1,3]] encoding circuit via Clifford tableau synthesis (10 CX gates, depth 16). Gate-level depolarizing noise model matching IBM hardware: p1q=0.0005, p2q=0.008, readout=0.015. 50k shots per circuit. Parity decoding for [[5,1,3]], majority (Z) / parity (X,Y) for 3-qubit.

**Technical note:** Y_L = -YYY for 3-qubit code (n=3), requiring inverted parity for Y-basis. Y_L = YYYYY for [[5,1,3]] (n=5), parity works directly. This sign from Y_L = i·(-iY)^n depends on n mod 4.

**Results:**

| Code | Z-Holevo | X-Holevo | Y-Holevo | Avg | Asymmetry | CX gates |
|------|----------|----------|----------|-----|-----------|----------|
| Uncoded | 0.887 | 0.887 | 0.885 | 0.886 | 0.002 | 0 |
| 3-qubit | 0.939 | 0.694 | 0.701 | 0.778 | 0.245 | 2 |
| [[5,1,3]] | 0.506 | 0.516 | 0.508 | 0.510 | 0.010 | 10 |

**Key findings:**
1. **[[5,1,3]] is 24x more isotropic than 3-qubit** (asymmetry 0.01 vs 0.24)
2. **Uncoded beats both codes** on average Holevo (0.886 > 0.778 > 0.510) — gate noise dominates
3. **3-qubit Z-basis is best single measurement** (0.939) but X/Y collapse to ~0.70
4. **[[5,1,3]] encoding overhead (10 CX)** costs ~10% error per CX, destroying its advantage
5. **The isotropy prediction survives** gate noise — [[5,1,3]] treats all bases equally even with circuit errors

### 25b: Real Hardware QPU

**Backend:** ibm_kingston (Heron, 156 qubits)
**QPU time:** 20 seconds (of 600s monthly budget)
**Shots:** 4096 per circuit, 18 circuits total

**T1/T2 for first 5 qubits:** T1 ~ 220-400 μs, T2 ~ 60-380 μs

| Code | Z-Holevo | X-Holevo | Y-Holevo | Avg | Asymmetry |
|------|----------|----------|----------|-----|-----------|
| Uncoded | 0.959 | 0.970 | 0.971 | 0.967 | 0.012 |
| 3-qubit | 0.976 | 0.815 | 0.722 | 0.838 | 0.254 |
| [[5,1,3]] | 0.482 | 0.522 | 0.492 | 0.499 | 0.040 |

**Key findings:**
1. **Isotropy prediction CONFIRMED**: [[5,1,3]] asymmetry 0.040 vs 3-qubit 0.254 — **6x more isotropic** on real hardware
2. **Hardware beats simulator model**: Uncoded Holevo 0.967 >> 0.886 predicted — IBM Heron has lower errors than our conservative noise model
3. **3-qubit asymmetry matches prediction**: HW 0.254 ≈ Sim 0.245 (3.7% difference)
4. **[[5,1,3]] asymmetry 4x worse on hardware** (0.040 vs 0.010): correlated noise partially breaks isotropy
5. **Uncoded still wins overall**: 0.967 > 0.838 > 0.499 — encoding overhead dominates at current error rates
6. **3-qubit Z-basis champion**: 0.976 — error rates so low that bit-flip correction provides real gain in Z-basis

### 25c: Gap Analysis

(Results pending)

## Key Findings

(To be filled)
