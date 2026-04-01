# Sprint 027 — Flag Fault-Tolerant Syndrome Extraction for [[5,1,3]]

**Date:** 2026-03-31
**Status:** In Progress

## Motivation

Sprint 026 proved that "bare" (non-fault-tolerant) syndrome extraction ALWAYS hurts — at every error rate from 0.0001 to 0.02. The root cause: a single ancilla fault during syndrome measurement propagates through CX/CZ gates to create weight-2 data errors that the distance-3 code cannot correct. The standard decoder mistakes these weight-2 errors for weight-1 errors and applies the wrong correction, creating weight-3 (logical) errors.

The fix: **flag qubits** (Chao & Reichardt, 2018). Add one flag qubit per stabilizer measurement, placed to detect when ancilla errors propagate to multi-qubit data errors. When a flag triggers, use a modified correction lookup that accounts for the specific weight-2 error pattern.

## Literature

- Chao & Reichardt, "Quantum Error Correction with Only Two Extra Qubits" (PRL 2018) — original flag qubit proposal
- Chao & Reichardt, "Flag Fault-Tolerant Error Correction for any Stabilizer Code" (PRX Quantum 2020) — generalization
- Effectiveness of flag syndrome extraction on IBM hardware (Quantum, 2025) — recent experimental validation

Key insight: For distance-3 codes with weight-4 stabilizers, one flag qubit per stabilizer suffices. The flag CNOT brackets the middle gates where weight-2 propagation occurs. This converts O(p) logical errors from propagation into O(p²) (detected and correctly handled).

## Experimental Plan

### 27a: Flag Circuit Design & Noiseless Verification
- Build flag-FT syndrome circuit: 5 data + 4 ancilla + 4 flag = 13 qubits
- Verify flag correctly detects ancilla errors that create weight-2 data errors
- Build modified correction lookup table for flagged syndromes
- Verify correction succeeds for all flagged error patterns

### 27b: Flag-FT vs Bare vs Passive Under Noise
- Compare three strategies at hardware noise level (p2q=0.008)
- Basis-averaged Holevo as figure of merit
- Key question: does the flag information overcome the extra gate overhead (24 vs 16 2Q gates)?

### 27c: Error Rate Sweep — Finding the Flag-FT Threshold
- Sweep p2q from 0.0001 to 0.02
- Find crossover where flag-FT beats passive encoding
- Compare scaling: bare is O(p) logical error, flag-FT should be O(p²)

---

## Results

*(to be filled as experiments complete)*
