# Sprint 026 — Active Syndrome Extraction: Can Correction Beat Encoding Overhead?

**Date:** 2026-03-31
**Status:** In progress

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

**Results:** (pending)

### 26b: Noisy Active Correction vs Passive Encoding
*Gate-level noise from Sprint 025. Compare: (1) uncoded, (2) encode-only, (3) encode + 1 round syndrome + correction. Basis-averaged Holevo.*

**Results:** (pending)

### 26c: Depth-Noise Tradeoff and Break-Even Analysis
*Sweep gate error rate. Find crossover where active correction beats passive. Analyze how many syndrome rounds help vs hurt.*

**Results:** (pending)
