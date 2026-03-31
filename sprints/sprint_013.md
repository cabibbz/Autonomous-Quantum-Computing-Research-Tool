# Sprint 013 — Hayden-Preskill Protocol: Information Recovery from Quantum Scrambling

**Date:** 2026-03-31
**Status:** In Progress

## Motivation

Sprint 012 established that scrambling delocalizes information — GHZ can be maximally entangled yet zero-scrambled. The natural next question: once information IS scrambled, can it be recovered? The Hayden-Preskill thought experiment (2007) says yes, and with surprisingly little effort.

**The setup:** Alice throws a qubit into a "black hole" (a scrambling unitary). Bob holds the purification of the black hole's initial state ("early radiation"). After scrambling, Bob captures a small amount of "late radiation" from the output. The miracle: if the unitary is a good scrambler, Bob can recover Alice's qubit from just his early radiation + 1 qubit of late radiation — regardless of the black hole's size.

**Protocol (8 qubits total, within our 10-qubit limit):**
- Reference R: 1 qubit (entangled with Alice's message A)
- Alice's message A: 1 qubit
- Black hole B: 3 qubits (initially entangled with early radiation B')
- Early radiation B': 3 qubits (Bob holds these)
- Scrambling unitary U acts on A+B (4 qubits)
- Output split into "late radiation" D and "remaining black hole"
- Success metric: MI(R : B'∪D) → 2.0 means full recovery possible

## Experiments

### 13a: MI Recovery vs Scrambling Depth

**Setup:** 8 qubits — R(7), A(0), BH system(0-3), B'(4-6). All-to-all scrambling on BH. Track MI(R : B'∪D) for D = 1 or 2 qubits of late radiation.

**Results:**
| Layer | MI(R:B') | MI(R:B'∪D1) | MI(R:B'∪D2) | MI(R:D1 only) |
|-------|----------|-------------|-------------|---------------|
| 0     | 0.000    | 2.000       | 2.000       | 2.000         |
| 1     | 0.000    | 1.549       | 1.563       | 0.402         |
| 3     | 0.000    | 1.368       | 1.717       | 0.079         |
| 6     | 0.000    | 1.399       | 1.837       | 0.056         |
| 12    | 0.000    | 1.446       | 1.874       | 0.034         |

**Key observations:**
- **Layer 0:** MI(R:B'∪D1) = 2.0 trivially — D1 IS Alice's qubit (qubit 0), directly Bell-paired with R
- **After scrambling:** MI drops to ~1.4 with 1 qubit, ~1.87 with 2 qubits. Alice's info is delocalized.
- **MI(R:B') = 0 always** — early radiation alone never reveals Alice's message
- **MI(R:D1 only) → 0** — late radiation alone is useless. Bob NEEDS early radiation.
- At the decoupling threshold |D|+|B'| = |BH| (1+3=4), recovery is partial (~1.4/2.0 = 70%). With |D|=2 (above threshold), recovery reaches ~94%.

### 13b: Minimum Radiation for Recovery
[Results pending]

### 13c: Geometry Comparison — 1D vs All-to-All
[Results pending]

## Analysis
[Pending]
