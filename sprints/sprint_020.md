# Sprint 020 — Topological Codes: The Toric Code's Entanglement Structure

**Date:** 2026-03-31
**Status:** In Progress

## Motivation

The toric code has been the top "next" item since Sprint 007 when we discovered that 2D cluster states show topological entanglement protection. Sprints 014-019 explored algebraic QEC codes ([[5,1,3]], Steane, Shor) and their information structure. Now we ask: how does a genuinely *topological* code differ?

The toric code on a 2×2 torus encodes 2 logical qubits in 8 physical qubits (distance 2). It's the simplest code with topological order — its logical operators are non-contractible loops on the torus, and no local measurement can distinguish the 4 ground states. This is qualitatively different from algebraic codes where logical operators are algebraic constructions.

**Key questions:**
1. How does the toric code's entanglement structure (MI, I3, negativity) compare to algebraic codes?
2. Can we measure the topological entanglement entropy (TEE) — a genuinely topological invariant?
3. Does the Page curve show topological features distinct from algebraic codes?

## Experiments

### 20a: Toric Code Entanglement Structure
**Goal:** Prepare [[8,2,2]] toric code ground state, compute MI matrix, I3, negativity spectrum. Compare with [[5,1,3]], Steane, Shor.

*(results pending)*

### 20b: Topological Entanglement Entropy
**Goal:** Compute Kitaev-Preskill TEE using the tripartite construction S_topo = S_A + S_B + S_C - S_AB - S_BC - S_AC + S_ABC. For toric code, expect S_topo = -ln(2) ≈ -0.693.

*(results pending)*

### 20c: Page Curve and Noise Fingerprint
**Goal:** Page curve for subset recovery of logical information. Noise response under depolarizing, amplitude damping, phase damping. Compare with algebraic code archetypes.

*(results pending)*
