# Sprint 019 — Quantum Channel Capacity: The Fundamental Limit of Error Correction

**Date:** 2026-03-31
**Status:** In Progress

## Motivation

We've established that the error correction threshold is a phase transition (Sprint 017), that different codes have different performance under structured noise (Sprint 016), and that syndrome information is lossy (Sprint 018). But we haven't asked the fundamental question: **what is the maximum rate at which quantum information can be reliably transmitted through a noisy channel?**

The answer is the **quantum channel capacity**, determined by the **coherent information** — the quantum generalization of mutual information. Coherent information tells us:
- Below what noise rate is error correction possible AT ALL?
- How much logical information can we protect per physical qubit?
- Why do different noise channels have fundamentally different correctable limits?

This connects channel capacity theory to everything we've built.

## Experiments

### 19a: Coherent Information of Quantum Channels
Compute coherent information I_coh = S(B) - S(AB) for depolarizing, amplitude damping, and phase damping channels at varying noise rates. The coherent information determines the quantum capacity Q = max_ρ I_coh(ρ, N).

### 19b: Channel Capacity vs Error Correction Threshold
Compare the theoretical channel capacity (where coherent information hits zero) to the error correction thresholds we measured in Sprint 017. The hashing bound should match the threshold.

### 19c: Code Rate vs Channel Capacity
For each code (3-qubit, [[5,1,3]], Steane, Shor), compare the achieved information rate (Holevo/n) to the channel capacity. How close do our codes get to the fundamental limit?

---

## Results

*(will be filled in as experiments complete)*
