# Sprint 019 — Quantum Channel Capacity: The Fundamental Limit of Error Correction

**Date:** 2026-03-31
**Status:** Complete (3/3 experiments)

## Motivation

We've established that the error correction threshold is a phase transition (Sprint 017), that different codes have different performance under structured noise (Sprint 016), and that syndrome information is lossy (Sprint 018). But we haven't asked the fundamental question: **what is the maximum rate at which quantum information can be reliably transmitted through a noisy channel?**

The answer is the **quantum channel capacity**, determined by the **coherent information** I_coh = S(B) - S(AB). This is the quantum generalization of Shannon's channel capacity theorem. It tells us the absolute physical limit on error correction.

## Experiments

### 19a: Coherent Information of Quantum Channels

Computed the quantum capacity (hashing bound) for three noise channels by optimizing coherent information over all input states.

**Channel capacity thresholds (where Q → 0):**
| Channel | Threshold | Theory |
|---------|-----------|--------|
| Depolarizing | p ≈ 0.20 | p ≈ 0.1893 |
| Amplitude damping | γ ≈ 0.50 | γ = 0.50 |
| Phase damping | λ ≈ 0.50 | Never reaches zero! |

**Key finding:** Maximally entangled input is near-optimal for depolarizing and phase damping (< 0.1% gap), but amplitude damping benefits slightly from optimized (non-maximally-entangled) input. This makes physical sense: amplitude damping breaks the symmetry between |0⟩ and |1⟩, so a biased input can extract more capacity.

**Capacity at p=0.02:**
- Amplitude damping: 0.919 bits/use (highest)
- Phase damping: 0.858 bits/use
- Depolarizing: 0.826 bits/use (lowest — most destructive)

### 19b: Channel Capacity vs Error Correction Threshold

Compared the channel capacity to actual code performance (fidelity and Holevo information) for the 3-qubit and [[5,1,3]] codes under depolarizing noise.

**Code efficiency at p=0.05:**
| Code | Rate k/n | Holevo/n | Fraction of Q(N) |
|------|----------|----------|-------------------|
| 3-qubit | 1/3 | 0.329 | 51% |
| [[5,1,3]] | 1/5 | 0.190 | 28% |

**Break-even points:**
- 3-qubit: never crosses below uncoded for depolarizing (Z-basis states)
- [[5,1,3]]: crosses below uncoded at p ≈ 0.14
- Channel capacity threshold: p ≈ 0.19

**The efficiency paradox:** At high noise (p > 0.12), the 3-qubit code's Holevo rate *exceeds* the channel capacity! This isn't a violation — the channel capacity bounds the *asymptotic* rate for *arbitrary* states. A fixed-length code protecting one specific basis can exceed this because:
1. It doesn't need to protect arbitrary states
2. It can't be concatenated to arbitrarily high fidelity
3. The hashing bound is only a lower bound for depolarizing (non-degradable channel)

### 19c: Channel Capacity Landscape

Full comparison of quantum capacity, entanglement-assisted capacity, and cross-channel behavior.

**Capacity ordering at p=0.1:**
```
Q(amplitude_damping) = 0.707 > Q(phase_damping) = 0.528 > Q(depolarizing) = 0.368
```

**Phase damping anomaly:** The capacity is *non-monotonic* — it drops to near zero around λ≈0.5 then RISES again! At λ=1.0, dephasing becomes deterministic Z, which is a unitary operation (fully reversible). This means:
- λ = 0.0: perfect channel (Q = 1.0)
- λ ≈ 0.5: worst case (Q ≈ 0)
- λ = 1.0: deterministic Z rotation (Q = 1.0!)

**Entanglement-assisted capacity (C_EA) at p=0.1:**
| Channel | Q | C_EA | Boost factor |
|---------|------|------|-------------|
| Depolarizing | 0.37 | 1.37 | 3.7x |
| Amplitude damping | 0.71 | 1.71 | 2.4x |
| Phase damping | 0.53 | 1.53 | 2.9x |

Pre-shared entanglement provides a 2-4x boost in communication rate. The boost is largest for the most noisy channel (depolarizing) — entanglement is most valuable when the channel is most destructive.

**Capacity ordering predicts code performance (Sprint 016):**
- Q(amp_damping) > Q(depolarizing) ✓ — codes do better under amplitude damping
- BUT Q(phase_damping) > Q(depolarizing) while codes do WORSE under phase damping!
- This anomaly reveals that our codes are mismatched to phase noise — capacity exists but current codes can't reach it

## Surprises

1. **Phase damping capacity is non-monotonic** — worst at λ≈0.5, perfect at both extremes. Maximum dephasing is fully reversible because it's deterministic.
2. **3-qubit code "exceeds" channel capacity** — a finite code protecting a specific basis can beat the asymptotic bound for arbitrary states. The hashing bound isn't tight for non-degradable channels.
3. **Entanglement boost is largest when channel is worst** — depolarizing gets 3.7x from pre-shared entanglement vs 2.4x for amplitude damping. Entanglement is more valuable when noise is more destructive.
4. **Capacity predicts code performance... except when it doesn't** — the ordering Q(amp) > Q(phase) > Q(depol) matches code fidelity for amp vs depol, but phase damping has higher capacity yet worse code performance. The gap is code architecture, not channel physics.
5. **Codes operate far below capacity** — even the 3-qubit code at 51% efficiency is remarkable; [[5,1,3]] at 28% shows enormous room for improvement. Shannon's theorem guarantees codes exist that approach 100%.

## Key Insight

The quantum channel capacity reveals three distinct regimes of error correction:

1. **The capacity gap** (Sprint 017 threshold vs channel capacity threshold): Our codes break even well below the fundamental limit. The 3-qubit code breaks even at all depolarizing rates (for Z-basis), the [[5,1,3]] at p≈0.14, but the channel allows correction up to p≈0.19. Better codes can close this gap.

2. **The code-channel mismatch** (Sprint 016 connection): Phase damping has *more* capacity than depolarizing, yet our codes perform *worse* under it. The capacity exists in the channel but our codes can't access it. This is the information-theoretic explanation for why Steane's correction "backfires" under phase noise — it's not that phase noise is worse, it's that the code's error model doesn't match the channel.

3. **The entanglement dividend**: Pre-shared entanglement always helps, but helps most when the channel is most destructive. This connects to the Hayden-Preskill protocol (Sprint 013) — early radiation (pre-shared entanglement) is what enables information recovery from a scrambled black hole.

**Unifying thread through Sprints 012-019:** Information has a *capacity to survive* through noisy quantum channels, determined by how much entanglement the channel preserves. Error correction codes are designed scramblers (Sprint 014) that spread information across qubits so the channel can't concentrate its damage. The threshold (Sprint 017) is where the channel's damage rate exceeds the code's spreading rate. The channel capacity is the ultimate speed limit — the maximum spreading rate any code could achieve.

## Next

- Toric/surface code entanglement structure (connects 2D cluster topology to QEC capacity)
- Combined T1+T2 noise model (realistic hardware noise)
- Capacity-achieving codes (how close can we get to Q(N) with small codes?)
- Real hardware comparison (the gap between simulated and real capacity)
