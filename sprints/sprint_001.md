# Sprint 001 — Bell States & CHSH Inequality

**Date:** 2026-03-31
**Status:** Partial (2 of 4 planned experiments completed; GHZ scaling timed out)

## Idea

Establish baselines: verify entanglement works on the simulator, measure CHSH inequality violation, and confirm tooling is functional.

## Experiments Run

### 1. Bell State Verification (4 states, 10k shots each)

Created all four Bell states and verified:

| State | Amplitudes | Counts (approx) | Entropy (bits) | Correlation |
|-------|-----------|-----------------|----------------|-------------|
| Phi+  | (|00>+|11>)/sqrt(2) | 00:5013, 11:4987 | 1.0 | +1.0 |
| Phi-  | (|00>-|11>)/sqrt(2) | 00:5012, 11:4988 | 1.0 | +1.0 |
| Psi+  | (|01>+|10>)/sqrt(2) | 01:4969, 10:5031 | 1.0 | -1.0 |
| Psi-  | (|01>-|10>)/sqrt(2) | 10:5026, 01:4974 | 1.0 | -1.0 |

All states maximally entangled (entropy = 1.0 bit exactly). Phi+/Phi- show perfect positive correlation (qubits always agree), Psi+/Psi- show perfect anti-correlation (qubits always disagree). Sampling distributions are 50/50 as expected.

### 2. CHSH Inequality Test (100k shots per setting)

Measured CHSH parameter S using optimal angles:
- Alice: a=0, a'=pi/2
- Bob: b=pi/4, b'=3*pi/4

| Setting | Correlation |
|---------|------------|
| E(a0,b0) | +0.7092 |
| E(a0,b1) | -0.7094 |
| E(a1,b0) | +0.7092 |
| E(a1,b1) | +0.7058 |

**S = 2.834** (classical bound: 2.0, Tsirelson bound: 2.828)

The simulator achieves the quantum maximum (slight overshoot is sampling noise from finite shots). This confirms genuine quantum correlations — no classical hidden variable model can produce S > 2.

## What Failed

- GHZ entropy scaling experiment (n=2..25) timed out. `partial_trace` on >15 qubits is too slow for CPU.
- Random vs GHZ entropy comparison also timed out for same reason.
- Running all experiments in one giant script meant losing results from completed experiments when later ones hung.

## Lessons Learned

1. **Keep qubit counts <= 10 when using partial_trace.** The cost is exponential.
2. **Save results after EACH experiment**, not at the end.
3. **Run experiments as separate small scripts** to avoid one timeout killing everything.
4. **Test timing on a single case before scaling up.**
5. The initial CHSH implementation had wrong Bob angles (b'=-pi/4 instead of 3*pi/4), giving S=0. Fixed by thinking through the correlation function E(a,b)=cos(a-b) for Phi+.

## What Surprised Me

- The simulator is extremely accurate for small circuits — CHSH S value matches Tsirelson bound to 3 decimal places with 100k shots.
- The CHSH angle bug was subtle: all 4 individual correlations were correct (~0.707) but the specific combination S cancelled to zero with wrong angles. Easy to miss if you only check individual measurements.

## Next Sprint Ideas

- Explore CHSH violation as a function of measurement angles (map the full landscape)
- Quantum random number generation: compare statistical properties of quantum vs pseudo-random
- Small-scale (<=8 qubit) entanglement entropy for different state families (GHZ, W, cluster)
- Noise effects: add depolarizing noise and see how CHSH degrades
