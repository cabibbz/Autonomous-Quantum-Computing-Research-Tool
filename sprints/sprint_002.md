# Sprint 002 — CHSH Violation Landscape & Noise Sensitivity

**Date:** 2026-03-31
**Status:** Complete

## Idea

Sprint 001 showed S=2.834 at the optimal CHSH angles. Two questions follow:
1. What does the full angle landscape look like? Where exactly does violation occur?
2. How fragile is the violation? How much noise kills it?

Both experiments use only 2 qubits — well within resource limits.

## Experiments

### Experiment 2a: CHSH Angle Sweep

Computed exact CHSH S over a 50x50 grid of (a1, b0) angles with a0=0, b1=b0+pi/2.

Used the analytic formula E(a,b)=cos(a-b) for Phi+ state — verified against statevector at 3 test points (all match to 10 decimal places).

**Results:**
- S ranges from **-2.825 to +2.828** across the landscape
- Maximum at a1=1.539 (pi/2), b0=0.769 (pi/4) — the known optimal
- **Only 19.9% of the angle space violates the classical bound |S|>2**
- The violation region forms distinct "islands" in angle space
- At a1=pi, b0=pi/2: S=2.000 exactly — right at the classical boundary
- The landscape has 4-fold symmetry (S and -S regions are mirror images)

**Insight:** Quantum advantage in CHSH is not generic — it requires specific angle choices. ~80% of measurement settings don't violate the classical bound at all. This makes the violation more remarkable: nature supports correlations stronger than classical physics allows, but only if you know where to look.

### Experiment 2b: CHSH Under Depolarizing Noise

Swept depolarizing noise rate p from 0% to 30% on all gates (1q and 2q), measuring CHSH S at optimal angles with 20k shots per setting.

**Results:**
| Noise rate | S | Violates? |
|-----------|------|-----------|
| 0.00 | 2.826 | Yes |
| 0.05 | 2.374 | Yes |
| 0.09 | 2.040 | Barely |
| **0.095** | **2.000** | **Critical point** |
| 0.10 | 1.963 | No |
| 0.15 | 1.619 | No |
| 0.20 | 1.303 | No |
| 0.30 | 0.811 | No |

**Critical noise rate: ~9.5% depolarizing error per gate.**

At this point, quantum correlations become indistinguishable from classical ones.

**Insights:**
- Quantum advantage in CHSH is surprisingly fragile — only ~10% noise kills it
- The degradation is roughly exponential: S decays as ~(1-p)^k where k relates to circuit depth
- At p=0.30, S=0.81 — well into classical territory, approaching the uncorrelated limit of 0
- This has implications for real hardware: IBM's current error rates (~0.1-1% per gate) should allow CHSH violation, but more complex protocols may not survive
- The critical threshold is a hard physical boundary — no error mitigation strategy changes when S *actually* violates the bound, only when we can *detect* the violation through noise

### Experiment 2c: Entanglement Zoo (GHZ vs W vs Cluster vs Product)

Compared entanglement structure of 4 state families at n=4, 6, 8 qubits using three metrics:
1. Half-cut entropy (bipartite entanglement)
2. Single-qubit entropy (how entangled is one qubit with the rest)
3. Entropy after losing one qubit (robustness)

**Results at n=8:**

| State | Half-cut | Single-qubit | After qubit loss |
|---------|----------|-------------|-----------------|
| GHZ | 1.000 | 1.000 | 1.000 |
| W | 1.000 | 0.544 | 0.954 |
| Cluster | 1.000 | 1.000 | **2.000** |
| Product | 0.000 | 0.000 | 0.000 |

**Key findings:**

1. **Half-cut entropy is identical (1.0) for ALL entangled states.** This means half-cut entropy alone cannot distinguish qualitatively different entanglement structures. It's too coarse a measure.

2. **W state single-qubit entropy decreases with n** (0.81 at n=4, 0.54 at n=8). The entanglement is "diluted" across more qubits — each individual qubit is less entangled with the rest.

3. **Cluster state entropy INCREASES after qubit loss** (1.0 -> 2.0 at n>=6). This is remarkable: removing a qubit from a cluster state makes the remaining system MORE entangled. This is because measuring/tracing a qubit in a cluster state projects neighbors into entangled pairs — it's the basis of measurement-based quantum computation.

4. **GHZ is completely rigid** — all three metrics stay at exactly 1.0 regardless of n. Its entanglement structure doesn't change, it just involves more qubits.

**Surprise:** The cluster state result is the most interesting. In classical systems, losing a component always reduces correlations. In quantum cluster states, losing a qubit *creates* entanglement between its neighbors. This is a genuinely non-classical phenomenon with no analog in classical information theory.

## Sprint Summary

Three experiments, all fast (<30s), all revealing:
- **2a:** Only 20% of measurement angles yield quantum advantage in CHSH
- **2b:** ~9.5% noise kills quantum advantage — it's fragile
- **2c:** Different entangled states have radically different structures, but half-cut entropy can't tell them apart. Cluster states get MORE entangled when you lose qubits.

## What Surprised Me Most

The cluster state qubit-loss result. This suggests that entanglement can be a *resource that is released by destruction* — like breaking a chemical bond releases energy. This connects to measurement-based quantum computing (MBQC) where computation happens by destroying entanglement through measurement.

## Next Sprint Ideas

1. **Go deeper on cluster states:** Map how entropy changes as you trace out qubits 1 by 1. Is there a "phase transition" in entanglement?
2. **Entanglement under noise by state type:** Which state family is most robust to depolarizing noise?
3. **Quantum randomness:** Compare statistical properties of measurement outcomes from entangled vs product states
4. **Try real hardware:** Run CHSH on IBM QPU — compare real noise to our simulated depolarizing model
