# Sprint 087 — Entanglement Spectrum Scaling via DMRG: Tail Weight Saturation?

## Status: In Progress

## Motivation
Sprint 084 revealed that walking breakdown = entropy concentration in the (q-1)-fold degenerate multiplet. But that was exact diag at n≤8. Key open question: **does the tail weight (levels ≥ 2) grow indefinitely with system size, or saturate?**

Sprint 086b showed q=7 tail grows 8× from n=6→12 (0.014%→0.11%). If tail keeps growing, it could eventually affect spectral observables too — the "entropy-only" breakdown would be temporary. If it saturates, entropy breakdown is genuinely decoupled from energy physics.

## Plan
- **087a**: DMRG entanglement spectrum at q=5, n=8,12,16,20,24. Track λ_max, (q-1) multiplet fraction, tail weight.
- **087b**: Same for q=7, n=8,10,12,16 (smaller sizes due to d=7 per site).
- **087c**: Analysis — saturation vs growth. Extrapolation. Scaling law for entropy fractions.

## Experiment 087a — q=5 entanglement spectrum scaling (n=8-24)
*(results pending)*

## Experiment 087b — q=7 entanglement spectrum scaling (n=8-16)
*(results pending)*

## Experiment 087c — Tail weight analysis and extrapolation
*(results pending)*
