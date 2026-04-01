# Sprint 035 — BW Size Scaling: Does Entanglement Hamiltonian Locality Sharpen with System Size?

**Date:** 2026-04-01
**Status:** In Progress

## Idea

Sprints 032-034 established that BW locality depends on the H/G-inv ratio:
- U(1) at d=2: 100% locality (XXZ, n=8)
- Z₂ at d=2: 91% locality (TFIM, n=8)
- S₃ at d=3: 76.5% locality (Potts, n=8)
- Z₃ at d=3: 69.7% locality (chiral clock, n=8)

But these are ALL n=8 results. The critical question: do these numbers converge to the H/G-inv prediction as n grows, or are they finite-size artifacts? TeNPy DMRG allows us to push to n=24+ for d=2 models.

**Literature:** Dalmonte et al. (2018) showed BW improves with system size for critical XXZ. Giudici et al. (2018) tested BW for Potts. Gap: no systematic BW locality vs system size across multiple symmetry groups.

## Experiments

### 35a: TFIM BW locality vs system size (n=8,12,16,20,24)
- Z₂ symmetry, d=2
- At criticality (h/J=1.0) and in ordered phase (h/J=0.5)
- Use exact diag for n≤10, DMRG for n>10
- Track BW locality = 1 - ||H_E - α·H_BW||²/||H_E||²

### 35b: XXZ BW locality vs system size (n=8,12,16,20,24)
- U(1) symmetry, d=2
- At Δ=0.5 (XY phase) and Δ=0.0 (XX point)
- Expect 100% locality to persist or improve

### 35c: MI-CV size scaling for TFIM and XXZ
- Track MI-CV order parameter at n=8,12,16,20,24
- Test whether MI-CV dome/inflection sharpens with system size
- Use DMRG 2-site RDMs for efficient MI computation

## Results

(To be filled as experiments complete)
