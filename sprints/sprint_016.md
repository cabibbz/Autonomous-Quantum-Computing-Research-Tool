# Sprint 016 — Structured Noise on QEC Codes: Breaking Topology-Blindness

**Date:** 2026-03-31
**Goal:** Test whether structured noise (amplitude damping, phase damping) breaks the topology-blindness observed under depolarizing noise in Sprint 015. Each code architecture should have a unique noise vulnerability profile.

**Hypothesis:** Depolarizing noise treats all codes equally (~0.03 per-qubit efficiency). Structured noise should discriminate:
- **Shor** (GHZ-like blocks): should survive dephasing (Sprint 009 showed GHZ Phi survives dephasing) but be vulnerable to amplitude damping
- **Steane** (CSS structure): symmetric X/Z stabilizers should give balanced noise response
- **[[5,1,3]]** (perfect code): should be most balanced — perfect code corrects all single-qubit errors equally

**Connection:** Sprint 009 (noise fingerprints on entanglement archetypes) + Sprint 015 (code families under depolarizing) = this sprint.

---

## Experiment 16a: Amplitude Damping on Code Families

*In progress...*
