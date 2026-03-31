# Sprint 005 — Quantum Discord: Separating Quantum from Classical Correlations

**Date:** 2026-03-31
**Status:** In progress

## Motivation

We've built an information-theoretic hierarchy: entropy < single-qubit entropy < MI < I3. Each level reveals structure the previous misses. The next natural measure is **quantum discord** — it decomposes mutual information into classical and quantum parts.

Key question: GHZ has positive I3 (redundant, "classical-like"). Does it have low discord too? If yes, GHZ correlations are genuinely classical in character despite being "maximally entangled." If discord is high, then I3 is misleading about classicality.

Discord formula for bipartite state rho_AB:
- D(A|B) = S(rho_A) - S(rho_AB) + min_{measurements on B} S(A | post-measurement)
- The min is over all projective measurements on B (Bloch sphere angles theta, phi)
- Discord = MI - Classical correlations

## Experiments

### 5a: Pairwise Quantum Discord
- Compute discord for all qubit pairs in GHZ, W, Cluster (n=6)
- Compare discord magnitude and spatial structure
- Expected: GHZ low discord (classical-like), Cluster high discord for neighbors

### 5b: Discord vs MI Ratio
- What fraction of MI is quantum (discord) vs classical?
- Compare across state types
- This gives us the "quantumness fraction" of each correlation

## Results

(to be filled in)
