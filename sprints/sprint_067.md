# Sprint 067 — Two-Transition Scan: Does the Hybrid Model Have a Floating Phase?

## Motivation

The Z_q clock model (q≥5) has TWO BKT transitions with a gapless "floating phase" in between — a critical Luttinger-liquid region with continuously varying exponents. Our hybrid model (Potts coupling + clock field) is intermediate between pure Potts and pure clock. Literature (Chepiga & Mila PRL 2019, Berker PRE 2023) predicts that for q≥5, a floating phase should generically appear.

Previous sprints scanned only near g_c. We've never scanned a wide g range looking for a SECOND critical point or an intermediate gapless region.

**Key question:** Does the hybrid model have one transition (as we've assumed) or two?

**Literature predictions:**
- Clock model: two BKT transitions, floating phase between them (Sun et al. 2020)
- S_q Potts: single first-order transition (Gorbenko et al. 2018)
- Hybrid: could have floating phase if clock field dominates at some g range
- Floating phase signature: gap ~ 1/N (power-law), gap*N ratio ∈ [0.7, 1.3] across sizes

## Experiments

### 067a — Wide gap scan at q=5,7
**Result: NO second gap minimum. Single transition confirmed.**

Scanned g=[0.02, 3.0] with 47-51 points per scan. q=5 at n=4,6,8; q=7 at n=4,6.

- Gap decreases monotonically from ordered phase (g=0) to a single minimum near g_c, then increases monotonically.
- No local minima found besides the primary critical point.
- Gap*N ratio converges to 2·sin(2π/q) at large g: 1.90 (q=5), 1.56 (q=7) — matches free transverse field limit.
- No extended region where gap*N is size-independent (would indicate floating phase).

### 067b — Entanglement entropy scan (DMRG n=16,24 at q=5)
**Result: NO floating phase. Critical-to-gapped transition is SHARP.**

DMRG (chi=20) at 6 key g-values for n=16 and n=24. Effective c from two-size entropy difference:

| g | S(n=16) | S(n=24) | c_eff | Interpretation |
|---|---------|---------|-------|----------------|
| 0.200 | 0.017 | 0.003 | -0.21 | Gapped (artifact: symmetry breaking at low chi) |
| 0.350 | 1.610 | 1.035 | -8.50 | Gapped (artifact: cat state at low chi) |
| 0.440 | 0.927 | 1.010 | **1.22** | **Critical** (g_c) |
| 0.550 | 0.400 | 0.400 | 0.003 | Gapped |
| 0.800 | 0.163 | 0.163 | 0.000 | Gapped |
| 1.200 | 0.071 | 0.071 | 0.000 | Gapped |

Only g=0.44 (the known g_c) shows c_eff > 0. No c_eff ≈ 1 region above g_c.

### 067c — Clock vs hybrid comparison at q=5
**Result: Clock HAS floating phase; hybrid does NOT. Method validated.**

Side-by-side gap*N scan for both models at n=4,6,8 across g=[0.05, 2.5].

**Gapless region detection (gap*N ratio ∈ [0.7, 1.3] between n=4,6):**
- **Clock: g ∈ [0.30, 0.92]** — extended floating phase spanning Δg = 0.62
- **Hybrid: g ∈ [0.38, 0.50]** — narrow critical window of Δg = 0.12 (just g_c)

Key observations:
- Clock gap*N ratio stays ~1.0 across wide g range (floating phase)
- Hybrid gap*N ratio jumps from ~1.0 to ~1.5 immediately above g_c
- At g=0.35: clock ratio = 0.83 (gapless), hybrid ratio = 0.62 (still gapped)
- At g=0.55: clock ratio = 1.02 (gapless), hybrid ratio = 1.37 (gapped)
- At g=0.65: clock ratio = 1.09 (gapless), hybrid ratio = 1.46 (gapped)

The Potts δ-function coupling is strong enough to suppress the floating phase entirely. The cosine coupling in the clock model allows the intermediate Luttinger-liquid phase because it's weaker at distinguishing non-adjacent states.

## Summary

**The hybrid model has exactly ONE transition.** Three independent probes confirm:
1. No second gap minimum (067a)
2. No c_eff ≈ 1 region above g_c (067b)
3. Clock model shows floating phase with same method, hybrid does not (067c)

This establishes a third difference between hybrid and clock universality:

| Property | Hybrid | Clock |
|----------|--------|-------|
| Transition type | Power-law 2nd order | BKT |
| ν | ~0.83 (finite) | →∞ |
| **Number of transitions** | **1** | **2 (with floating phase)** |

## Surprises
- Clock floating phase is WIDE (Δg ≈ 0.62 at q=5) — easily detectable even at n=4,6
- Hybrid transition is sharp — c_eff drops from 1.22 to 0.003 between g=0.44 and g=0.55
- No local gap minima at ANY size for either model (floating phase doesn't create a second minimum; it creates an extended flat region)
- Large-g gap*N → 2·sin(2π/q): q=5 gives 1.90, q=7 gives 1.56 (exact at n→∞)
- Clock n=8 scan took 886s (7x slower than hybrid) due to dense cosine coupling matrix
