# Sprint 067 — Two-Transition Scan: Does the Hybrid Model Have a Floating Phase?

## Motivation

The Z_q clock model (q≥5) has TWO BKT transitions with a gapless "floating phase" in between — a critical Luttinger-liquid region with continuously varying exponents. Our hybrid model (Potts coupling + clock field) is intermediate between pure Potts and pure clock. Literature (Chepiga & Mila PRL 2019, Berker PRE 2023) predicts that for q≥5, a floating phase should generically appear.

Previous sprints scanned only near g_c. We've never scanned a wide g range looking for a SECOND critical point or an intermediate gapless region.

**Key question:** Does the hybrid model have one transition (as we've assumed) or two?

**Predictions from literature:**
- Clock model: two BKT transitions, floating phase between them
- S_q Potts: single first-order transition
- Hybrid: could have floating phase if clock field dominates at some g range
- Floating phase signature: gap ~ 1/N (power-law, not exponential), continuously varying exponents

## Experiments

### 067a — Wide gap scan at q=5,7
**Result: NO second gap minimum. Single transition confirmed.**

Scanned g=[0.02, 3.0] with 47-51 points per scan. q=5 at n=4,6,8; q=7 at n=4,6.

- Gap decreases monotonically from ordered phase (g=0) to a single minimum near g_c, then increases monotonically.
- No local minima found besides the primary critical point.
- Gap*N ratio converges to 2·sin(2π/q) at large g: 1.90 (q=5), 1.56 (q=7) — matches free transverse field limit.
- No extended region where gap*N is size-independent (would indicate floating phase).

The Potts coupling (δ function) suppresses the floating phase entirely. Unlike the clock model's cosine coupling, the δ-function ordering is strong enough for a direct transition.

### 067b — Entanglement entropy across full phase diagram
[Results pending]

### 067c — Gap scaling in candidate floating region
[Results pending]

## Results
[Pending]

## Surprises
[Pending]
