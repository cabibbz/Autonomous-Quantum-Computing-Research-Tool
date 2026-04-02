# Sprint 109 — α(q) Correction Form: Power-Law vs Logarithmic

**Date:** 2026-04-02
**Status:** In progress

## Motivation

Sprint 108 fitted pairwise α drift vs 1/ln(N) and found q=3 α_∞ = 1.22, q=4 α_∞ = 1.66, q=5 α_∞ = 2.08. But for q=3 (continuous 3-state Potts, exact ν=5/6), the CFT prediction is α = 2/ν - 1 = 7/5 = 1.4. The extrapolated 1.22 is 13% below — a significant discrepancy.

The issue: 1/ln(N) corrections are appropriate for BKT (q=4) where multiplicative log corrections are known. But q=3 is a continuous transition with NO log corrections — its FSS corrections should be power-law: α_pair(N) = α_∞ + c/N^ω, where ω is the leading irrelevant exponent.

Key questions:
1. Does a power-law correction fit q=3 better than 1/ln(N)?
2. Does the power-law extrapolation recover the exact ν=5/6 (α_∞ = 1.4)?
3. Does q=4 (BKT) genuinely prefer the 1/ln(N) form?
4. Does q=5 (walking) remain correction-free under both models?

## Literature context

- q=3 Potts: ν = 5/6, ω_irr = 4/5 (leading irrelevant). FSS corrections ∝ 1/N^{4/5}.
- q=4 Potts: BKT with multiplicative log corrections. Gap ~ 1/(N·(ln N)^p).
- q≥5 Potts: complex CFT walking regime. Sprint 108: zero log corrections.

## Experiments

### 109a — q=3 GPU extension to n=12, 14

*(results pending)*

### 109b — q=4 n=11 and q=2 n=14 GPU extension

*(results pending)*

### 109c — Power-law vs logarithmic correction fits

*(results pending)*
