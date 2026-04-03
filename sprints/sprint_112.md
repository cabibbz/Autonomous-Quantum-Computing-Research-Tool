# Sprint 112 — q=12 χ_F Spectral: Power-Law vs Linear α(q) Discrimination

**Status:** In progress

## Motivation
Sprint 111 found power-law α(q)=0.69·q^0.69 slightly preferred over linear α(q)=0.260q+0.827 (ΔAIC=4.1), but discrimination was weak at q≤10. At q=12, the models diverge significantly:
- Linear: α = 0.260×12 + 0.827 = 3.947
- Power-law: α = 0.69×12^0.69 ≈ 4.50
- Difference: ~0.55 (vs ~0.15 at q=10)

q=12: n=4 (dim=20,736), n=5 (dim=248,832), n=6 (dim=2,985,984). All GPU-feasible.

## Experiments

### 112a — q=12 χ_F spectral decomposition at n=4,5,6
