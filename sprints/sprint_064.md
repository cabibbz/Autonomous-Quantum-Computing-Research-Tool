# Sprint 064 — Entanglement Hamiltonian at q>4 Potts: BW Locality in the Novel CFT

## Motivation
Sprints 032-034 established that BW locality is controlled by the H/G-inv ratio (Hamiltonian operator dimension / symmetry-invariant dimension). Measured for d=2 (Z₂, U(1)) and d=3 (S₃, Z₃). The novel q>4 Potts CFT (Sprints 054-063) has d=q with Z_q symmetry — the H/G-inv ratio should decrease with q, predicting degrading BW locality.

**Key question:** Does the entanglement Hamiltonian remain local at the novel q>4 critical points?

**Predictions:**
- BW locality should decrease monotonically with q (from 76.5% at q=3)
- H/G-inv ratio should predict the ordering
- The entanglement temperature profile β(x) should show Unruh-like gradient for all q

## Experiments

### Experiment 064a — BW locality for q=4 Potts at g_c=0.392
*Status:* Complete

**System:** n=8, n_A=4, dim_full=65536, dim_A=256. Sweep g=0.05-1.0 (15 points).

**Results:**
- **Peak BW locality: 56.0% at g=0.30** (ordered phase, NOT at g_c=0.392)
- At g_c: locality=50.1%, BW[linear]=48.6%
- Linear envelope wins at ALL g values (sin_inv slightly worse, uniform worst)
- Coupling profile at g_c: Unruh-like gradient preserved (delta bonds 5.80→4.10→2.21, XXd sites 1.83→1.50→0.90→0.27)
- Entropy at g_c: S=0.933 bits (vs q=3 S=1.585 at ordered phase)

**Key finding:** q=4 peak locality (56%) is dramatically below q=3 (76.5%). 20 percentage point drop from q=3→4. The Unruh-like gradient persists but captures a much smaller fraction of H_E.

**Surprise:** Peak locality at g=0.30, well BELOW g_c=0.392. At q=3, peak was near g_c. The ordered phase deep in the 4-fold degeneracy region has more BW-like structure than the critical point itself.

### Experiment 064b — BW locality for q=5,7,10 Potts
*Status:* Complete

**Systems:** q=5 n=8 (dim_A=625), q=7 n=6 (dim_A=343), q=10 n=6 (dim_A=1000).

**Peak BW locality:**
| q | n | n_A | Peak locality | Peak g | At g_c |
|---|---|-----|--------------|--------|--------|
| 5 | 8 | 4 | 42.3% | 0.30 | 36.2% |
| 7 | 6 | 3 | 69.4% | 0.20 | 60.2% |
| 10 | 6 | 3 | 51.2% | 0.30 | 38.7% |

**At same n_A=4 (fair comparison): q=2(91%) > q=3(76.5%) > q=4(56%) > q=5(42.3%) — monotonic decrease.**

**At same n_A=3: q=7(69.4%) > q=10(51.2%) — also monotonic.**

The n_A=3 values are inflated relative to n_A=4 (fewer sites = fewer non-Hamiltonian operators). Cannot directly compare q=5 (n_A=4) to q=7 (n_A=3).

**Coupling profiles at g_c preserve Unruh-like gradient for ALL q.** Bond couplings decrease monotonically toward entanglement cut. Site couplings also decrease. Gradient becomes steeper with q.

**Surprises:**
- Peak locality always in ordered phase (g < g_c), never at criticality
- Linear envelope wins for q≥4 (sin_inv won for q=2,3 in Sprint 032-033)
- At n_A=4: locality drops ~15% per unit increase in q (91→76→56→42)

### Experiment 064c — H/G-inv ratio analysis and scaling with q
*Status:* Complete

**Method:** Count Z_q-invariant traceless operators on n_A sites with local dim q. An operator is Z_q-invariant iff it commutes with X_total = X^{⊗n_A}. Equivalently, it must be block-diagonal in the q charge sectors. Verified by explicit eigendecomposition for all q=2-10.

**Key result: Z_q-invariant fraction = exactly 1/q.**
| q | n_A | G-inv/Total | H_terms | H/G-inv |
|---|-----|-------------|---------|---------|
| 2 | 4 | 49.8% (1/2) | 7 | 0.0551 |
| 3 | 4 | 33.3% (1/3) | 7 | 0.0032 |
| 4 | 4 | 25.0% (1/4) | 7 | 0.00043 |
| 5 | 4 | 20.0% (1/5) | 7 | 0.00009 |
| 7 | 3 | 14.3% (1/7) | 5 | 0.00030 |
| 10 | 3 | 10.0% (1/10) | 5 | 0.00005 |

**H/G-inv scales as ~5/q^5 at n_A=3** (power law fit: exponent = -5.016, coefficient = 5.157, max error 1.2%). The ratio plummets because G-inv operators grow as q^{2n_A-1} while H terms are fixed at 2n_A-1.

**H/G-inv perfectly predicts BW locality ordering.** At fixed n_A=4: q=2 (0.055, 91%) > q=3 (0.003, 76.5%) > q=4 (0.0004, 56%) > q=5 (0.0001, 42.3%). Monotonic.

**Physical explanation:** As q grows, the Z_q symmetry constrains less (only 1/q of operators are forbidden), leaving a vast space of non-Hamiltonian operators that H_E can populate. The physical Hamiltonian terms become a vanishingly small fraction of the symmetry-allowed space.

**Note:** These Z_q H/G-inv values differ from Sprint 034's because Sprint 034 used different symmetry groups (physical symmetry of TFIM is Z_2 but acted differently on Pauli space; Potts q=3 had S_3 not Z_3). The current computation is self-consistent for pure Z_q cyclic symmetry.

**Surprise:** G-inv fraction = exactly 1/q — a clean algebraic identity. Each Z_q charge sector has dimension q^{n_A}/q = q^{n_A-1}, so inv = q·(q^{n_A-1})^2 = q^{2n_A-1}, fraction = 1/q.

## Summary

Sprint 064 extends the Bisognano-Wichmann entanglement Hamiltonian analysis (Sprints 032-034) to the novel q>4 Potts CFT. Three experiments measured BW locality for q=4,5,7,10 and analyzed the H/G-inv operator ratio.

**All three predictions confirmed:**
1. BW locality decreases monotonically with q (at fixed n_A): 91% → 76.5% → 56% → 42.3%
2. H/G-inv ratio predicts the ordering perfectly
3. Unruh-like entanglement temperature gradient persists for ALL q

**Grand BW locality table (at fixed n_A=4 where available):**
| q | Symmetry | BW Peak | At g_c | H/G-inv |
|---|----------|---------|--------|---------|
| 2 | Z₂ | 91.0% | 90.6% | 0.055 |
| 3 | S₃ | 76.5% | ~76% | 0.003 |
| 4 | Z₄ | 56.0% | 50.1% | 0.0004 |
| 5 | Z₅ | 42.3% | 36.2% | 0.00009 |

**Additional (n_A=3, not directly comparable):**
| 7 | Z₇ | 69.4% | 60.2% | 0.0003 |
| 10 | Z₁₀ | 51.2% | 38.7% | 0.00005 |

**Surprises:**
1. **Peak BW locality shifts into ordered phase for q≥4.** For q=2,3, peak is near g_c; for q≥4, peak is at g ≈ 0.30 regardless of g_c value. The ordered phase has more BW-like structure.
2. **Linear envelope wins for q≥4** (vs sin_inv for q=2,3). The finite-size BW envelope changes character at the q=4 boundary.
3. **~15% drop per unit q** at n_A=4: (91→76→56→42). If this continues linearly, BW locality approaches 0 near q~8 at n_A=4.
4. **G-inv fraction = exactly 1/q** — the Z_q symmetry eliminates exactly 1-1/q of operators. Simple but previously uncomputed for general q.
5. **H/G-inv ~ 5/q^5** — the physical Hamiltonian terms become a vanishingly small fraction of symmetry-allowed operators as q grows.

**Implications for the novel CFT:** The entanglement Hamiltonian at q>4 critical points is NOT well-approximated by the physical Hamiltonian. H_E contains substantial non-local corrections that are symmetry-allowed but absent from H. This does not mean the CFT is non-local — it means the BW approximation breaks down because the operator space is too large for Z_q to constrain effectively.

**Connection to c>1:** The growing central charge c(q)~ln(q) means more effective degrees of freedom at criticality. These additional DOF manifest as non-Hamiltonian operators in H_E, degrading BW locality. The correlation between decreasing BW locality and increasing c(q) suggests a fundamental tradeoff: richer CFTs have less local entanglement Hamiltonians.
