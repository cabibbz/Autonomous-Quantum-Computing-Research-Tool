# Sprint 123 — QPU Hardware Validation: TFIM Critical State

**Date:** 2026-04-08
**Thread:** Hardware validation
**Status:** Complete (simulation only — QPU credentials not configured)

## Motivation

580s QPU budget unused for 97 sprints. Goal: test simplest simulator prediction on real hardware.

**QPU justification:**
1. **Prediction:** TFIM n=4 at criticality: <ZZ>_nn = 0.653, <X>/n = 0.653
2. **Surprising if:** Hardware value deviates >20% (beyond noise model)
3. **Action:** Validates quantum critical state preparation on hardware

## Literature Search

Not needed — TFIM is well-understood. Focus on circuit design.

## Experiments

### 123a — Exact TFIM chi_F for n=4-8

TFIM at J=h=1 (critical). Exact chi_F scaling: alpha=1.15-1.29 for n=4-8 (approaching 1.0 with 1/N^2 corrections). <ZZ>_nn = 0.653, <X>/n = 0.653 at n=4. Gap*N = 1.59 (approaching pi*v_s = pi at large N).

### 123b — VQE optimization (FAILED)

Hardware-efficient ansatz (Ry-CNOT ladder) with 2-4 layers achieves only 0.50-0.56 fidelity for the critical ground state. The HEA lacks expressibility for the entangled critical state. Manual statevector simulation also had bugs in the CNOT gate implementation.

### 123c — Trotterized adiabatic prep (FAILED)

Even with 30 noiseless Trotter steps (240 CNOTs), the adiabatic circuit produces <ZZ>=0.33 vs exact 0.65. The total evolution time T=1 is far too short relative to 1/gap^2 ~ 6 at the critical point.

### 123d — Exact state preparation via initialize (SUCCESS)

**Key insight:** Qiskit's `initialize` decomposes the exact ground state into a circuit with only **11 CX gates** for n=4. This is far shallower than VQE or adiabatic approaches.

| Metric | Value |
|--------|-------|
| CX gates | 11 |
| Circuit depth | 32 |
| Signal at 1% CX error | 88% |
| Signal at 2% CX error | 79% |

### 123e — QPU submission (BLOCKED)

Prepared 6 circuits (3 coupling values x 2 measurement bases):
- h=0.5 (ordered): <ZZ>=0.922, <X>=0.291 — strong ZZ correlations
- h=1.0 (critical): <ZZ>=0.653, <X>=0.653 — equal correlations (self-dual point)
- h=2.0 (paramagnetic): <ZZ>=0.291, <X>=0.922 — strong X magnetization

Noisy simulation with 1% CX depolarizing error:

| h | Phase | <ZZ> noisy | <ZZ> exact | Signal |
|---|-------|-----------|-----------|--------|
| 0.5 | ordered | 0.82 | 0.92 | 89% |
| 1.0 | critical | 0.58 | 0.65 | 88% |
| 2.0 | paramagnetic | 0.27 | 0.29 | 91% |

**QPU submission blocked:** ~/.qiskit/qiskit-ibm.json is empty. Need to save IBM Quantum credentials before submitting.

## Key Findings

### 1. Exact state preparation is the right approach

VQE and adiabatic methods fail for the critical ground state at n=4. The `initialize` decomposition gives 11 CX gates — far better than VQE (9-20 CX at 50-56% fidelity) or adiabatic (40-240 CX at poor fidelity).

### 2. Phase transition is resolvable on noisy hardware

Even with 1% CX error, the noisy simulation clearly distinguishes the three phases:
- Ordered (h=0.5): <ZZ> >> <X>
- Critical (h=1.0): <ZZ> = <X> (self-duality)
- Paramagnetic (h=2.0): <X> >> <ZZ>

The self-dual symmetry <ZZ> = <X> at h_c is a hardware-robust signature of criticality.

### 3. QPU credentials needed

The IBM account credentials file is empty. To proceed:
```python
from qiskit_ibm_runtime import QiskitRuntimeService
QiskitRuntimeService.save_account(channel='ibm_quantum_platform', token='YOUR_TOKEN')
```

## Circuits Ready for QPU

The 6 circuits (3 h-values x 2 bases) are saved and tested. When credentials are configured, submission requires ~30s QPU time (well within 580s budget).

## Next Steps

1. **Configure IBM credentials** and submit the prepared circuits
2. **n=6 exact prep** — check CX count for `initialize` at n=6. If <30 CX, add to QPU run.
3. **Return to main thread** — the chi_F model comparison is nearly complete. Consider writing up the hybrid model results.
