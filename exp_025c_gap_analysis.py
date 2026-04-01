"""
Sprint 025c: Gap analysis — simulator vs hardware
Quantify what our noise model gets right, what it gets wrong, and what
the deviation reveals about correlated noise on real hardware.
"""

import numpy as np
import json, time

start = time.time()

# Load results
with open('results/sprint_025a_hardware_baseline.json') as f:
    sim_data = json.load(f)
with open('results/sprint_025b_hardware_qpu.json') as f:
    hw_data = json.load(f)

sim = sim_data['results']
hw = hw_data['results']

# ============================================================
# 1. Direct comparison table
# ============================================================
print("=" * 80)
print("SIMULATOR vs HARDWARE: Direct Comparison")
print("=" * 80)
print(f"{'Code':>12s} {'Basis':>5s} | {'Sim Hol':>8s} {'HW Hol':>8s} {'Gap':>8s} | {'Sim Err':>8s} {'HW Err':>8s} {'Ratio':>6s}")
print("-" * 72)

for code in ['uncoded', '3-qubit', '[[5,1,3]]']:
    for basis in ['Z', 'X', 'Y']:
        sh = sim[code][basis]['holevo']
        hh = hw[code][basis]['holevo']
        se = sim[code][basis]['avg_err']
        he = hw[code][basis]['avg_err']
        gap = hh - sh
        ratio = he / se if se > 0.001 else float('inf')
        print(f"{code:>12s} {basis:>5s} | {sh:8.4f} {hh:8.4f} {gap:+8.4f} | {se:8.5f} {he:8.5f} {ratio:6.2f}")
    print()

# ============================================================
# 2. Asymmetry comparison
# ============================================================
print("\n" + "=" * 60)
print("ASYMMETRY: Simulator vs Hardware")
print("=" * 60)
print(f"{'Code':>12s} | {'Sim Asym':>9s} | {'HW Asym':>9s} | {'Ratio':>6s} | Interpretation")
print("-" * 70)

for code in ['uncoded', '3-qubit', '[[5,1,3]]']:
    sa = sim[code]['asymmetry']
    ha = hw[code]['asymmetry']
    ratio = ha / sa if sa > 0.001 else float('inf')

    if code == 'uncoded':
        interp = "HW slightly more asymmetric (qubit-specific noise)"
    elif code == '3-qubit':
        interp = "Nearly identical — structural asymmetry dominates"
    else:
        interp = "4x degradation — correlated noise breaks isotropy"

    print(f"{code:>12s} | {sa:9.4f} | {ha:9.4f} | {ratio:6.1f}x | {interp}")

# ============================================================
# 3. Infer effective noise parameters from hardware data
# ============================================================
print("\n" + "=" * 60)
print("NOISE PARAMETER INFERENCE")
print("=" * 60)

# From uncoded error rates, infer effective depolarizing parameter
# For uncoded with depolarizing noise p: error rate = p * (2/3) for each Pauli direction
# Readout error p_ro: error rate ≈ p_ro (for very low gate noise)
# Total: err ≈ p_gate + p_ro

# From hardware: uncoded errors Z=0.0044, X=0.0031, Y=0.0029
# Average: 0.0035 (this is mostly readout error for depth-1/2 circuits)
hw_uncoded_avg_err = np.mean([hw['uncoded'][b]['avg_err'] for b in ['Z', 'X', 'Y']])
print(f"Uncoded avg error (hardware): {hw_uncoded_avg_err:.5f}")
print(f"  → Dominated by readout error (~{hw_uncoded_avg_err*100:.2f}%)")

# From simulator: uncoded errors all ~0.015 (with p_ro=0.015)
sim_uncoded_avg_err = np.mean([sim['uncoded'][b]['avg_err'] for b in ['Z', 'X', 'Y']])
print(f"Uncoded avg error (simulator): {sim_uncoded_avg_err:.5f}")
print(f"  → Hardware has ~{sim_uncoded_avg_err/hw_uncoded_avg_err:.1f}x lower base error")

# Effective readout error from hardware
p_ro_hw = hw_uncoded_avg_err
print(f"\nInferred hardware readout error: ~{p_ro_hw*100:.2f}%")

# From 3-qubit Z-basis: majority vote + bit-flip correction
# Expected error for 3-qubit under depolarizing: 3p² (for small p)
# HW: 0.0023 (Z-basis), so 3p² ≈ 0.0023 → p ≈ 0.028
# But this also includes readout (3 readouts, majority vote)
# For majority vote with readout error r: P(error) ≈ 3r²(1-r) + r³ ≈ 3r²
# With r=0.0035: 3*(0.0035)² ≈ 0.000037 — much less than 0.0023
# So the 3-qubit Z-error is dominated by physical gate errors, not readout
hw_3q_z_err = hw['3-qubit']['Z']['avg_err']
print(f"\n3-qubit Z-basis error: {hw_3q_z_err:.5f}")
# If error = 3*p_gate²: p_gate ≈ sqrt(0.0023/3) ≈ 0.028
# But need to account for 2 CX gates in encoding: ~2*p_2q
# If majority vote corrects 1 error: residual ≈ 3*p²
# p = total per-qubit error ≈ 2*p_2q + p_ro
# So 3*p² ≈ 0.0023 → p ≈ 0.028 → p_2q ≈ (0.028 - 0.0035)/2 ≈ 0.012

# From [[5,1,3]]: 10 CX gates → ~10*p_2q total gate error
# Without error correction: err ≈ 10*p_2q + 5*p_ro
# HW: ~0.11 → 10*p_2q ≈ 0.11 - 5*0.0035 ≈ 0.093 → p_2q ≈ 0.0093
hw_513_avg_err = np.mean([hw['[[5,1,3]]'][b]['avg_err'] for b in ['Z', 'X', 'Y']])
p_2q_inf = (hw_513_avg_err - 5 * p_ro_hw) / 10
print(f"[[5,1,3]] avg error: {hw_513_avg_err:.5f}")
print(f"  Inferred 2Q gate error: ~{p_2q_inf*100:.2f}% (vs model {sim_data['noise_model']['p2q']*100:.2f}%)")

# ============================================================
# 4. What does the asymmetry degradation tell us?
# ============================================================
print("\n" + "=" * 60)
print("WHAT THE GAP REVEALS")
print("=" * 60)

# [[5,1,3]] asymmetry on hardware: 0.040 vs simulator 0.010
# The 4x degradation means hardware noise isn't purely symmetric depolarizing
# Possible sources:
# 1. Qubit-specific T1/T2 → different qubits experience different noise
# 2. Crosstalk → correlated errors between neighboring qubits
# 3. Coherent errors → systematic rotations, not random
# 4. Readout asymmetry → different readout fidelity for |0⟩ vs |1⟩

print("[[5,1,3]] isotropy degradation: 4x worse on hardware")
print()

# Check per-state asymmetry in [[5,1,3]]
for basis in ['Z', 'X', 'Y']:
    e0 = hw['[[5,1,3]]'][basis]['err_0']
    e1 = hw['[[5,1,3]]'][basis]['err_1']
    print(f"  {basis}-basis: err(|0⟩_L) = {e0:.4f}, err(|1⟩_L) = {e1:.4f}, diff = {abs(e1-e0):.4f}")

print()
print("State asymmetry (|0⟩_L vs |1⟩_L) within each basis:")
print("  If large → readout-dominated or T1-dominated (asymmetric)")
print("  If small → gate-noise dominated (symmetric)")

# Same for 3-qubit
print()
for basis in ['Z', 'X', 'Y']:
    e0 = hw['3-qubit'][basis]['err_0']
    e1 = hw['3-qubit'][basis]['err_1']
    print(f"  3-qubit {basis}-basis: err(|0⟩_L) = {e0:.4f}, err(|1⟩_L) = {e1:.4f}, diff = {abs(e1-e0):.4f}")

# ============================================================
# 5. The key question: would [[5,1,3]] win with lower gate error?
# ============================================================
print("\n" + "=" * 60)
print("BREAK-EVEN ANALYSIS")
print("=" * 60)

# At what 2Q gate error does [[5,1,3]] break even with uncoded?
# Uncoded error ≈ p_ro ≈ 0.0035 (mostly readout)
# [[5,1,3]] error ≈ 10*p_2q + 5*p_ro (encoding overhead)
# Break-even: 10*p_2q + 5*0.0035 = 0.0035
# 10*p_2q = -4*0.0035 = -0.014 → impossible without error correction!

# With error correction (distance 3, corrects 1 error):
# [[5,1,3]] logical error ≈ C(5,2)*p_2q^2 (dominant correctable errors)
# = 10*p_2q^2
# Break-even with uncoded: 10*p_2q^2 = p_ro
# p_2q = sqrt(p_ro / 10) = sqrt(0.0035/10) = 0.019

print(f"Current hardware 2Q gate error: ~{p_2q_inf*100:.2f}%")
print(f"Break-even for [[5,1,3]] (without active correction): impossible")
print(f"  (encoding adds more noise than it protects against)")
print(f"Break-even for [[5,1,3]] (WITH active correction): p_2q < {np.sqrt(p_ro_hw/10)*100:.2f}%")
print(f"  Current hardware: {p_2q_inf*100:.2f}% → {'below' if p_2q_inf < np.sqrt(p_ro_hw/10) else 'above'} threshold")
print()
print(f"3-qubit Z-basis break-even: already winning (0.976 > 0.959 uncoded)")
print(f"  3-qubit gains from majority vote error correction")
print(f"  But ONLY in Z-basis — X/Y basis worse than uncoded")

# ============================================================
# 6. Summary metrics
# ============================================================
print("\n" + "=" * 60)
print("SUMMARY METRICS")
print("=" * 60)

metrics = {
    'isotropy_ratio': round(hw['3-qubit']['asymmetry'] / hw['[[5,1,3]]']['asymmetry'], 1),
    'hw_vs_sim_isotropy_513': round(hw['[[5,1,3]]']['asymmetry'] / sim['[[5,1,3]]']['asymmetry'], 1),
    'hw_vs_sim_isotropy_3q': round(hw['3-qubit']['asymmetry'] / sim['3-qubit']['asymmetry'], 2),
    'best_single_measurement': '3-qubit Z-basis (0.976)',
    'best_average': 'uncoded (0.967)',
    'best_isotropic_code': '[[5,1,3]] (0.040 asymmetry)',
    'inferred_readout_error': round(p_ro_hw * 100, 3),
    'inferred_2q_error': round(p_2q_inf * 100, 3),
}

for k, v in metrics.items():
    print(f"  {k}: {v}")

elapsed = time.time() - start
print(f"\nTotal time: {elapsed:.1f}s")

# Save
output = {
    'sprint': '025c',
    'description': 'Gap analysis between simulator and hardware',
    'metrics': metrics,
    'noise_inference': {
        'readout_error': round(p_ro_hw, 5),
        'two_qubit_gate_error': round(p_2q_inf, 5),
        'break_even_with_correction': round(np.sqrt(p_ro_hw / 10), 5),
    },
    'comparison': {},
    'elapsed_seconds': round(elapsed, 1),
}

for code in ['uncoded', '3-qubit', '[[5,1,3]]']:
    output['comparison'][code] = {
        'sim_avg_holevo': sim[code]['avg_holevo'],
        'hw_avg_holevo': hw[code]['avg_holevo'],
        'sim_asymmetry': sim[code]['asymmetry'],
        'hw_asymmetry': hw[code]['asymmetry'],
        'holevo_gap': round(hw[code]['avg_holevo'] - sim[code]['avg_holevo'], 4),
        'asymmetry_ratio': round(hw[code]['asymmetry'] / max(sim[code]['asymmetry'], 0.001), 2),
    }

with open('results/sprint_025c_gap_analysis.json', 'w') as f:
    json.dump(output, f, indent=2)
print("Saved to results/sprint_025c_gap_analysis.json")
