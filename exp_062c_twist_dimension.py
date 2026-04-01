"""Sprint 062c: Test twist-energy = scaling-dimension correspondence.

Key prediction: In a CFT, twisted BC ground state energy is:
  ΔE(twist k) = (2π/L) × x_twist_k

So ΔE(k)/ΔE(1) should equal x_{twist_k}/x_{twist_1}.

If the twist field IS σᵏ, then ΔE(k)/ΔE(1) = x_σᵏ/x_σ = R(σᵏ)/R(σ)
which we already measured in Sprints 057, 061.

Also test:
1. Does ΔE(1) × L/(2π) → x₁ as L→∞? (extract x_twist directly)
2. q=10 at n=4 with all twist sectors k=1..5
3. ρ_s·L trend: extrapolate to L→∞
"""
import numpy as np
from scipy.sparse import csr_matrix, kron as sp_kron, eye as sp_eye
from scipy.sparse.linalg import eigsh
import json, time

def potts_hamiltonian_twisted(n, q, g, twist=0):
    dim = q**n
    eye_q = np.eye(q)
    X = np.zeros((q, q))
    for s in range(q):
        X[(s+1) % q, s] = 1.0
    XpXd = X + X.T
    H = csr_matrix((dim, dim))

    # Bulk bonds
    delta_2 = np.zeros((q**2, q**2))
    for s in range(q):
        delta_2[s*q+s, s*q+s] = 1.0
    for i in range(n - 1):
        j = i + 1
        left = q**i if i > 0 else 1
        right = q**(n-j-1) if n-j-1 > 0 else 1
        op = csr_matrix(delta_2)
        if i > 0:
            op = sp_kron(sp_eye(left), op, format='csr')
        if n-j-1 > 0:
            op = sp_kron(op, sp_eye(right), format='csr')
        H = H - op

    # Boundary bond with twist
    for s in range(q):
        sn = (s + twist) % q
        proj_0 = np.zeros((q, q)); proj_0[s, s] = 1.0
        proj_n = np.zeros((q, q)); proj_n[sn, sn] = 1.0
        ops = [eye_q] * n; ops[0] = proj_0; ops[n-1] = proj_n
        r = csr_matrix(ops[0])
        for op in ops[1:]:
            r = sp_kron(r, csr_matrix(op), format='csr')
        H = H - r

    # Transverse field
    for i in range(n):
        ops = [eye_q] * n; ops[i] = XpXd
        r = csr_matrix(ops[0])
        for op in ops[1:]:
            r = sp_kron(r, csr_matrix(op), format='csr')
        H = H - g * r
    return H

def g_c_formula(q):
    return 0.2 * np.sqrt(q - 1) + 0.05

# Known harmonic ratios from Sprint 057/061 (periodic BC spectrum)
harmonic_ratios = {
    # R(σ²)/R(σ) from charge-resolved spectrum
    4: {2: 1.68},     # σ² at R=1.68
    5: {2: 2.41},     # σ² at R=2.41
    7: {2: 3.32, 3: 5.95},  # σ², σ³
    10: {2: 3.66, 3: 7.69}, # σ², σ³
}

print("=" * 70)
print("TWIST-ENERGY = SCALING DIMENSION TEST")
print("=" * 70)

gc_vals = {2: 0.250, 3: 0.333, 4: 0.392, 5: 0.441, 7: 0.535, 10: 0.684}
all_results = {}

# Part 1: Extract x_twist from ΔE·L/(2π) and compare to x₁ from spectrum
print("\n--- Part 1: x_twist = ΔE·L/(2πv) ---")
print("Need velocity v from CFT: v = Δ₁·L/(2π·x₁)")
print("Known x₁: q=2→0.125, q=3→0.133, q=4→0.117, q=5→0.101, q=7→0.086")

x1_known = {2: 0.125, 3: 0.1337, 4: 0.1172, 5: 0.1015, 7: 0.0860, 10: 0.0831}

for q in [2, 3, 4, 5, 7, 10]:
    gc = gc_vals[q]
    max_twists = min(q // 2, 5)

    print(f"\n{'='*60}")
    print(f"q={q}, g_c={gc}")
    print(f"{'='*60}")

    results_q = {"q": q, "g_c": gc, "sizes": {}}

    for n in [4, 6, 8]:
        dim = q**n
        if dim > 2_000_000:
            continue

        t0 = time.time()

        # Periodic BC ground state + gap
        H0 = potts_hamiltonian_twisted(n, q, gc, twist=0)
        n_eig = min(6, dim - 2)
        evals0, _ = eigsh(H0, k=n_eig, which='SA')
        evals0 = np.sort(evals0)
        E0_periodic = evals0[0]
        delta1 = evals0[1] - evals0[0]

        # Velocity from gap
        x1 = x1_known[q]
        v = delta1 * n / (2 * np.pi * x1)

        # All twist sectors
        twist_data = {}
        for k in range(1, max_twists + 1):
            Hk = potts_hamiltonian_twisted(n, q, gc, twist=k)
            evalsk, _ = eigsh(Hk, k=2, which='SA')
            E0_twist = np.sort(evalsk)[0]
            dE = E0_twist - E0_periodic

            # Extract x_twist = ΔE · L / (2π v)
            x_twist = dE * n / (2 * np.pi * v)

            # Also: ΔE·n is scale-invariant at criticality
            dE_n = dE * n

            twist_data[k] = {
                "dE": float(dE),
                "dE_n": float(dE_n),
                "x_twist": float(x_twist),
            }

        t_elapsed = time.time() - t0

        # Twist ratios
        dE_1 = twist_data[1]["dE"]
        print(f"\n  n={n} (dim={dim}, {t_elapsed:.1f}s):")
        print(f"    Δ₁ = {delta1:.8f}, v = {v:.6f}")
        print(f"    {'k':>3} {'ΔE(k)':>12} {'ΔE(k)/ΔE(1)':>12} {'x_twist':>10} {'x₁·k²':>10}")
        for k, td in sorted(twist_data.items()):
            ratio = td["dE"] / dE_1 if dE_1 > 1e-15 else None
            x_boson = x1 * k**2  # free boson prediction
            r_str = f"{ratio:.4f}" if ratio else "N/A"
            print(f"    {k:>3} {td['dE']:>12.8f} {r_str:>12} {td['x_twist']:>10.5f} {x_boson:>10.5f}")

        # Compare to harmonic ratios from spectrum
        if q in harmonic_ratios:
            print(f"    Spectrum harmonic ratios (Sprint 057/061):")
            for k, R_ratio in harmonic_ratios[q].items():
                twist_ratio = twist_data[k]["dE"] / dE_1 if k in twist_data and dE_1 > 1e-15 else None
                if twist_ratio:
                    diff = abs(twist_ratio - R_ratio) / R_ratio * 100
                    print(f"      σ^{k}: twist={twist_ratio:.4f}, spectrum={R_ratio:.2f}, diff={diff:.1f}%")

        results_q["sizes"][str(n)] = {
            "dim": dim, "delta1": float(delta1), "v": float(v),
            "twists": twist_data, "time_s": round(t_elapsed, 2),
        }

    all_results[f"q={q}"] = results_q

# Part 2: q=10 all twist sectors at n=4
print("\n" + "=" * 70)
print("Part 2: q=10 complete twist analysis at n=4")
print("=" * 70)

q = 10; gc = 0.684; n = 4; dim = 10000
H0 = potts_hamiltonian_twisted(n, q, gc, twist=0)
evals0, _ = eigsh(H0, k=6, which='SA')
evals0 = np.sort(evals0)
E0 = evals0[0]; delta1 = evals0[1] - evals0[0]

print(f"\nq=10, n=4, dim={dim}")
print(f"E0={E0:.8f}, Δ₁={delta1:.8f}")
print(f"\n{'k':>3} {'ΔE(k)':>12} {'ΔE(k)/ΔE(1)':>12} {'ΔE(k)/k²/ΔE(1)':>16}")

dE_all = {}
for k in range(1, 6):  # k=1 to 5 (k and 10-k are conjugate)
    Hk = potts_hamiltonian_twisted(n, q, gc, twist=k)
    evalsk, _ = eigsh(Hk, k=2, which='SA')
    E0k = np.sort(evalsk)[0]
    dE = E0k - E0
    dE_all[k] = dE

dE_1 = dE_all[1]
for k in range(1, 6):
    ratio = dE_all[k] / dE_1
    ratio_k2 = ratio / k**2
    print(f"{k:>3} {dE_all[k]:>12.8f} {ratio:>12.4f} {ratio_k2:>16.4f}")

# k and q-k should give same energy (conjugate twists)
print("\nConjugate twist test (k ↔ q-k):")
for k in range(1, 6):
    Hconj = potts_hamiltonian_twisted(n, q, gc, twist=q-k)
    econj, _ = eigsh(Hconj, k=2, which='SA')
    E0conj = np.sort(econj)[0]
    dE_conj = E0conj - E0
    print(f"  twist {k}: {dE_all[k]:.8f}  |  twist {q-k}: {dE_conj:.8f}  |  diff: {abs(dE_all[k]-dE_conj):.2e}")

# Part 3: Summary comparison - twist vs spectrum
print("\n" + "=" * 70)
print("SUMMARY: Twist ΔE(k)/ΔE(1) vs Spectrum R(σᵏ)/R(σ)")
print("=" * 70)

print(f"\n{'q':>3} {'k':>3} {'n':>3} {'Twist ratio':>12} {'Spectrum R':>11} {'Diff(%)':>8}")
for q_key in sorted(all_results.keys(), key=lambda x: all_results[x]["q"]):
    q = all_results[q_key]["q"]
    if q not in harmonic_ratios:
        continue
    # Use largest n available
    sizes = all_results[q_key]["sizes"]
    largest_n = max(sizes.keys(), key=int)
    sd = sizes[largest_n]
    dE_1 = sd["twists"][1]["dE"] if 1 in sd["twists"] else sd["twists"]["1"]["dE"]
    for k_key, td in sorted(sd["twists"].items(), key=lambda x: int(x) if isinstance(x, str) else x):
        k = int(k_key) if isinstance(k_key, str) else k_key
        if k in harmonic_ratios[q]:
            twist_ratio = td["dE"] / dE_1
            spec_ratio = harmonic_ratios[q][k]
            diff = (twist_ratio - spec_ratio) / spec_ratio * 100
            print(f"{q:>3} {k:>3} {largest_n:>3} {twist_ratio:>12.4f} {spec_ratio:>11.2f} {diff:>8.1f}")

# Part 4: c·x₁ product analysis (from 062a finding)
print("\n" + "=" * 70)
print("BONUS: c·x₁ ≈ constant for q≥3")
print("=" * 70)
c_data = {2: 0.500, 3: 0.800, 4: 1.000, 5: 1.10, 7: 1.30, 10: 1.40}
print(f"{'q':>3} {'c':>6} {'x₁':>7} {'c·x₁':>7} {'1/9':>7} {'diff':>7}")
for q in sorted(c_data.keys()):
    c = c_data[q]
    x1 = x1_known[q]
    cx1 = c * x1
    ninth = 1/9
    diff = (cx1 - ninth) / ninth * 100
    print(f"{q:>3} {c:>6.3f} {x1:>7.4f} {cx1:>7.5f} {ninth:>7.5f} {diff:>6.1f}%")

with open("results/sprint_062c_twist_dimension.json", "w") as f:
    json.dump(all_results, f, indent=2)
print("\nResults saved to results/sprint_062c_twist_dimension.json")
