"""Sprint 062b: Twisted boundary conditions for Potts model.

Insert a Z_q twist at the boundary: replace delta(s_n, s_1) with delta(s_n, (s_1+k) mod q).
This is equivalent to threading a Z_q flux through the ring.

For a Luttinger liquid (compact boson at radius R):
  ΔE(twist k) = E_0(twist k) - E_0(no twist) = (2π v / L) × (k/(qR))²/2

  But k/q is the fraction of the full twist 2π/q^... actually let me be more careful.

  A Z_q twist by 1 unit corresponds to shifting the boundary condition by 2π/q.
  For a compact boson with U(1) symmetry, threading flux φ:
    ΔE(φ) = v·φ²/(4πL) × K
  where K is the Luttinger parameter.

  The Z_q twist by k units: φ = 2πk/q, so:
    ΔE(k) = v·(2πk/q)²/(4πL) × K = v·π·k²·K/(q²·L)

We can extract v from the gap: Δ₁ = 2πv·x₁/L
So ΔE(k)/Δ₁ = k²·K/(2q²·x₁)
→ K = 2q²·x₁·ΔE(k)/(k²·Δ₁)

For Ising (q=2, K=1/2): K = 2·4·0.125·ΔE/(Δ₁) = ΔE/Δ₁

Also: for a Luttinger liquid, K determines x_σ = 1/(2K), so we can check consistency.
"""
import numpy as np
from scipy.sparse import csr_matrix, kron as sp_kron, eye as sp_eye
from scipy.sparse.linalg import eigsh
import json, time

def potts_hamiltonian_twisted(n, q, g, twist=0):
    """Potts Hamiltonian with twisted boundary condition.

    H = -J Σ δ(s_i, s_{i+1}) - g Σ (X + X†)

    The boundary bond (n-1, 0) uses δ(s_{n-1}, (s_0 + twist) mod q).
    twist=0 is periodic BC. twist=k threads k units of Z_q flux.
    """
    dim = q**n
    eye_q = np.eye(q)
    X = np.zeros((q, q))
    for s in range(q):
        X[(s+1) % q, s] = 1.0
    XpXd = X + X.T

    H = csr_matrix((dim, dim))

    # Nearest-neighbor coupling: -Σ δ(s_i, s_{i+1})
    delta_2 = np.zeros((q**2, q**2))
    for s in range(q):
        delta_2[s*q+s, s*q+s] = 1.0

    # Bulk bonds (i, i+1) for i = 0, ..., n-2
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

    # Boundary bond (n-1, 0) with twist: δ(s_{n-1}, (s_0 + twist) mod q)
    delta_twist = np.zeros((q, q))  # acts on (site 0, site n-1) basis
    for s in range(q):
        # s_{n-1} = (s_0 + twist) mod q → s_0 = (s_{n-1} - twist) mod q
        s0 = s  # site 0 state
        sn = (s0 + twist) % q  # site n-1 state
        delta_twist[s0 * 1 + 0, 0] = 0  # will build full projector below

    # Actually build the boundary bond properly
    # |s_0, s_1, ..., s_{n-1}>
    # delta(s_{n-1}, (s_0+twist) mod q) projects onto states where s_{n-1} = (s_0+twist)%q
    # Build as: Σ_s |s, ..., (s+twist)%q><s, ..., (s+twist)%q|  over site 0 and n-1
    for s in range(q):
        sn = (s + twist) % q
        # Projector: |s>_0 <s|_0  ⊗  I_middle  ⊗  |sn>_{n-1} <sn|_{n-1}
        proj_0 = np.zeros((q, q)); proj_0[s, s] = 1.0
        proj_n = np.zeros((q, q)); proj_n[sn, sn] = 1.0

        ops = [eye_q] * n
        ops[0] = proj_0
        ops[n-1] = proj_n
        r = csr_matrix(ops[0])
        for op in ops[1:]:
            r = sp_kron(r, csr_matrix(op), format='csr')
        H = H - r

    # Transverse field: -g Σ (X + X†)
    for i in range(n):
        ops = [eye_q] * n
        ops[i] = XpXd
        r = csr_matrix(ops[0])
        for op in ops[1:]:
            r = sp_kron(r, csr_matrix(op), format='csr')
        H = H - g * r

    return H

def g_c_formula(q):
    return 0.2 * np.sqrt(q - 1) + 0.05

# First: validate on q=2 (Ising) where everything is known
print("=" * 70)
print("TWISTED BOUNDARY CONDITIONS: SPIN STIFFNESS")
print("=" * 70)

print("\n--- Validation: q=2 Ising at g_c=0.250 ---")
print("Ising is a Luttinger liquid with K=1/2, c=1/2")
print("Prediction: ΔE(twist 1) = v·π·K/(q²·L) = v·π/(8L)")
print("With Δ₁ = 2πv·x₁/L, x₁=1/8:")
print("ΔE/Δ₁ = K/(2q²·x₁) = 0.5/(2·4·0.125) = 0.5")
print()

all_results = {}

for q in [2, 3, 4, 5, 7]:
    gc_vals = {2: 0.250, 3: 0.333, 4: 0.392, 5: 0.441, 7: 0.535}
    gc = gc_vals[q]

    print(f"\n{'='*60}")
    print(f"q={q}, g_c={gc}")
    print(f"{'='*60}")

    results_q = {"q": q, "g_c": gc, "sizes": {}}

    for n in [4, 6, 8]:
        dim = q**n
        if dim > 2_000_000:
            print(f"  n={n}: dim={dim} too large, skipping")
            continue

        t0 = time.time()

        # Untwisted (periodic BC)
        H0 = potts_hamiltonian_twisted(n, q, gc, twist=0)
        n_eig = min(6, dim - 2)
        evals0, _ = eigsh(H0, k=n_eig, which='SA')
        evals0 = np.sort(evals0)
        E0 = evals0[0]
        delta1 = evals0[1] - evals0[0]  # First gap

        # Twist by 1 unit
        H1 = potts_hamiltonian_twisted(n, q, gc, twist=1)
        evals1, _ = eigsh(H1, k=2, which='SA')
        E1_twist = np.sort(evals1)[0]
        dE1 = E1_twist - E0

        # Twist by 2 units (if q > 2)
        dE2 = None
        if q > 2:
            H2 = potts_hamiltonian_twisted(n, q, gc, twist=2)
            evals2, _ = eigsh(H2, k=2, which='SA')
            E2_twist = np.sort(evals2)[0]
            dE2 = E2_twist - E0

        t_elapsed = time.time() - t0

        # Extract stiffness: ΔE(k) = v·π·k²·K/(q²·L)
        # ΔE/Δ₁ = k²·K/(2q²·x₁)
        # → K = 2q²·x₁·(ΔE/Δ₁)/k²
        # But we need x₁. Use model-independent ratio ΔE·N/(Δ₁·N) = ΔE/Δ₁
        ratio1 = dE1 / delta1 if delta1 > 1e-15 else None
        ratio2 = (dE2 / delta1) if (dE2 is not None and delta1 > 1e-15) else None

        # Spin stiffness: ρ_s = L · ΔE(k=1) / (2πk/q)² = L · dE1 · q² / (4π²)
        # Actually: ΔE = ρ_s · (2π/q)² / (2L) for lowest twist
        # → ρ_s = ΔE · 2L · q² / (4π²) = ΔE · L · q² / (2π²)
        rho_s_times_L = dE1 * n * q**2 / (2 * np.pi**2) if dE1 is not None else None

        # Check quadratic scaling: ΔE(2)/ΔE(1) should = 4
        ratio_21 = dE2 / dE1 if (dE2 is not None and dE1 > 1e-15) else None

        print(f"\n  n={n} (dim={dim}, {t_elapsed:.1f}s):")
        print(f"    E0 = {E0:.8f}")
        print(f"    Δ₁ = {delta1:.8f}")
        print(f"    ΔE(twist=1) = {dE1:.8f}")
        print(f"    ΔE/Δ₁ = {ratio1:.6f}" if ratio1 else "    ΔE/Δ₁ = N/A")
        print(f"    ρ_s·L = {rho_s_times_L:.6f}" if rho_s_times_L else "    ρ_s·L = N/A")
        if dE2 is not None:
            print(f"    ΔE(twist=2) = {dE2:.8f}")
            print(f"    ΔE(2)/ΔE(1) = {ratio_21:.4f} (boson predicts 4.000)")

        results_q["sizes"][str(n)] = {
            "dim": dim,
            "E0": float(E0),
            "delta1": float(delta1),
            "dE_twist1": float(dE1),
            "dE_twist2": float(dE2) if dE2 is not None else None,
            "ratio_dE_delta": float(ratio1) if ratio1 else None,
            "ratio_dE2_dE1": float(ratio_21) if ratio_21 else None,
            "rho_s_L": float(rho_s_times_L) if rho_s_times_L else None,
            "time_s": round(t_elapsed, 2),
        }

    all_results[f"q={q}"] = results_q

# Summary
print("\n\n" + "=" * 70)
print("SUMMARY: TWISTED BC RESULTS")
print("=" * 70)

print(f"\n{'q':>3} {'n':>3} {'ΔE(1)/Δ₁':>10} {'ΔE(2)/ΔE(1)':>12} {'ρ_s·L':>10}")
for key in sorted(all_results.keys(), key=lambda x: all_results[x]["q"]):
    r = all_results[key]
    q = r["q"]
    for ns, sd in sorted(r["sizes"].items(), key=lambda x: int(x[0])):
        n = int(ns)
        ratio = sd["ratio_dE_delta"]
        r21 = sd.get("ratio_dE2_dE1")
        rho = sd.get("rho_s_L")
        r_str = f"{ratio:.6f}" if ratio else "N/A"
        r21_str = f"{r21:.4f}" if r21 else "—"
        rho_str = f"{rho:.6f}" if rho else "N/A"
        print(f"{q:>3} {n:>3} {r_str:>10} {r21_str:>12} {rho_str:>10}")

# Check if ΔE(2)/ΔE(1) = 4 (quadratic, free boson)
print("\n--- Quadratic twist test: ΔE(k) ∝ k² ---")
print("Free boson/Luttinger liquid predicts ΔE(2)/ΔE(1) = 4")
for key in sorted(all_results.keys(), key=lambda x: all_results[x]["q"]):
    r = all_results[key]
    q = r["q"]
    for ns, sd in sorted(r["sizes"].items(), key=lambda x: int(x[0])):
        if sd.get("ratio_dE2_dE1") is not None:
            dev = abs(sd["ratio_dE2_dE1"] - 4) / 4 * 100
            print(f"  q={q}, n={int(ns)}: ΔE(2)/ΔE(1) = {sd['ratio_dE2_dE1']:.4f} ({dev:.1f}% from 4)")

# Check ρ_s·L convergence (should be constant for Luttinger liquid)
print("\n--- ρ_s·L convergence (constant = Luttinger liquid) ---")
for key in sorted(all_results.keys(), key=lambda x: all_results[x]["q"]):
    r = all_results[key]
    q = r["q"]
    vals = [(int(ns), sd["rho_s_L"]) for ns, sd in r["sizes"].items() if sd.get("rho_s_L")]
    vals.sort()
    if len(vals) >= 2:
        ns = [v[0] for v in vals]
        rhos = [v[1] for v in vals]
        drift = abs(rhos[-1] - rhos[0]) / abs(rhos[0]) * 100
        print(f"  q={q}: ρ_s·L = {', '.join(f'{v:.4f}' for v in rhos)} (drift {drift:.1f}%)")
    elif vals:
        print(f"  q={q}: ρ_s·L = {vals[0][1]:.4f} (single size)")

with open("results/sprint_062b_twisted_bc.json", "w") as f:
    json.dump(all_results, f, indent=2)
print("\nResults saved to results/sprint_062b_twisted_bc.json")
