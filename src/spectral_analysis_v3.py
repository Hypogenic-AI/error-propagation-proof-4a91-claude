"""
Spectral analysis v3: Lazy random walks (aperiodic) and Diaconis upper bound.

Key insight from v2: periodicity confounds abelian/non-abelian comparison.
Fix: use lazy walks P' = (1/2)I + (1/2)P to ensure aperiodicity.

Main theoretical result: Diaconis's upper bound lemma gives
  ‖μ^{*k} - U‖² ≤ (1/4) Σ_{ρ nontrivial} d_ρ · ‖ρ̂(μ)^k‖_F²

For abelian groups: all d_ρ = 1 → tighter bound
For non-abelian groups: some d_ρ > 1 → potentially looser bound
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import json
import os
from itertools import permutations

np.random.seed(42)

# ============================================================================
# Group Constructions
# ============================================================================

def cyclic_group_cayley(n):
    return np.array([[(i + j) % n for j in range(n)] for i in range(n)])

def dihedral_group_cayley(n):
    order = 2 * n
    table = np.zeros((order, order), dtype=int)
    for i in range(order):
        for j in range(order):
            if i < n and j < n:
                table[i, j] = (i + j) % n
            elif i < n and j >= n:
                table[i, j] = n + (j - n - i) % n
            elif i >= n and j < n:
                table[i, j] = n + (i - n + j) % n
            else:
                table[i, j] = (i - n - (j - n)) % n
    return table

def symmetric_group_cayley(n):
    elems = list(permutations(range(n)))
    elem_to_idx = {e: i for i, e in enumerate(elems)}
    order = len(elems)
    table = np.zeros((order, order), dtype=int)
    for i, g in enumerate(elems):
        for j, h in enumerate(elems):
            gh = tuple(g[h[k]] for k in range(n))
            table[i, j] = elem_to_idx[gh]
    return table, elems, elem_to_idx

def quaternion_group_cayley():
    return np.array([
        [0,1,2,3,4,5,6,7],[1,0,3,2,5,4,7,6],[2,3,1,0,6,7,5,4],[3,2,0,1,7,6,4,5],
        [4,5,7,6,1,0,2,3],[5,4,6,7,0,1,3,2],[6,7,4,5,3,2,1,0],[7,6,5,4,2,3,0,1],
    ])

def klein_four_cayley():
    """Z_2 x Z_2 (abelian, order 4)."""
    return np.array([[0,1,2,3],[1,0,3,2],[2,3,0,1],[3,2,1,0]])

def direct_product_cayley(t1, t2):
    n1, n2 = len(t1), len(t2)
    order = n1 * n2
    table = np.zeros((order, order), dtype=int)
    for i in range(order):
        for j in range(order):
            table[i, j] = t1[i//n2, j//n2] * n2 + t2[i%n2, j%n2]
    return table

# ============================================================================
# Operator and Analysis Tools
# ============================================================================

def random_walk_operator(cayley_table, generators, weights=None):
    n = len(cayley_table)
    if weights is None:
        weights = np.ones(len(generators)) / len(generators)
    P = np.zeros((n, n))
    for g_idx, w in zip(generators, weights):
        for i in range(n):
            P[i, cayley_table[i, g_idx]] += w
    return P

def lazy_walk(P, alpha=0.5):
    """Lazy random walk: P' = alpha*I + (1-alpha)*P. Ensures aperiodicity."""
    n = len(P)
    return alpha * np.eye(n) + (1 - alpha) * P

def spectral_radius_nontrivial(P):
    evals = np.abs(np.linalg.eigvals(P))
    evals_sorted = sorted(evals, reverse=True)
    return evals_sorted[1] if len(evals_sorted) > 1 else 0.0

def is_abelian(t):
    n = len(t)
    return all(t[i, j] == t[j, i] for i in range(n) for j in range(i+1, n))

def compute_inverses(t):
    n = len(t)
    inv = np.zeros(n, dtype=int)
    for i in range(n):
        for j in range(n):
            if t[i, j] == 0:
                inv[i] = j
                break
    return inv

def commutator_subgroup_size(t):
    n = len(t)
    inv = compute_inverses(t)
    comms = set()
    for g in range(n):
        for h in range(n):
            comms.add(t[t[g, h], t[inv[g], inv[h]]])
    subgroup = set(comms)
    changed = True
    while changed:
        changed = False
        new = set()
        for a in subgroup:
            for b in subgroup:
                p = t[a, b]
                if p not in subgroup:
                    new.add(p)
                    changed = True
        subgroup |= new
    return len(subgroup)

def total_variation_distance(P_n, n_group):
    """TV distance of P^n rows from uniform."""
    uniform = np.ones(n_group) / n_group
    return np.max(0.5 * np.sum(np.abs(P_n - uniform), axis=1))

# ============================================================================
# Diaconis Upper Bound Computation
# ============================================================================

def diaconis_upper_bound_eigenvalues(P, k):
    """
    Compute Diaconis upper bound on TV distance at step k.

    For the regular representation, the bound is:
    ‖μ^{*k} - U‖_TV² ≤ (1/4) Σ_{λ ≠ 1} m_λ · |λ|^{2k}

    where m_λ is the multiplicity of eigenvalue λ.
    For abelian groups: each m_λ = 1
    For non-abelian groups: some m_λ > 1 (from higher-dim irreps)
    """
    evals = np.linalg.eigvals(P)
    n = len(P)

    # Group eigenvalues by approximate value
    bound = 0.0
    for ev in evals:
        if abs(abs(ev) - 1.0) > 1e-10 or abs(ev.imag) > 1e-10 or ev.real < 0.99:
            # Non-trivial eigenvalue
            bound += abs(ev) ** (2 * k)

    return 0.25 * bound

# ============================================================================
# Main Experiments
# ============================================================================

def experiment_1_lazy_walk_comparison():
    """
    Compare LAZY random walks on abelian vs non-abelian groups.
    Lazy walks guarantee aperiodicity, isolating the group structure effect.
    """
    print("=" * 75)
    print("EXPERIMENT 1: Lazy Walk Spectral Comparison")
    print("=" * 75)

    groups = []

    # Order 6
    z6 = cyclic_group_cayley(6)
    d3 = dihedral_group_cayley(3)
    for name, t, gens in [
        ("Z_6", z6, [1, 5]),
        ("D_3=S_3", d3, [1, 3]),  # {r, s}
        ("D_3(r,r⁻¹,s)", d3, [1, 2, 3]),
    ]:
        P = lazy_walk(random_walk_operator(t, gens))
        groups.append({"name": name, "order": len(t), "abelian": is_abelian(t),
                       "n_gens": len(gens), "rho": spectral_radius_nontrivial(P),
                       "gap": 1-spectral_radius_nontrivial(P)})

    # Order 8
    z8 = cyclic_group_cayley(8)
    d4 = dihedral_group_cayley(4)
    q8 = quaternion_group_cayley()
    z2z2z2 = direct_product_cayley(direct_product_cayley(cyclic_group_cayley(2), cyclic_group_cayley(2)), cyclic_group_cayley(2))
    for name, t, gens in [
        ("Z_8", z8, [1, 7]),
        ("Z2³", z2z2z2, [1, 2, 4]),  # three generators
        ("D_4", d4, [1, 3, 4]),
        ("Q_8", q8, [2, 3, 4, 5]),
    ]:
        P = lazy_walk(random_walk_operator(t, gens))
        groups.append({"name": name, "order": len(t), "abelian": is_abelian(t),
                       "n_gens": len(gens), "rho": spectral_radius_nontrivial(P),
                       "gap": 1-spectral_radius_nontrivial(P)})

    # Order 24
    z24 = cyclic_group_cayley(24)
    s4_table, s4_elems, s4_idx = symmetric_group_cayley(4)
    s4_gens = []
    identity = list(range(4))
    for i in range(3):
        perm = list(identity)
        perm[i], perm[i+1] = perm[i+1], perm[i]
        s4_gens.append(s4_idx[tuple(perm)])
    for name, t, gens in [
        ("Z_24", z24, [1, 23]),
        ("S_4", s4_table, s4_gens),
    ]:
        P = lazy_walk(random_walk_operator(t, gens))
        groups.append({"name": name, "order": len(t), "abelian": is_abelian(t),
                       "n_gens": len(gens), "rho": spectral_radius_nontrivial(P),
                       "gap": 1-spectral_radius_nontrivial(P)})

    print(f"\n{'Group':<15} {'|G|':>4} {'Abel':>5} {'#Gen':>5} {'ρ(P_lazy)':>10} {'Gap':>10}")
    print("-" * 55)
    for g in groups:
        print(f"{g['name']:<15} {g['order']:>4} {'Y' if g['abelian'] else 'N':>5} "
              f"{g['n_gens']:>5} {g['rho']:>10.6f} {g['gap']:>10.6f}")

    return groups

def experiment_2_convergence_curves():
    """
    Plot convergence curves for lazy walks on matched abelian/non-abelian pairs.
    Also plot Diaconis upper bound.
    """
    print("\n" + "=" * 75)
    print("EXPERIMENT 2: Convergence Curves (Lazy Walks)")
    print("=" * 75)

    n_steps = 80
    pairs = [
        ("Z_6", cyclic_group_cayley(6), [1, 5],
         "D_3", dihedral_group_cayley(3), [1, 2, 3]),
        ("Z_8", cyclic_group_cayley(8), [1, 7],
         "D_4", dihedral_group_cayley(4), [1, 3, 4]),
        ("Z_24", cyclic_group_cayley(24), [1, 23],
         "S_4", symmetric_group_cayley(4)[0],
         []),  # will fill in
    ]

    # Fix S_4 generators
    s4_table, s4_elems, s4_idx = symmetric_group_cayley(4)
    s4_gens = []
    for i in range(3):
        perm = list(range(4))
        perm[i], perm[i+1] = perm[i+1], perm[i]
        s4_gens.append(s4_idx[tuple(perm)])
    pairs[2] = ("Z_24", cyclic_group_cayley(24), [1, 23], "S_4", s4_table, s4_gens)

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    for idx, (a_name, a_table, a_gens, na_name, na_table, na_gens) in enumerate(pairs):
        Pa = lazy_walk(random_walk_operator(a_table, a_gens))
        Pna = lazy_walk(random_walk_operator(na_table, na_gens))
        n_a, n_na = len(a_table), len(na_table)

        # Actual convergence
        tv_a, tv_na = [], []
        bound_a, bound_na = [], []
        Pan, Pnan = np.eye(n_a), np.eye(n_na)

        for step in range(n_steps):
            Pan = Pan @ Pa
            Pnan = Pnan @ Pna
            tv_a.append(total_variation_distance(Pan, n_a))
            tv_na.append(total_variation_distance(Pnan, n_na))
            bound_a.append(np.sqrt(diaconis_upper_bound_eigenvalues(Pa, step+1)))
            bound_na.append(np.sqrt(diaconis_upper_bound_eigenvalues(Pna, step+1)))

        ax = axes[idx]
        steps = range(1, n_steps+1)
        ax.semilogy(steps, tv_a, 'b-', linewidth=2, label=f'{a_name} (TV)')
        ax.semilogy(steps, tv_na, 'r-', linewidth=2, label=f'{na_name} (TV)')
        ax.semilogy(steps, bound_a, 'b--', linewidth=1, alpha=0.6, label=f'{a_name} (UB)')
        ax.semilogy(steps, bound_na, 'r--', linewidth=1, alpha=0.6, label=f'{na_name} (UB)')
        ax.set_xlabel("Step k")
        ax.set_ylabel("Total Variation Distance")
        ax.set_title(f"Order {n_a}: {a_name} vs {na_name}")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
        ax.set_ylim(bottom=1e-8)

        # Print comparison
        rho_a = spectral_radius_nontrivial(Pa)
        rho_na = spectral_radius_nontrivial(Pna)
        print(f"\n{a_name} vs {na_name} (order {n_a}):")
        print(f"  ρ({a_name}) = {rho_a:.6f}, ρ({na_name}) = {rho_na:.6f}")
        print(f"  TV@20: {a_name}={tv_a[19]:.6e}, {na_name}={tv_na[19]:.6e}")
        print(f"  TV@50: {a_name}={tv_a[49]:.6e}, {na_name}={tv_na[49]:.6e}")

    plt.tight_layout()
    plt.savefig("results/plots/lazy_walk_convergence.png", dpi=150)
    plt.close()
    print("\nSaved: results/plots/lazy_walk_convergence.png")

def experiment_3_systematic_comparison():
    """
    Systematic comparison across many group sizes with matched generator counts.
    Use lazy walks with exactly k generators (+ their inverses for symmetry).
    """
    print("\n" + "=" * 75)
    print("EXPERIMENT 3: Systematic Lazy Walk Comparison")
    print("=" * 75)

    results = []

    for n in range(3, 13):
        order = 2 * n

        # Z_{2n} with {1, 2n-1} (generator + inverse), lazy
        z = cyclic_group_cayley(order)
        Pz = lazy_walk(random_walk_operator(z, [1, order-1]))
        rho_z = spectral_radius_nontrivial(Pz)

        # D_n with {r, r^{-1}, s} (3 generators), lazy
        d = dihedral_group_cayley(n)
        Pd = lazy_walk(random_walk_operator(d, [1, n-1, n]))
        rho_d = spectral_radius_nontrivial(Pd)

        comm = commutator_subgroup_size(d)

        results.append({
            "n": n, "order": order,
            "rho_Z": float(rho_z), "gap_Z": float(1-rho_z),
            "rho_D": float(rho_d), "gap_D": float(1-rho_d),
            "comm_size": comm,
            "comm_ratio": comm/order,
            "D_faster": rho_d < rho_z,
        })

    print(f"\n{'n':>3} {'|G|':>5} {'ρ(Z_{2n})':>10} {'ρ(D_n)':>10} {'gap_Z':>8} {'gap_D':>8} "
          f"{'|[D,D]|':>7} {'Faster':>8}")
    print("-" * 70)
    for r in results:
        faster = "D_n" if r["D_faster"] else "Z_{2n}"
        print(f"{r['n']:>3} {r['order']:>5} {r['rho_Z']:>10.6f} {r['rho_D']:>10.6f} "
              f"{r['gap_Z']:>8.4f} {r['gap_D']:>8.4f} {r['comm_size']:>7} {faster:>8}")

    n_d_faster = sum(1 for r in results if r["D_faster"])
    print(f"\nD_n (non-abelian) has smaller ρ in {n_d_faster}/{len(results)} cases")
    print("→ Non-abelian groups with MORE generators can converge FASTER")
    print("→ The key factor is the SPECTRAL GAP, which depends on both")
    print("  group structure AND generator set")

    return results

def experiment_4_fixed_generators_fourier():
    """
    Control for generator count by using exactly the same number of generators.
    For order 2n: Z_{2n} with {1, 2n-1}, D_n with {r, r⁻¹} (both 2-gen, lazy).
    """
    print("\n" + "=" * 75)
    print("EXPERIMENT 4: Fixed 2-Generator Lazy Walk Comparison")
    print("=" * 75)

    results = []

    for n in range(3, 16):
        order = 2 * n

        # Z_{2n}: {1, 2n-1}
        z = cyclic_group_cayley(order)
        Pz = lazy_walk(random_walk_operator(z, [1, order-1]))
        rho_z = spectral_radius_nontrivial(Pz)

        # D_n: {r, r⁻¹} = {1, n-1} (stays in rotation subgroup ≅ Z_n, so walk
        # is reducible!) Instead use {r, s} = {1, n}
        d = dihedral_group_cayley(n)
        Pd = lazy_walk(random_walk_operator(d, [1, n]))
        rho_d = spectral_radius_nontrivial(Pd)

        comm = commutator_subgroup_size(d)

        results.append({
            "n": n, "order": order,
            "rho_Z": float(rho_z), "rho_D": float(rho_d),
            "gap_Z": float(1-rho_z), "gap_D": float(1-rho_d),
            "gap_ratio": float((1-rho_d)/(1-rho_z)) if rho_z < 1 else float('inf'),
            "comm_size": comm,
        })

    print(f"\n{'n':>3} {'|G|':>5} {'ρ(Z_{2n})':>10} {'ρ(D_n)':>10} {'gap_Z':>8} "
          f"{'gap_D':>8} {'gap_D/gap_Z':>11} {'|[D,D]|':>7}")
    print("-" * 75)
    for r in results:
        print(f"{r['n']:>3} {r['order']:>5} {r['rho_Z']:>10.6f} {r['rho_D']:>10.6f} "
              f"{r['gap_Z']:>8.4f} {r['gap_D']:>8.4f} {r['gap_ratio']:>11.4f} "
              f"{r['comm_size']:>7}")

    # Analyze the pattern
    print("\n--- Analysis ---")
    gap_ratios = [r["gap_ratio"] for r in results]
    print(f"Mean gap_D/gap_Z: {np.mean(gap_ratios):.4f}")
    print(f"Min gap_D/gap_Z:  {np.min(gap_ratios):.4f}")
    print(f"Max gap_D/gap_Z:  {np.max(gap_ratios):.4f}")

    # Plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    ns = [r["n"] for r in results]
    gaps_z = [r["gap_Z"] for r in results]
    gaps_d = [r["gap_D"] for r in results]

    ax1.plot(ns, gaps_z, 'bo-', linewidth=2, markersize=6, label='Z_{2n} (abelian)')
    ax1.plot(ns, gaps_d, 'rs-', linewidth=2, markersize=6, label='D_n (non-abelian)')
    ax1.set_xlabel("n")
    ax1.set_ylabel("Spectral Gap")
    ax1.set_title("Spectral Gap: Z_{2n} vs D_n\n(2-generator lazy walk)")
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    ax2.plot(ns, gap_ratios, 'g^-', linewidth=2, markersize=8)
    ax2.axhline(y=1.0, color='k', linestyle='--', alpha=0.5)
    ax2.set_xlabel("n")
    ax2.set_ylabel("gap(D_n) / gap(Z_{2n})")
    ax2.set_title("Gap Ratio (>1 means non-abelian converges faster)")
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig("results/plots/spectral_gap_comparison.png", dpi=150)
    plt.close()
    print("Saved: results/plots/spectral_gap_comparison.png")

    return results

def experiment_5_diaconis_bound_analysis():
    """
    Key theorem verification: The Diaconis upper bound lemma
    gives a bound that depends on irrep dimensions.

    For abelian G: bound = (1/4) Σ_{χ} |χ̂(μ)|^{2k}
    For non-abelian G: bound = (1/4) Σ_{ρ} d_ρ ‖ρ̂(μ)‖_F^{2k}

    The factor d_ρ > 1 for non-abelian groups makes the bound LOOSER,
    even when the spectral radius is the same.
    """
    print("\n" + "=" * 75)
    print("EXPERIMENT 5: Diaconis Upper Bound Analysis")
    print("=" * 75)

    # Compare Z_6 vs D_3 (S_3) with 2-generator lazy walk
    z6 = cyclic_group_cayley(6)
    d3 = dihedral_group_cayley(3)

    Pz = lazy_walk(random_walk_operator(z6, [1, 5]))
    Pd = lazy_walk(random_walk_operator(d3, [1, 3]))

    # Eigenvalue analysis
    evals_z = np.linalg.eigvals(Pz)
    evals_d = np.linalg.eigvals(Pd)

    print("\nZ_6 eigenvalues:")
    for i, ev in enumerate(sorted(evals_z, key=lambda x: -abs(x))):
        print(f"  λ_{i} = {ev.real:+.6f} + {ev.imag:+.6f}i  (|λ| = {abs(ev):.6f})")

    print("\nD_3 eigenvalues:")
    for i, ev in enumerate(sorted(evals_d, key=lambda x: -abs(x))):
        print(f"  λ_{i} = {ev.real:+.6f} + {ev.imag:+.6f}i  (|λ| = {abs(ev):.6f})")

    # Compute bounds vs actual TV
    steps = range(1, 61)
    tv_z, tv_d, ub_z, ub_d = [], [], [], []
    Pzn, Pdn = np.eye(6), np.eye(6)
    for k in steps:
        Pzn = Pzn @ Pz
        Pdn = Pdn @ Pd
        tv_z.append(total_variation_distance(Pzn, 6))
        tv_d.append(total_variation_distance(Pdn, 6))
        ub_z.append(np.sqrt(diaconis_upper_bound_eigenvalues(Pz, k)))
        ub_d.append(np.sqrt(diaconis_upper_bound_eigenvalues(Pd, k)))

    print(f"\nStep  TV(Z_6)     TV(D_3)     UB(Z_6)     UB(D_3)")
    print("-" * 55)
    for k in [5, 10, 15, 20, 30, 40, 50]:
        print(f"{k:>4}  {tv_z[k-1]:.6e}  {tv_d[k-1]:.6e}  "
              f"{ub_z[k-1]:.6e}  {ub_d[k-1]:.6e}")

    return {"tv_z": tv_z, "tv_d": tv_d, "ub_z": ub_z, "ub_d": ub_d}

def experiment_6_proof_model():
    """
    Refined proof refinement model.

    Model proof state as vector in R^d. At each step:
    1. Apply a "correction operator" C that moves toward origin
    2. Apply an "error operator" E drawn from an error group

    Abelian model: E = diag(1 + ε_i) — errors along each coordinate are independent
    Non-abelian model: E = rotation by small angle — errors couple coordinates

    Key prediction: Non-abelian errors cause coordinate coupling that can
    amplify or dampen errors depending on alignment.
    """
    print("\n" + "=" * 75)
    print("EXPERIMENT 6: Proof Refinement with Structured Error Models")
    print("=" * 75)

    d = 6  # proof state dimension
    n_steps = 200
    n_trials = 2000
    correction_rate = 0.15  # how much we move toward origin each step

    error_magnitudes = [0.01, 0.03, 0.05, 0.08, 0.10, 0.12]

    plt.figure(figsize=(12, 6))

    all_results = {}

    for eps in error_magnitudes:
        dists_abelian = np.zeros((n_trials, n_steps + 1))
        dists_nonabelian = np.zeros((n_trials, n_steps + 1))

        for trial in range(n_trials):
            x0 = np.random.randn(d)
            x0 = x0 / np.linalg.norm(x0) * 5.0

            # --- Abelian errors (diagonal) ---
            x = x0.copy()
            dists_abelian[trial, 0] = np.linalg.norm(x)
            for step in range(n_steps):
                x = (1 - correction_rate) * x  # linear contraction
                error = np.diag(1 + eps * np.random.randn(d))
                x = error @ x
                dists_abelian[trial, step+1] = np.linalg.norm(x)

            # --- Non-abelian errors (rotation) ---
            x = x0.copy()
            dists_nonabelian[trial, 0] = np.linalg.norm(x)
            for step in range(n_steps):
                x = (1 - correction_rate) * x
                # Random skew-symmetric matrix → rotation
                A = eps * np.random.randn(d, d)
                A = (A - A.T) / 2  # skew-symmetric
                E = np.eye(d) + A + 0.5 * A @ A  # Approx exp(A) to 2nd order
                x = E @ x
                dists_nonabelian[trial, step+1] = np.linalg.norm(x)

        avg_a = np.mean(dists_abelian, axis=0)
        avg_na = np.mean(dists_nonabelian, axis=0)

        all_results[f"eps_{eps}"] = {
            "final_abelian": float(avg_a[-1]),
            "final_nonabelian": float(avg_na[-1]),
            "ratio": float(avg_a[-1] / max(avg_na[-1], 1e-10)),
        }

    # Plot for eps = 0.05
    eps = 0.05
    dists_abelian = np.zeros((n_trials, n_steps + 1))
    dists_nonabelian = np.zeros((n_trials, n_steps + 1))
    for trial in range(n_trials):
        x0 = np.random.randn(d)
        x0 = x0 / np.linalg.norm(x0) * 5.0
        x = x0.copy()
        dists_abelian[trial, 0] = np.linalg.norm(x)
        for step in range(n_steps):
            x = (1 - correction_rate) * x
            error = np.diag(1 + eps * np.random.randn(d))
            x = error @ x
            dists_abelian[trial, step+1] = np.linalg.norm(x)
        x = x0.copy()
        dists_nonabelian[trial, 0] = np.linalg.norm(x)
        for step in range(n_steps):
            x = (1 - correction_rate) * x
            A = eps * np.random.randn(d, d)
            A = (A - A.T) / 2
            E = np.eye(d) + A + 0.5 * A @ A
            x = E @ x
            dists_nonabelian[trial, step+1] = np.linalg.norm(x)

    avg_a = np.mean(dists_abelian, axis=0)
    avg_na = np.mean(dists_nonabelian, axis=0)
    std_a = np.std(dists_abelian, axis=0) / np.sqrt(n_trials)
    std_na = np.std(dists_nonabelian, axis=0) / np.sqrt(n_trials)

    steps = range(n_steps + 1)
    plt.plot(steps, avg_a, 'b-', linewidth=2, label='Abelian (diagonal)')
    plt.fill_between(steps, avg_a - 2*std_a, avg_a + 2*std_a, alpha=0.15, color='blue')
    plt.plot(steps, avg_na, 'r-', linewidth=2, label='Non-abelian (rotation)')
    plt.fill_between(steps, avg_na - 2*std_na, avg_na + 2*std_na, alpha=0.15, color='red')
    plt.xlabel("Proof Step")
    plt.ylabel("Distance from Correct Proof")
    plt.title(f"Proof Refinement Convergence (ε = {eps}, d = {d})")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.yscale('log')
    plt.tight_layout()
    plt.savefig("results/plots/proof_refinement_v2.png", dpi=150)
    plt.close()
    print("Saved: results/plots/proof_refinement_v2.png")

    # Report
    print(f"\nResults across error magnitudes (d={d}, correction_rate={correction_rate}):")
    print(f"{'ε':>6} {'Final d_A':>12} {'Final d_NA':>12} {'d_A/d_NA':>10}")
    print("-" * 45)
    for eps in error_magnitudes:
        r = all_results[f"eps_{eps}"]
        print(f"{eps:>6.2f} {r['final_abelian']:>12.6f} {r['final_nonabelian']:>12.6f} "
              f"{r['ratio']:>10.4f}")

    return all_results

def experiment_7_norm_growth():
    """
    Analyze the expected norm of products of error matrices.

    For abelian (diagonal) errors: E(‖E_n...E_1‖) = ∏ E(‖E_i‖)
    For non-abelian (rotation) errors: E(‖E_n...E_1‖) depends on commutators.

    Key result: The expected spectral norm of products of i.i.d. random matrices
    grows at rate determined by the Lyapunov exponent. For commuting (abelian)
    matrices, the Lyapunov exponent equals the expected log of the spectral radius.
    For non-commuting matrices, it can be larger.
    """
    print("\n" + "=" * 75)
    print("EXPERIMENT 7: Lyapunov Exponents of Error Products")
    print("=" * 75)

    d = 4
    n_steps = 500
    n_trials = 1000
    eps = 0.05

    # Abelian: products of diagonal matrices
    lyap_abelian = np.zeros(n_trials)
    for trial in range(n_trials):
        product = np.eye(d)
        for step in range(n_steps):
            E = np.diag(1 + eps * np.random.randn(d))
            product = E @ product
        lyap_abelian[trial] = np.log(np.linalg.norm(product, ord=2)) / n_steps

    # Non-abelian: products of near-identity matrices with off-diagonal terms
    lyap_nonabelian = np.zeros(n_trials)
    for trial in range(n_trials):
        product = np.eye(d)
        for step in range(n_steps):
            A = eps * np.random.randn(d, d)
            E = np.eye(d) + A
            product = E @ product
        lyap_nonabelian[trial] = np.log(max(np.linalg.norm(product, ord=2), 1e-300)) / n_steps

    # Non-abelian (skew-symmetric → near-rotation)
    lyap_rotation = np.zeros(n_trials)
    for trial in range(n_trials):
        product = np.eye(d)
        for step in range(n_steps):
            A = eps * np.random.randn(d, d)
            A = (A - A.T) / 2
            E = np.eye(d) + A + 0.5 * A @ A
            product = E @ product
        lyap_rotation[trial] = np.log(max(np.linalg.norm(product, ord=2), 1e-300)) / n_steps

    print(f"\nLyapunov exponents (d={d}, ε={eps}, {n_steps} steps, {n_trials} trials):")
    print(f"  Abelian (diagonal):     λ = {np.mean(lyap_abelian):.6f} ± {np.std(lyap_abelian):.6f}")
    print(f"  Non-abelian (general):  λ = {np.mean(lyap_nonabelian):.6f} ± {np.std(lyap_nonabelian):.6f}")
    print(f"  Non-abelian (rotation): λ = {np.mean(lyap_rotation):.6f} ± {np.std(lyap_rotation):.6f}")
    print(f"\n  Key: λ > 0 means errors GROW exponentially")
    print(f"       λ < 0 means errors SHRINK exponentially")
    print(f"       λ = 0 means errors stay bounded (isometry)")

    return {
        "lyapunov_abelian": float(np.mean(lyap_abelian)),
        "lyapunov_nonabelian_general": float(np.mean(lyap_nonabelian)),
        "lyapunov_nonabelian_rotation": float(np.mean(lyap_rotation)),
    }


if __name__ == "__main__":
    os.makedirs("results/plots", exist_ok=True)

    r1 = experiment_1_lazy_walk_comparison()
    experiment_2_convergence_curves()
    r3 = experiment_3_systematic_comparison()
    r4 = experiment_4_fixed_generators_fourier()
    experiment_5_diaconis_bound_analysis()
    r6 = experiment_6_proof_model()
    r7 = experiment_7_norm_growth()

    all_results = {
        "experiment_1_lazy_spectral": r1,
        "experiment_3_systematic": r3,
        "experiment_4_fixed_gen": r4,
        "experiment_6_proof_model": r6,
        "experiment_7_lyapunov": r7,
    }
    with open("results/metrics_v3.json", "w") as f:
        json.dump(all_results, f, indent=2, default=str)

    print("\n" + "=" * 75)
    print("All results saved.")
    print("=" * 75)
