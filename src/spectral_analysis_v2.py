"""
Spectral analysis v2: Using natural/restricted generators to reveal
how group structure affects error propagation convergence.

Key insight: The uniform walk on ALL elements trivially gives ρ = 1/(n-1).
The group structure manifests through the Cayley graph with specific generators.
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
# Group Constructions (same as v1)
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
        [0,1,2,3,4,5,6,7],
        [1,0,3,2,5,4,7,6],
        [2,3,1,0,6,7,5,4],
        [3,2,0,1,7,6,4,5],
        [4,5,7,6,1,0,2,3],
        [5,4,6,7,0,1,3,2],
        [6,7,4,5,3,2,1,0],
        [7,6,5,4,2,3,0,1],
    ])

def direct_product_cayley(table1, table2):
    """Cayley table for G1 × G2."""
    n1, n2 = len(table1), len(table2)
    order = n1 * n2
    table = np.zeros((order, order), dtype=int)
    for i in range(order):
        for j in range(order):
            i1, i2 = i // n2, i % n2
            j1, j2 = j // n2, j % n2
            table[i, j] = table1[i1, j1] * n2 + table2[i2, j2]
    return table

# ============================================================================
# Operator Construction
# ============================================================================

def random_walk_operator(cayley_table, generators, weights=None):
    n = len(cayley_table)
    if weights is None:
        weights = np.ones(len(generators)) / len(generators)
    P = np.zeros((n, n))
    for g_idx, w in zip(generators, weights):
        for i in range(n):
            j = cayley_table[i, g_idx]
            P[i, j] += w
    return P

def spectral_radius_nontrivial(P):
    eigenvalues = np.linalg.eigvals(P)
    abs_evals = sorted(np.abs(eigenvalues), reverse=True)
    return abs_evals[1] if len(abs_evals) > 1 else 0.0

def all_eigenvalues(P):
    return np.linalg.eigvals(P)

def spectral_gap(P):
    return 1.0 - spectral_radius_nontrivial(P)

def is_abelian(cayley_table):
    n = len(cayley_table)
    for i in range(n):
        for j in range(i+1, n):
            if cayley_table[i, j] != cayley_table[j, i]:
                return False
    return True

def compute_inverses(cayley_table):
    n = len(cayley_table)
    inv = np.zeros(n, dtype=int)
    for i in range(n):
        for j in range(n):
            if cayley_table[i, j] == 0:
                inv[i] = j
                break
    return inv

def commutator_subgroup(cayley_table):
    n = len(cayley_table)
    inv = compute_inverses(cayley_table)
    comms = set()
    for g in range(n):
        for h in range(n):
            gh = cayley_table[g, h]
            ginv_hinv = cayley_table[inv[g], inv[h]]
            comms.add(cayley_table[gh, ginv_hinv])
    # Generate subgroup
    subgroup = set(comms)
    changed = True
    while changed:
        changed = False
        new = set()
        for a in subgroup:
            for b in subgroup:
                p = cayley_table[a, b]
                if p not in subgroup:
                    new.add(p)
                    changed = True
        subgroup |= new
    return subgroup

# ============================================================================
# Natural Generator Sets
# ============================================================================

def cyclic_generators(n):
    """Z_n: generators {1, n-1} (i.e., +1 and -1)."""
    return [1, n-1]

def dihedral_generators(n):
    """D_n: generators {r, r^{-1}, s} = {1, n-1, n}."""
    return [1, n-1, n]

def symmetric_generators_adjacent(n, elem_to_idx):
    """S_n: adjacent transpositions (i, i+1)."""
    gens = []
    identity = list(range(n))
    for i in range(n-1):
        perm = list(identity)
        perm[i], perm[i+1] = perm[i+1], perm[i]
        gens.append(elem_to_idx[tuple(perm)])
    return gens

def quaternion_generators():
    """Q_8: generators {i, i^{-1}, j, j^{-1}} = {2, 3, 4, 5}."""
    return [2, 3, 4, 5]

# ============================================================================
# Main Experiments
# ============================================================================

def experiment_1_natural_generators():
    """
    Compare spectral radii using natural (minimal) generator sets.
    This reveals the true impact of group structure on convergence.
    """
    print("=" * 75)
    print("EXPERIMENT 1: Spectral Radius with Natural Generators")
    print("=" * 75)
    print("\nUsing natural/minimal generators (not all elements).\n")

    results = []

    # Z_6 vs D_3 (order 6)
    z6 = cyclic_group_cayley(6)
    d3 = dihedral_group_cayley(3)
    Pz6 = random_walk_operator(z6, cyclic_generators(6))
    Pd3 = random_walk_operator(d3, dihedral_generators(3))
    results.append({"name": "Z_6", "order": 6, "abelian": True, "n_gens": 2,
                     "rho": spectral_radius_nontrivial(Pz6), "gap": spectral_gap(Pz6)})
    results.append({"name": "D_3", "order": 6, "abelian": False, "n_gens": 3,
                     "rho": spectral_radius_nontrivial(Pd3), "gap": spectral_gap(Pd3)})

    # Z_8 vs Q_8 vs D_4 (order 8)
    z8 = cyclic_group_cayley(8)
    q8 = quaternion_group_cayley()
    d4 = dihedral_group_cayley(4)
    Pz8 = random_walk_operator(z8, cyclic_generators(8))
    Pq8 = random_walk_operator(q8, quaternion_generators())
    Pd4 = random_walk_operator(d4, dihedral_generators(4))
    results.append({"name": "Z_8", "order": 8, "abelian": True, "n_gens": 2,
                     "rho": spectral_radius_nontrivial(Pz8), "gap": spectral_gap(Pz8)})
    results.append({"name": "Q_8", "order": 8, "abelian": False, "n_gens": 4,
                     "rho": spectral_radius_nontrivial(Pq8), "gap": spectral_gap(Pq8)})
    results.append({"name": "D_4", "order": 8, "abelian": False, "n_gens": 3,
                     "rho": spectral_radius_nontrivial(Pd4), "gap": spectral_gap(Pd4)})

    # Z_10 vs D_5 (order 10)
    z10 = cyclic_group_cayley(10)
    d5 = dihedral_group_cayley(5)
    Pz10 = random_walk_operator(z10, cyclic_generators(10))
    Pd5 = random_walk_operator(d5, dihedral_generators(5))
    results.append({"name": "Z_10", "order": 10, "abelian": True, "n_gens": 2,
                     "rho": spectral_radius_nontrivial(Pz10), "gap": spectral_gap(Pz10)})
    results.append({"name": "D_5", "order": 10, "abelian": False, "n_gens": 3,
                     "rho": spectral_radius_nontrivial(Pd5), "gap": spectral_gap(Pd5)})

    # S_3 (=D_3, order 6) vs S_4 (order 24)
    s3_table, s3_elems, s3_idx = symmetric_group_cayley(3)
    s4_table, s4_elems, s4_idx = symmetric_group_cayley(4)
    s3_gens = symmetric_generators_adjacent(3, s3_idx)
    s4_gens = symmetric_generators_adjacent(4, s4_idx)
    Ps3 = random_walk_operator(s3_table, s3_gens)
    Ps4 = random_walk_operator(s4_table, s4_gens)
    results.append({"name": "S_3", "order": 6, "abelian": False, "n_gens": 2,
                     "rho": spectral_radius_nontrivial(Ps3), "gap": spectral_gap(Ps3)})
    results.append({"name": "S_4", "order": 24, "abelian": False, "n_gens": 3,
                     "rho": spectral_radius_nontrivial(Ps4), "gap": spectral_gap(Ps4)})

    # Z_2 x Z_3 (abelian, order 6) vs S_3 (non-abelian, order 6)
    z2 = cyclic_group_cayley(2)
    z3 = cyclic_group_cayley(3)
    z2xz3 = direct_product_cayley(z2, z3)
    # generators for Z2xZ3: (1,0) and (0,1) → indices 3 and 1
    z2xz3_gens = [3, 1]  # (1,0)=1*3+0=3, (0,1)=0*3+1=1
    Pz2xz3 = random_walk_operator(z2xz3, z2xz3_gens)
    results.append({"name": "Z2×Z3", "order": 6, "abelian": True, "n_gens": 2,
                     "rho": spectral_radius_nontrivial(Pz2xz3), "gap": spectral_gap(Pz2xz3)})

    print(f"{'Group':<12} {'Order':>6} {'Abelian':>8} {'#Gens':>6} {'ρ(P)':>10} {'Gap':>10}")
    print("-" * 60)
    for r in results:
        print(f"{r['name']:<12} {r['order']:>6} {'Yes' if r['abelian'] else 'No':>8} "
              f"{r['n_gens']:>6} {r['rho']:>10.6f} {r['gap']:>10.6f}")

    return results

def experiment_2_same_generators_different_structure():
    """
    CRITICAL EXPERIMENT: Compare groups with the SAME number of generators
    but different algebraic structure.

    Use 2-generator random walks uniformly:
    - Z_n with generators {1, n-1} (abelian)
    - D_{n/2} with generators {r, s} (non-abelian)
    """
    print("\n" + "=" * 75)
    print("EXPERIMENT 2: Same #Generators, Different Structure")
    print("=" * 75)
    print("\nAll walks use exactly 2 generators with equal weight.\n")

    results = []

    for n in [3, 4, 5, 6, 7, 8, 10, 12]:
        order = 2 * n

        # Abelian: Z_{2n} with gens {1, 2n-1}
        z = cyclic_group_cayley(order)
        Pz = random_walk_operator(z, [1, order-1])
        rho_z = spectral_radius_nontrivial(Pz)

        # Non-abelian: D_n with gens {r, s} = {1, n}
        d = dihedral_group_cayley(n)
        Pd = random_walk_operator(d, [1, n])
        rho_d = spectral_radius_nontrivial(Pd)

        comm = commutator_subgroup(d)
        comm_ratio = len(comm) / order

        results.append({
            "order": order,
            "abelian_name": f"Z_{order}",
            "nonabelian_name": f"D_{n}",
            "rho_abelian": rho_z,
            "rho_nonabelian": rho_d,
            "gap_abelian": 1 - rho_z,
            "gap_nonabelian": 1 - rho_d,
            "comm_ratio": comm_ratio,
            "abelian_faster": rho_z < rho_d,
        })

    print(f"{'Order':>6} {'Z_{2n}':>12} {'ρ(Z)':>10} {'D_n':>12} {'ρ(D)':>10} "
          f"{'|[G,G]|/|G|':>12} {'Z faster?':>10}")
    print("-" * 80)
    for r in results:
        faster = "YES ✓" if r["abelian_faster"] else "NO ✗"
        print(f"{r['order']:>6} {r['abelian_name']:>12} {r['rho_abelian']:>10.6f} "
              f"{r['nonabelian_name']:>12} {r['rho_nonabelian']:>10.6f} "
              f"{r['comm_ratio']:>12.4f} {faster:>10}")

    n_faster = sum(1 for r in results if r["abelian_faster"])
    print(f"\nAbelian faster in {n_faster}/{len(results)} cases")
    return results

def experiment_3_random_generator_sampling():
    """
    For groups of the same order, sample random k-element generator sets
    and compare spectral radii. This averages over generator choice.
    """
    print("\n" + "=" * 75)
    print("EXPERIMENT 3: Random Generator Sampling")
    print("=" * 75)

    n_samples = 200
    results = {}

    for n in [3, 4, 5, 6]:
        order = 2 * n
        z = cyclic_group_cayley(order)
        d = dihedral_group_cayley(n)

        rho_z_list = []
        rho_d_list = []

        for _ in range(n_samples):
            k = 2  # use 2 generators
            gens_z = list(np.random.choice(range(1, order), size=k, replace=False))
            gens_d = list(np.random.choice(range(1, order), size=k, replace=False))

            Pz = random_walk_operator(z, gens_z)
            Pd = random_walk_operator(d, gens_d)

            rho_z_list.append(spectral_radius_nontrivial(Pz))
            rho_d_list.append(spectral_radius_nontrivial(Pd))

        rho_z_arr = np.array(rho_z_list)
        rho_d_arr = np.array(rho_d_list)

        results[f"order_{order}"] = {
            "abelian_mean": float(np.mean(rho_z_arr)),
            "abelian_std": float(np.std(rho_z_arr)),
            "nonabelian_mean": float(np.mean(rho_d_arr)),
            "nonabelian_std": float(np.std(rho_d_arr)),
        }

        print(f"\nOrder {order}: Z_{order} vs D_{n} (k=2 generators, {n_samples} samples)")
        print(f"  Z_{order}: ρ = {np.mean(rho_z_arr):.4f} ± {np.std(rho_z_arr):.4f}")
        print(f"  D_{n}:  ρ = {np.mean(rho_d_arr):.4f} ± {np.std(rho_d_arr):.4f}")
        print(f"  Mean difference (ρ_D - ρ_Z): {np.mean(rho_d_arr) - np.mean(rho_z_arr):+.4f}")

    return results

def experiment_4_convergence_distance():
    """
    Directly simulate the convergence of the n-th convolution power to uniform.
    P^n → (1/|G|) J as n → ∞. Measure ‖P^n - uniform‖ at each step.
    """
    print("\n" + "=" * 75)
    print("EXPERIMENT 4: Power Iteration Convergence ‖P^n - uniform‖")
    print("=" * 75)

    n_steps = 60

    groups = [
        ("Z_6", cyclic_group_cayley(6), [1, 5]),
        ("D_3", dihedral_group_cayley(3), [1, 3]),
        ("Z_8", cyclic_group_cayley(8), [1, 7]),
        ("D_4", dihedral_group_cayley(4), [1, 4]),
        ("Z_10", cyclic_group_cayley(10), [1, 9]),
        ("D_5", dihedral_group_cayley(5), [1, 5]),
    ]

    convergence_data = {}

    plt.figure(figsize=(12, 5))

    for name, table, gens in groups:
        order = len(table)
        P = random_walk_operator(table, gens)
        uniform = np.ones((order, order)) / order

        distances = []
        Pn = np.eye(order)
        for step in range(n_steps):
            Pn = Pn @ P
            dist = np.max(np.abs(Pn - uniform))  # L∞ distance
            distances.append(dist)

        convergence_data[name] = distances
        abelian = is_abelian(table)
        style = '-' if abelian else '--'
        color = 'blue' if abelian else 'red'
        plt.semilogy(range(1, n_steps+1), distances, style, label=name,
                    linewidth=2, alpha=0.8)

    plt.xlabel("Step n")
    plt.ylabel("‖P^n - uniform‖_∞")
    plt.title("Convergence of Error Propagation: Abelian (solid) vs Non-Abelian (dashed)")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig("results/plots/convergence_comparison.png", dpi=150)
    plt.close()
    print("Saved: results/plots/convergence_comparison.png")

    # Report convergence rates
    print(f"\n{'Group':<8} {'Abelian':>8} {'ρ(P)':>10} {'‖P^10‖':>10} {'‖P^30‖':>10}")
    print("-" * 52)
    for name, table, gens in groups:
        P = random_walk_operator(table, gens)
        rho = spectral_radius_nontrivial(P)
        abelian = is_abelian(table)
        d10 = convergence_data[name][9]
        d30 = convergence_data[name][29]
        print(f"{name:<8} {'Yes' if abelian else 'No':>8} {rho:>10.6f} {d10:>10.6e} {d30:>10.6e}")

    return convergence_data

def experiment_5_spectral_decomposition():
    """
    Full spectral decomposition to understand WHY abelian groups converge differently.
    """
    print("\n" + "=" * 75)
    print("EXPERIMENT 5: Full Spectral Decomposition Analysis")
    print("=" * 75)

    groups = [
        ("Z_6", cyclic_group_cayley(6), [1, 5]),
        ("D_3", dihedral_group_cayley(3), [1, 3]),
        ("Z_8", cyclic_group_cayley(8), [1, 7]),
        ("D_4", dihedral_group_cayley(4), [1, 4]),
        ("Z_12", cyclic_group_cayley(12), [1, 11]),
        ("D_6", dihedral_group_cayley(6), [1, 6]),
    ]

    results = []

    for name, table, gens in groups:
        order = len(table)
        P = random_walk_operator(table, gens)
        evals = np.linalg.eigvals(P)
        evals_abs = sorted(np.abs(evals), reverse=True)
        abelian = is_abelian(table)

        # Count distinct eigenvalue magnitudes
        unique_mags = [evals_abs[0]]
        for m in evals_abs[1:]:
            if abs(m - unique_mags[-1]) > 1e-8:
                unique_mags.append(m)

        # Multiplicities
        multiplicities = {}
        for m in evals_abs:
            key = round(m, 6)
            multiplicities[key] = multiplicities.get(key, 0) + 1

        rho = evals_abs[1] if len(evals_abs) > 1 else 0
        max_mult = max(multiplicities.values())

        print(f"\n{name} (order {order}, {'abelian' if abelian else 'non-abelian'}):")
        print(f"  ρ(P) = {rho:.6f}")
        print(f"  # distinct |eigenvalues|: {len(unique_mags)}")
        print(f"  Max multiplicity: {max_mult}")
        print(f"  Eigenvalue magnitudes: {[round(m, 4) for m in unique_mags]}")
        print(f"  Multiplicities: {multiplicities}")

        results.append({
            "name": name, "order": order, "abelian": abelian,
            "rho": float(rho), "n_distinct": len(unique_mags),
            "max_multiplicity": max_mult,
            "eigenvalues": [complex(e) for e in sorted(evals, key=lambda x: -abs(x))],
        })

    return results

def experiment_6_proof_step_simulation():
    """
    Simulate a concrete proof refinement process.

    Model: Proof state is a vector in R^d. Correct proof = origin.
    Each step applies a "reasoning operator" that should move toward origin,
    but with error drawn from a group-structured error model.

    Abelian model: errors are diagonal matrices (commuting)
    Non-abelian model: errors are rotation matrices (non-commuting)
    """
    print("\n" + "=" * 75)
    print("EXPERIMENT 6: Simulated Proof Refinement")
    print("=" * 75)

    d = 4  # proof state dimension
    n_steps = 100
    n_trials = 1000
    epsilon = 0.05  # error magnitude

    results = {"abelian": [], "nonabelian": []}

    for trial in range(n_trials):
        # Initial proof state (away from origin)
        x0 = np.random.randn(d)
        x0 = x0 / np.linalg.norm(x0) * 5.0  # start at distance 5

        # ---- Abelian error model ----
        x = x0.copy()
        distances_a = [np.linalg.norm(x)]
        for step in range(n_steps):
            # Correct step: move toward origin
            direction = -x / max(np.linalg.norm(x), 1e-10)
            step_size = 0.1
            x_correct = x + step_size * direction

            # Abelian error: diagonal perturbation (commuting)
            error = np.diag(1 + epsilon * np.random.randn(d))
            x = error @ x_correct
            distances_a.append(np.linalg.norm(x))

        results["abelian"].append(distances_a)

        # ---- Non-abelian error model ----
        x = x0.copy()
        distances_na = [np.linalg.norm(x)]
        for step in range(n_steps):
            direction = -x / max(np.linalg.norm(x), 1e-10)
            step_size = 0.1
            x_correct = x + step_size * direction

            # Non-abelian error: random rotation (non-commuting)
            # Generate random rotation via QR decomposition
            A = np.eye(d) + epsilon * np.random.randn(d, d)
            Q, R = np.linalg.qr(A)
            # Ensure det = +1
            Q = Q * np.sign(np.diag(R))
            x = Q @ x_correct
            distances_na.append(np.linalg.norm(x))

        results["nonabelian"].append(distances_na)

    # Average over trials
    avg_a = np.mean(results["abelian"], axis=0)
    avg_na = np.mean(results["nonabelian"], axis=0)
    std_a = np.std(results["abelian"], axis=0)
    std_na = np.std(results["nonabelian"], axis=0)

    plt.figure(figsize=(10, 6))
    steps = range(n_steps + 1)
    plt.plot(steps, avg_a, 'b-', linewidth=2, label='Abelian errors (diagonal)')
    plt.fill_between(steps, avg_a - std_a, avg_a + std_a, alpha=0.2, color='blue')
    plt.plot(steps, avg_na, 'r--', linewidth=2, label='Non-abelian errors (rotation)')
    plt.fill_between(steps, avg_na - std_na, avg_na + std_na, alpha=0.2, color='red')
    plt.xlabel("Proof Step")
    plt.ylabel("Distance from Correct Proof")
    plt.title("Proof Refinement: Abelian vs Non-Abelian Error Models")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig("results/plots/proof_refinement_simulation.png", dpi=150)
    plt.close()
    print("Saved: results/plots/proof_refinement_simulation.png")

    # Compute convergence rates
    # Fit exponential decay: d(t) ≈ d(0) * exp(-r*t)
    from scipy.optimize import curve_fit

    def exp_decay(t, d0, r):
        return d0 * np.exp(-r * t)

    t = np.arange(len(avg_a))
    try:
        popt_a, _ = curve_fit(exp_decay, t, avg_a, p0=[5, 0.1], maxfev=5000)
        popt_na, _ = curve_fit(exp_decay, t, avg_na, p0=[5, 0.1], maxfev=5000)
        rate_a = popt_a[1]
        rate_na = popt_na[1]
    except:
        # Fallback: linear fit on log
        valid_a = avg_a[avg_a > 0.01]
        valid_na = avg_na[avg_na > 0.01]
        if len(valid_a) > 5:
            slope_a, _ = np.polyfit(np.arange(len(valid_a)), np.log(valid_a), 1)
            rate_a = -slope_a
        else:
            rate_a = 0
        if len(valid_na) > 5:
            slope_na, _ = np.polyfit(np.arange(len(valid_na)), np.log(valid_na), 1)
            rate_na = -slope_na
        else:
            rate_na = 0

    print(f"\nConvergence rates (fitted):")
    print(f"  Abelian errors:     r = {rate_a:.6f}")
    print(f"  Non-abelian errors: r = {rate_na:.6f}")
    print(f"  Ratio (r_A / r_NA): {rate_a / max(rate_na, 1e-10):.4f}")
    print(f"  Abelian converges {'FASTER' if rate_a > rate_na else 'SLOWER'}")

    return {
        "rate_abelian": float(rate_a),
        "rate_nonabelian": float(rate_na),
        "final_distance_abelian": float(avg_a[-1]),
        "final_distance_nonabelian": float(avg_na[-1]),
    }


if __name__ == "__main__":
    os.makedirs("results/plots", exist_ok=True)

    r1 = experiment_1_natural_generators()
    r2 = experiment_2_same_generators_different_structure()
    r3 = experiment_3_random_generator_sampling()
    r4 = experiment_4_convergence_distance()
    r5 = experiment_5_spectral_decomposition()
    r6 = experiment_6_proof_step_simulation()

    all_results = {
        "experiment_1": r1,
        "experiment_2": [{k: v for k, v in r.items()} for r in r2],
        "experiment_3": r3,
        "experiment_6_proof_simulation": r6,
    }

    with open("results/metrics_v2.json", "w") as f:
        json.dump(all_results, f, indent=2, default=str)

    print("\n" + "=" * 75)
    print("All results saved to results/metrics_v2.json")
    print("Plots saved to results/plots/")
    print("=" * 75)
