"""
Spectral analysis of error propagation operators on finite groups.

Computes spectral radii, mixing times, and convergence rates for
abelian vs non-abelian error propagation structures.
"""

import numpy as np
from itertools import product as iterproduct
import json
import os

np.random.seed(42)

# ============================================================================
# Group Constructions
# ============================================================================

def cyclic_group_cayley(n):
    """Cayley table for Z_n (abelian)."""
    return np.array([[(i + j) % n for j in range(n)] for i in range(n)])

def dihedral_group_cayley(n):
    """
    Cayley table for D_n (dihedral group of order 2n, non-abelian for n>=3).
    Elements: r^0, r^1, ..., r^{n-1}, s, sr, sr^2, ..., sr^{n-1}
    where r = rotation, s = reflection.
    Multiplication: r^a * r^b = r^{(a+b)%n}, r^a * sr^b = sr^{(b-a)%n},
                    sr^a * r^b = sr^{(a+b)%n}, sr^a * sr^b = r^{(a-b)%n}
    """
    order = 2 * n
    table = np.zeros((order, order), dtype=int)
    for i in range(order):
        for j in range(order):
            if i < n and j < n:
                # r^i * r^j = r^{(i+j)%n}
                table[i, j] = (i + j) % n
            elif i < n and j >= n:
                # r^i * sr^{j-n} = sr^{(j-n-i)%n}
                table[i, j] = n + (j - n - i) % n
            elif i >= n and j < n:
                # sr^{i-n} * r^j = sr^{(i-n+j)%n}
                table[i, j] = n + (i - n + j) % n
            else:
                # sr^{i-n} * sr^{j-n} = r^{(i-n-(j-n))%n}
                table[i, j] = (i - n - (j - n)) % n
    return table

def symmetric_group_elements(n):
    """Generate all permutations of {0,...,n-1} as tuples."""
    from itertools import permutations
    return list(permutations(range(n)))

def symmetric_group_cayley(n):
    """Cayley table for S_n (non-abelian for n>=3)."""
    elems = symmetric_group_elements(n)
    elem_to_idx = {e: i for i, e in enumerate(elems)}
    order = len(elems)
    table = np.zeros((order, order), dtype=int)
    for i, g in enumerate(elems):
        for j, h in enumerate(elems):
            # composition: (g*h)(x) = g(h(x))
            gh = tuple(g[h[k]] for k in range(n))
            table[i, j] = elem_to_idx[gh]

    return table

def quaternion_group_cayley():
    """Cayley table for Q_8 = {1,-1,i,-i,j,-j,k,-k} (non-abelian, order 8)."""
    # Elements: 0=1, 1=-1, 2=i, 3=-i, 4=j, 5=-j, 6=k, 7=-k
    # Multiplication rules: i^2=j^2=k^2=-1, ij=k, ji=-k, etc.
    table = np.array([
        [0,1,2,3,4,5,6,7],  # 1*
        [1,0,3,2,5,4,7,6],  # -1*
        [2,3,1,0,6,7,5,4],  # i*
        [3,2,0,1,7,6,4,5],  # -i*
        [4,5,7,6,1,0,2,3],  # j*
        [5,4,6,7,0,1,3,2],  # -j*
        [6,7,4,5,3,2,1,0],  # k*
        [7,6,5,4,2,3,0,1],  # -k*
    ])
    return table

# ============================================================================
# Error Propagation Operator Construction
# ============================================================================

def random_walk_operator(cayley_table, generators=None, weights=None):
    """
    Construct the transition matrix for a random walk on a group.

    If generators is None, use all non-identity elements uniformly.
    The operator P has P[i,j] = prob of transitioning from state i to state j.
    """
    n = len(cayley_table)
    if generators is None:
        generators = list(range(1, n))  # all non-identity
    if weights is None:
        weights = np.ones(len(generators)) / len(generators)

    P = np.zeros((n, n))
    for g_idx, w in zip(generators, weights):
        for i in range(n):
            j = cayley_table[i, g_idx]
            P[i, j] += w
    return P

def spectral_radius_excluding_trivial(P):
    """
    Compute spectral radius of P restricted to non-trivial part.

    For a random walk on a group, the trivial eigenvalue is 1.
    The spectral radius we care about is the second largest |eigenvalue|.
    """
    eigenvalues = np.linalg.eigvals(P)
    # Sort by absolute value
    abs_evals = np.abs(eigenvalues)
    sorted_evals = np.sort(abs_evals)[::-1]

    # The largest should be ~1 (trivial). Return second largest.
    if len(sorted_evals) > 1:
        return sorted_evals[1]
    return 0.0

def spectral_gap(P):
    """Spectral gap = 1 - spectral_radius_excluding_trivial."""
    return 1.0 - spectral_radius_excluding_trivial(P)

def mixing_time_bound(P, epsilon=0.01):
    """Upper bound on mixing time: t_mix <= (1/gamma) * log(n/epsilon)."""
    n = len(P)
    gamma = spectral_gap(P)
    if gamma <= 1e-12:
        return float('inf')
    return (1.0 / gamma) * np.log(n / epsilon)

# ============================================================================
# Convergence Simulation
# ============================================================================

def simulate_convergence(cayley_table, generators, n_steps=200, n_trials=1000):
    """
    Simulate random walk convergence to uniform distribution.

    Returns total variation distance at each step.
    """
    n = len(cayley_table)
    uniform = np.ones(n) / n

    # Track distribution evolution
    tv_distances = np.zeros(n_steps)

    for trial in range(n_trials):
        current = 0  # start at identity
        counts = np.zeros(n)

        for step in range(n_steps):
            g = np.random.choice(generators)
            current = cayley_table[current, g]
            counts[current] += 1

            if (step + 1) % max(1, n_steps // 100) == 0 or step < 20:
                # Compute empirical distribution
                empirical = counts / (step + 1)
                tv = 0.5 * np.sum(np.abs(empirical - uniform))
                tv_distances[step] += tv

    tv_distances /= n_trials
    return tv_distances

def compute_convergence_rate_empirical(tv_distances, start=10, end=None):
    """Fit exponential decay to estimate convergence rate."""
    if end is None:
        end = len(tv_distances)

    # Find non-zero entries
    valid = []
    for i in range(start, end):
        if tv_distances[i] > 1e-10:
            valid.append((i, tv_distances[i]))

    if len(valid) < 5:
        return 0.0

    # Linear regression on log(TV)
    x = np.array([v[0] for v in valid])
    y = np.log(np.array([v[1] for v in valid]))

    # Filter out non-finite
    mask = np.isfinite(y)
    if np.sum(mask) < 3:
        return 0.0

    x, y = x[mask], y[mask]
    slope, intercept = np.polyfit(x, y, 1)
    return -slope  # positive = faster convergence

# ============================================================================
# Commutator Analysis
# ============================================================================

def compute_commutator_subgroup(cayley_table):
    """Compute [G,G] = {ghg^{-1}h^{-1} : g,h in G}."""
    n = len(cayley_table)

    # First compute inverses
    inverses = np.zeros(n, dtype=int)
    for i in range(n):
        for j in range(n):
            if cayley_table[i, j] == 0:  # identity is element 0
                inverses[i] = j
                break

    commutators = set()
    for g in range(n):
        for h in range(n):
            # [g,h] = g*h*g^{-1}*h^{-1}
            gh = cayley_table[g, h]
            ginv_hinv = cayley_table[inverses[g], inverses[h]]
            comm = cayley_table[gh, ginv_hinv]
            commutators.add(comm)

    # Generate the subgroup from commutators
    subgroup = set(commutators)
    changed = True
    while changed:
        changed = False
        new = set()
        for a in subgroup:
            for b in subgroup:
                prod = cayley_table[a, b]
                if prod not in subgroup:
                    new.add(prod)
                    changed = True
        subgroup |= new

    return subgroup

def is_abelian(cayley_table):
    """Check if group is abelian."""
    n = len(cayley_table)
    for i in range(n):
        for j in range(i+1, n):
            if cayley_table[i, j] != cayley_table[j, i]:
                return False
    return True

# ============================================================================
# Main Experiments
# ============================================================================

def experiment_1_spectral_comparison():
    """
    Compare spectral radii across abelian and non-abelian groups of similar sizes.
    """
    print("=" * 70)
    print("EXPERIMENT 1: Spectral Radius Comparison (Abelian vs Non-Abelian)")
    print("=" * 70)

    results = []

    # Groups to test: (name, cayley_table, is_abelian)
    groups = [
        ("Z_6 (abelian, order 6)", cyclic_group_cayley(6), True),
        ("D_3 (non-abelian, order 6)", dihedral_group_cayley(3), False),
        ("Z_8 (abelian, order 8)", cyclic_group_cayley(8), True),
        ("Q_8 (non-abelian, order 8)", quaternion_group_cayley(), False),
        ("D_4 (non-abelian, order 8)", dihedral_group_cayley(4), False),
        ("Z_10 (abelian, order 10)", cyclic_group_cayley(10), True),
        ("D_5 (non-abelian, order 10)", dihedral_group_cayley(5), False),
        ("Z_24 (abelian, order 24)", cyclic_group_cayley(24), True),
        ("S_4 (non-abelian, order 24)", symmetric_group_cayley(4), False),
    ]

    print(f"\n{'Group':<30} {'Order':>6} {'Abelian':>8} {'ρ(P)':>10} {'Gap':>10} {'t_mix':>10} {'|[G,G]|':>8}")
    print("-" * 90)

    for name, table, expected_abelian in groups:
        order = len(table)
        verified_abelian = is_abelian(table)
        assert verified_abelian == expected_abelian, f"Abelian check failed for {name}"

        # Uniform random walk on all non-identity elements
        P = random_walk_operator(table)
        rho = spectral_radius_excluding_trivial(P)
        gap = 1.0 - rho
        tmix = mixing_time_bound(P)
        comm_size = len(compute_commutator_subgroup(table))

        result = {
            "name": name, "order": order, "abelian": verified_abelian,
            "spectral_radius": float(rho), "spectral_gap": float(gap),
            "mixing_time_bound": float(tmix), "commutator_size": comm_size
        }
        results.append(result)

        print(f"{name:<30} {order:>6} {'Yes' if verified_abelian else 'No':>8} "
              f"{rho:>10.6f} {gap:>10.6f} {tmix:>10.2f} {comm_size:>8}")

    # Analysis: Compare matched pairs
    print("\n--- Matched Pair Analysis ---")
    pairs = [
        ("Z_6 (abelian, order 6)", "D_3 (non-abelian, order 6)"),
        ("Z_8 (abelian, order 8)", "Q_8 (non-abelian, order 8)"),
        ("Z_8 (abelian, order 8)", "D_4 (non-abelian, order 8)"),
        ("Z_10 (abelian, order 10)", "D_5 (non-abelian, order 10)"),
        ("Z_24 (abelian, order 24)", "S_4 (non-abelian, order 24)"),
    ]

    result_map = {r["name"]: r for r in results}
    print(f"\n{'Abelian':<30} {'Non-Abelian':<30} {'ρ_A':>8} {'ρ_NA':>8} {'ρ_A < ρ_NA?':>12}")
    print("-" * 95)
    for a_name, na_name in pairs:
        ra = result_map[a_name]
        rna = result_map[na_name]
        comparison = "YES ✓" if ra["spectral_radius"] < rna["spectral_radius"] else "NO ✗"
        print(f"{a_name:<30} {na_name:<30} {ra['spectral_radius']:>8.4f} {rna['spectral_radius']:>8.4f} {comparison:>12}")

    return results

def experiment_2_convergence_simulation():
    """
    Simulate convergence of random walks on abelian vs non-abelian groups.
    """
    print("\n" + "=" * 70)
    print("EXPERIMENT 2: Convergence Rate Simulation")
    print("=" * 70)

    n_steps = 500
    results = {}

    test_groups = [
        ("Z_6", cyclic_group_cayley(6)),
        ("D_3", dihedral_group_cayley(3)),
        ("Z_8", cyclic_group_cayley(8)),
        ("Q_8", quaternion_group_cayley()),
        ("Z_24", cyclic_group_cayley(24)),
        ("S_4", symmetric_group_cayley(4)),
    ]

    for name, table in test_groups:
        order = len(table)
        generators = list(range(1, order))
        print(f"\nSimulating {name} (order {order})...")
        tv = simulate_convergence(table, generators, n_steps=n_steps, n_trials=500)
        rate = compute_convergence_rate_empirical(tv, start=20)
        results[name] = {"tv_distances": tv.tolist(), "empirical_rate": float(rate)}
        print(f"  Empirical convergence rate: {rate:.6f}")

    return results

def experiment_3_commutator_prediction():
    """
    Test whether commutator subgroup size predicts convergence penalty.
    """
    print("\n" + "=" * 70)
    print("EXPERIMENT 3: Commutator Subgroup as Convergence Predictor")
    print("=" * 70)

    groups = []
    for n in range(3, 8):
        # Cyclic group Z_{2n}
        table = cyclic_group_cayley(2*n)
        P = random_walk_operator(table)
        groups.append({
            "name": f"Z_{2*n}",
            "order": 2*n,
            "abelian": True,
            "spectral_radius": spectral_radius_excluding_trivial(P),
            "commutator_ratio": 1.0 / (2*n),  # |[G,G]|/|G| = 1/|G| for abelian
        })

        # Dihedral group D_n
        table = dihedral_group_cayley(n)
        P = random_walk_operator(table)
        comm = compute_commutator_subgroup(table)
        groups.append({
            "name": f"D_{n}",
            "order": 2*n,
            "abelian": False,
            "spectral_radius": spectral_radius_excluding_trivial(P),
            "commutator_ratio": len(comm) / (2*n),
        })

    print(f"\n{'Group':<12} {'Order':>6} {'Abelian':>8} {'ρ(P)':>10} {'|[G,G]|/|G|':>12}")
    print("-" * 55)
    for g in groups:
        print(f"{g['name']:<12} {g['order']:>6} {'Yes' if g['abelian'] else 'No':>8} "
              f"{g['spectral_radius']:>10.6f} {g['commutator_ratio']:>12.4f}")

    # Correlation between commutator ratio and spectral radius
    non_abelian = [g for g in groups if not g["abelian"]]
    if len(non_abelian) >= 3:
        comm_ratios = np.array([g["commutator_ratio"] for g in non_abelian])
        spec_radii = np.array([g["spectral_radius"] for g in non_abelian])
        corr = np.corrcoef(comm_ratios, spec_radii)[0, 1]
        print(f"\nCorrelation (commutator ratio vs spectral radius, non-abelian): {corr:.4f}")

    return groups

def experiment_4_representation_theory():
    """
    Verify the representation-theoretic explanation:
    - Abelian groups: all irreps are 1-dimensional
    - Non-abelian groups: have irreps of dimension > 1
    - The spectral radius is determined by Fourier coefficients at irreps
    """
    print("\n" + "=" * 70)
    print("EXPERIMENT 4: Representation-Theoretic Analysis")
    print("=" * 70)

    # For small groups, compute the full spectral decomposition
    groups = [
        ("Z_6", cyclic_group_cayley(6)),
        ("D_3", dihedral_group_cayley(3)),
        ("Z_8", cyclic_group_cayley(8)),
        ("Q_8", quaternion_group_cayley()),
    ]

    for name, table in groups:
        order = len(table)
        P = random_walk_operator(table)
        eigenvalues = np.linalg.eigvals(P)

        # Sort by real part then imaginary part for cleaner output
        evals_sorted = sorted(eigenvalues, key=lambda x: (-abs(x), -x.real))

        abelian = is_abelian(table)
        print(f"\n{name} (order {order}, {'abelian' if abelian else 'non-abelian'}):")
        print(f"  Eigenvalues of P (sorted by |λ|):")

        # Count distinct eigenvalue magnitudes
        mags = np.sort(np.abs(eigenvalues))[::-1]
        unique_mags = []
        for m in mags:
            if not unique_mags or abs(m - unique_mags[-1]) > 1e-8:
                unique_mags.append(m)

        for i, ev in enumerate(evals_sorted[:min(10, len(evals_sorted))]):
            print(f"    λ_{i} = {ev.real:+.6f} {'+' if ev.imag >= 0 else ''}{ev.imag:.6f}i  "
                  f"(|λ| = {abs(ev):.6f})")

        print(f"  Number of distinct eigenvalue magnitudes: {len(unique_mags)}")
        print(f"  Spectral radius (excl. trivial): {spectral_radius_excluding_trivial(P):.6f}")

        if abelian:
            # For abelian: eigenvalues should be exactly the characters evaluated at the generator distribution
            print(f"  → All {order} eigenvalues correspond to 1-dim irreps (characters)")
        else:
            # Count multiplicities
            mult_count = {}
            for m in mags:
                m_round = round(m, 6)
                mult_count[m_round] = mult_count.get(m_round, 0) + 1
            print(f"  → Eigenvalue multiplicities: {dict(sorted(mult_count.items(), reverse=True))}")
            print(f"  → Multiplicities > 1 indicate higher-dim irreps")


if __name__ == "__main__":
    # Create results directory
    os.makedirs("results", exist_ok=True)

    # Run all experiments
    results1 = experiment_1_spectral_comparison()
    results2 = experiment_2_convergence_simulation()
    results3 = experiment_3_commutator_prediction()
    experiment_4_representation_theory()

    # Save results
    all_results = {
        "experiment_1_spectral_comparison": results1,
        "experiment_2_convergence_simulation": {
            k: {"empirical_rate": v["empirical_rate"]}
            for k, v in results2.items()
        },
        "experiment_3_commutator_prediction": results3,
    }

    with open("results/metrics.json", "w") as f:
        json.dump(all_results, f, indent=2, default=str)

    print("\n" + "=" * 70)
    print("All results saved to results/metrics.json")
    print("=" * 70)
