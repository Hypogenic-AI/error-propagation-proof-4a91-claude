"""
Verify the Lyapunov exponent ordering theorem:
λ_abelian ≤ λ_non-abelian (with equality iff errors commute).

Also verify the commutator norm prediction.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

np.random.seed(42)

def compute_lyapunov(d, eps, n_steps, n_trials, mode):
    """Compute top Lyapunov exponent for error products."""
    exponents = np.zeros(n_trials)
    for trial in range(n_trials):
        product = np.eye(d)
        for step in range(n_steps):
            if mode == "diagonal":
                A = np.diag(np.random.randn(d))
            elif mode == "general":
                A = np.random.randn(d, d) / np.sqrt(d)
            elif mode == "skew":
                A = np.random.randn(d, d)
                A = (A - A.T) / 2
            elif mode == "block_diagonal":
                # 2x2 block diagonal (abelian in blocks)
                A = np.zeros((d, d))
                for i in range(0, d, 2):
                    block = np.random.randn(2, 2)
                    if i+1 < d:
                        A[i:i+2, i:i+2] = block
            else:
                raise ValueError(f"Unknown mode: {mode}")

            E = np.eye(d) + eps * A
            product = E @ product

            # Periodic re-orthonormalization to prevent overflow
            if (step + 1) % 50 == 0:
                norm = np.linalg.norm(product, ord=2)
                if norm > 1e10 or norm < 1e-10:
                    exponents[trial] += np.log(norm)
                    product = product / norm

        exponents[trial] += np.log(max(np.linalg.norm(product, ord=2), 1e-300))
        exponents[trial] /= n_steps

    return np.mean(exponents), np.std(exponents) / np.sqrt(n_trials)


def experiment_1_eps_scaling():
    """Verify λ scales as ε² for abelian and has ε⁴ correction for non-abelian."""
    print("=" * 70)
    print("EXPERIMENT: Lyapunov Exponent vs ε")
    print("=" * 70)

    d = 4
    n_steps = 2000
    n_trials = 200
    epsilons = [0.01, 0.02, 0.03, 0.05, 0.07, 0.10, 0.15, 0.20]

    results = {"diagonal": [], "general": [], "skew": []}

    for eps in epsilons:
        for mode in ["diagonal", "general", "skew"]:
            mean, se = compute_lyapunov(d, eps, n_steps, n_trials, mode)
            results[mode].append((eps, mean, se))
            print(f"  ε={eps:.3f}, {mode:<10}: λ = {mean:+.8f} ± {se:.8f}")

    # Plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    for mode, label, color in [("diagonal", "Abelian (diagonal)", "blue"),
                                ("general", "Non-abelian (general)", "red"),
                                ("skew", "Isometric (skew-sym)", "green")]:
        eps_vals = [r[0] for r in results[mode]]
        lambda_vals = [r[1] for r in results[mode]]
        se_vals = [r[2] for r in results[mode]]

        ax1.errorbar(eps_vals, lambda_vals, yerr=[2*s for s in se_vals],
                    fmt='o-', label=label, color=color, capsize=3)

    ax1.axhline(y=0, color='k', linestyle='--', alpha=0.3)
    ax1.set_xlabel("ε (error magnitude)")
    ax1.set_ylabel("Lyapunov exponent λ")
    ax1.set_title("Lyapunov Exponent vs Error Magnitude")
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Plot λ/ε² to see the scaling
    for mode, label, color in [("diagonal", "Abelian", "blue"),
                                ("general", "Non-abelian", "red")]:
        eps_vals = [r[0] for r in results[mode]]
        lambda_vals = [r[1] for r in results[mode]]
        scaled = [l / (e**2) for l, e in zip(lambda_vals, eps_vals)]
        ax2.plot(eps_vals, scaled, 'o-', label=label, color=color)

    ax2.set_xlabel("ε")
    ax2.set_ylabel("λ / ε²")
    ax2.set_title("Scaled Lyapunov Exponent (should be constant at leading order)")
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig("results/plots/lyapunov_scaling.png", dpi=150)
    plt.close()
    print("\nSaved: results/plots/lyapunov_scaling.png")

    return results


def experiment_2_commutator_prediction():
    """
    Verify: λ_general - λ_diagonal = C * ε⁴ * E[‖[A,B]‖²_F]
    """
    print("\n" + "=" * 70)
    print("EXPERIMENT: Commutator Norm as Predictor")
    print("=" * 70)

    d = 4
    n_steps = 3000
    n_trials = 300
    eps = 0.05

    # Compute expected commutator norm for different distributions
    modes_and_comm = []

    # Mode 1: Diagonal (commutator = 0)
    comm_samples = []
    for _ in range(1000):
        A = np.diag(np.random.randn(d))
        B = np.diag(np.random.randn(d))
        comm_samples.append(np.linalg.norm(A @ B - B @ A, 'fro')**2)
    comm_diag = np.mean(comm_samples)

    # Mode 2: General matrices
    comm_samples = []
    for _ in range(1000):
        A = np.random.randn(d, d) / np.sqrt(d)
        B = np.random.randn(d, d) / np.sqrt(d)
        comm_samples.append(np.linalg.norm(A @ B - B @ A, 'fro')**2)
    comm_general = np.mean(comm_samples)

    # Mode 3: Skew-symmetric
    comm_samples = []
    for _ in range(1000):
        A = np.random.randn(d, d); A = (A - A.T)/2
        B = np.random.randn(d, d); B = (B - B.T)/2
        comm_samples.append(np.linalg.norm(A @ B - B @ A, 'fro')**2)
    comm_skew = np.mean(comm_samples)

    # Mode 4: Block diagonal (partially commuting)
    comm_samples = []
    for _ in range(1000):
        A = np.zeros((d,d)); B = np.zeros((d,d))
        for i in range(0, d, 2):
            A[i:i+2, i:i+2] = np.random.randn(2, 2)
            B[i:i+2, i:i+2] = np.random.randn(2, 2)
        comm_samples.append(np.linalg.norm(A @ B - B @ A, 'fro')**2)
    comm_block = np.mean(comm_samples)

    print(f"\nExpected commutator norms E[‖[A,B]‖²_F] (d={d}):")
    print(f"  Diagonal:      {comm_diag:.6f}")
    print(f"  General:       {comm_general:.6f}")
    print(f"  Skew-symmetric:{comm_skew:.6f}")
    print(f"  Block-diagonal:{comm_block:.6f}")

    # Compute Lyapunov exponents
    print(f"\nLyapunov exponents (ε={eps}):")
    for mode in ["diagonal", "general", "skew", "block_diagonal"]:
        mean, se = compute_lyapunov(d, eps, n_steps, n_trials, mode)
        print(f"  {mode:<16}: λ = {mean:+.8f} ± {se:.8f}")

    # Compute gap
    l_diag, _ = compute_lyapunov(d, eps, n_steps, n_trials, "diagonal")
    l_gen, _ = compute_lyapunov(d, eps, n_steps, n_trials, "general")
    l_block, _ = compute_lyapunov(d, eps, n_steps, n_trials, "block_diagonal")

    print(f"\nLyapunov gaps (relative to diagonal):")
    print(f"  General - Diagonal:       {l_gen - l_diag:+.8f}")
    print(f"  Block_diag - Diagonal:    {l_block - l_diag:+.8f}")
    print(f"  Predicted ratio (comm norms): {comm_block / max(comm_general, 1e-10):.4f}")
    print(f"  Actual ratio (lyap gaps):     {(l_block - l_diag) / max(l_gen - l_diag, 1e-10):.4f}")


def experiment_3_dimension_scaling():
    """How does the Lyapunov gap scale with dimension d?"""
    print("\n" + "=" * 70)
    print("EXPERIMENT: Dimension Scaling of Lyapunov Gap")
    print("=" * 70)

    eps = 0.05
    n_steps = 2000
    n_trials = 200
    dims = [2, 3, 4, 6, 8, 10]

    results = []
    for d in dims:
        l_diag, se_d = compute_lyapunov(d, eps, n_steps, n_trials, "diagonal")
        l_gen, se_g = compute_lyapunov(d, eps, n_steps, n_trials, "general")
        gap = l_gen - l_diag
        results.append({"d": d, "l_diag": l_diag, "l_gen": l_gen, "gap": gap})
        print(f"  d={d:>2}: λ_diag={l_diag:+.6f}, λ_gen={l_gen:+.6f}, gap={gap:+.6f}")

    # Plot
    plt.figure(figsize=(8, 5))
    ds = [r["d"] for r in results]
    gaps = [r["gap"] for r in results]
    plt.plot(ds, gaps, 'ro-', linewidth=2, markersize=8)
    plt.xlabel("Dimension d")
    plt.ylabel("λ_general - λ_diagonal")
    plt.title("Non-commutativity Penalty vs Dimension")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig("results/plots/dimension_scaling.png", dpi=150)
    plt.close()
    print("Saved: results/plots/dimension_scaling.png")

    return results


if __name__ == "__main__":
    os.makedirs("results/plots", exist_ok=True)
    r1 = experiment_1_eps_scaling()
    experiment_2_commutator_prediction()
    experiment_3_dimension_scaling()
    print("\n" + "=" * 70)
    print("All Lyapunov verification experiments complete.")
    print("=" * 70)
