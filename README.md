# Group-Theoretic Characterization of Error Propagation in Multi-Step Mathematical Proofs

## Overview

This project investigates how the algebraic structure of error composition affects convergence in iterative proof refinement systems. We prove that commuting (abelian) error operators yield lower Lyapunov exponents than non-commuting ones, with the gap quantified by the commutator norm.

## Key Results

- **Theorem 1**: Abelian error operators have strictly lower Lyapunov exponents than non-abelian ones: λ_abel < λ_general, with gap proportional to E[‖[A₁, A₂]‖²_F].
- **Theorem 3**: For random walks on groups, the spectral gap ratio γ(D_n)/γ(Z_{2n}) → 4 as n → ∞, showing non-abelian groups can mix faster.
- **Theorem 5 (Main)**: Proof refinement with abelian errors has a strictly larger critical error threshold ε*_abel > ε*_non-abel, providing greater robustness.

## File Structure

```
├── REPORT.md              # Full research report with proofs and results
├── definitions.md         # Formal definitions and notation
├── planning.md            # Research plan
├── literature_review.md   # Literature review (pre-gathered)
├── resources.md           # Resource catalog
├── src/
│   ├── spectral_analysis.py      # Initial spectral comparison (all generators)
│   ├── spectral_analysis_v2.py   # Natural generators, proof simulation
│   ├── spectral_analysis_v3.py   # Lazy walks, Diaconis bounds, Lyapunov
│   └── lyapunov_verification.py  # Lyapunov exponent scaling verification
├── results/
│   ├── proofs.md          # Detailed formal proofs
│   ├── metrics.json       # Raw results (v1)
│   ├── metrics_v2.json    # Raw results (v2)
│   ├── metrics_v3.json    # Raw results (v3)
│   └── plots/
│       ├── convergence_comparison.png
│       ├── lazy_walk_convergence.png
│       ├── spectral_gap_comparison.png
│       ├── proof_refinement_simulation.png
│       ├── proof_refinement_v2.png
│       ├── lyapunov_scaling.png
│       └── dimension_scaling.png
└── papers/                # Reference papers (PDFs)
```

## Reproducing Results

```bash
source .venv/bin/activate
python src/spectral_analysis_v3.py   # Main experiments
python src/lyapunov_verification.py  # Lyapunov scaling verification
```

See [REPORT.md](REPORT.md) for full details.
