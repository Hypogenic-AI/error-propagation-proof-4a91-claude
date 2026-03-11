# Resources Catalog

## Summary
This document catalogs all resources gathered for the research project: "Group-Theoretic Characterization of Error Propagation in Multi-Step Mathematical Proofs."

## Papers
Total papers downloaded: 10

| Title | Authors | Year | File | Key Results |
|-------|---------|------|------|-------------|
| IEKF as Stable Observer | Barrau, Bonnabel | 2014 | papers/barrau_bonnabel_2014_IEKF.pdf | Log-linear error property on Lie groups; group affine characterization |
| Stabilizer Codes | Gottesman | 1996 | papers/gottesman_1996_stabilizer.pdf | Stabilizer must be abelian; quantum Hamming bound |
| Non-stabilizer from Abelian Subgroups | Arvind et al. | 2002 | papers/arvind_2002_abelian_subgroups.pdf | Error codes from abelian subgroups of non-abelian error group |
| IMU Error Propagation Framework | Barrau, Bonnabel | 2020 | papers/barrau_bonnabel_2020_framework.pdf | Log-linearity on SE_2(3); exact preintegration |
| Linear Systems on Groups | Barrau, Bonnabel | 2019 | papers/barrau_2019_linear_systems_groups.pdf | Generalization of linear systems to Lie groups |
| Uncertainty on Unimodular Groups | Ye, Chirikjian | 2025 | papers/ye_chirikjian_2025_uncertainty_unimodular.pdf | Exact mean/covariance propagation on Lie groups |
| Transversal Gates Error Propagation | Bombin | 2018 | papers/bombin_2018_transversal.pdf | Error propagation through sequential gate operations |
| Random Walks on Groups | Diaconis | 2003 | papers/diaconis_2003_random_walks_groups.pdf | Spectral gap governs mixing time on groups |
| Group-theoretic Error Mitigation | Zhao, Miyake | 2023 | papers/zhao_miyake_2023_group_error_mitigation.pdf | Group symmetries for error mitigation |
| Error Dynamics Affine Groups | Li, Chen et al. | 2023 | papers/li_chen_2023_error_dynamics_affine.pdf | Stochastic error dynamics on affine group systems |

See papers/README.md for detailed descriptions.

## Prior Results Catalog
Key theorems and lemmas available for our proofs:

| Result | Source | Statement Summary | Used For |
|--------|--------|-------------------|----------|
| Log-Linear Error Property (Thm 2) | Barrau & Bonnabel, 2014 | Lie logarithm of invariant error satisfies linear ODE for group affine systems | Converting nonlinear group error to linear spectral analysis |
| Group Affine Characterization (Thm 1) | Barrau & Bonnabel, 2014 | f_u(ab) = f_u(a)b + af_u(b) - af_u(Id)b ⟺ trajectory-independent error propagation | Identifying which proof systems have autonomous error dynamics |
| Abelian Stabilizer (Sec II) | Gottesman, 1996 | Error-correcting stabilizer H must be abelian, H ≅ (Z₂)^a | Establishing abelian structure as necessary for error correction |
| Error Detection via Anticommutation | Gottesman, 1996 | Error E detectable iff ∃M ∈ H with {E,M}=0 | Connecting commutativity structure to error detectability |
| Spectral Radius Convergence | Standard (functional analysis) | ρ(A) < 1 ⟹ A^n → 0; rate governed by ρ(A) | Characterizing convergence rate of error propagation |
| Spectral Gap Mixing Bound | Diaconis et al. | t_mix ≤ C/γ · log(1/ε) for spectral gap γ | Bounding convergence time |
| Abelian Representation Theory | Standard (group theory) | Abelian groups have only 1-dim irreps | Explaining simpler spectral structure for abelian errors |

## Computational Tools

| Tool | Purpose | Location | Notes |
|------|---------|----------|-------|
| SymPy | Symbolic computation | pip: sympy | For algebraic verification of group properties, commutator computations |
| NumPy | Numerical computation | pip: numpy | For spectral radius computation, matrix operations |

### Relevant SymPy Capabilities for This Research
- `sympy.combinatorics.PermutationGroup`: Group operations, subgroup testing
- `sympy.matrices`: Matrix Lie algebra computations, eigenvalue/spectral radius
- `sympy.physics.quantum`: Pauli matrices, commutator/anticommutator

## Resource Gathering Notes

### Search Strategy
1. Primary search via paper-finder service with queries covering:
   - Error propagation + group theory + proofs
   - Spectral radius + operators + proof complexity
   - Non-abelian groups + convergence + composition
   - Random walks + groups + spectral gap
   - Group-theoretic error correction (quantum)
   - Semigroup operator composition
   - LLM mathematical reasoning errors
2. 8+ separate searches conducted with "diligent" mode
3. 85+ papers reviewed from Semantic Scholar results
4. Top 10 downloaded based on relevance scoring

### Selection Criteria
- Papers providing mathematical framework adaptable to proof state spaces
- Papers establishing connections between group structure (abelian vs non-abelian) and convergence/error properties
- Papers with precise theorem statements and proof techniques reusable in our setting
- Foundational papers in Lie group error propagation (Barrau, Bonnabel, Chirikjian)
- Foundational papers in group-theoretic error correction (Gottesman, Arvind)

### Key Observation
The research topic is genuinely novel — no existing work directly applies group-theoretic error propagation to mathematical proof theory. The closest analogies are:
1. **Lie group error propagation** (robotics): Non-abelian group structure governs error dynamics, log-linear property enables spectral analysis
2. **Quantum error correction**: Group-theoretic (specifically abelian subgroup) structure determines error correctability
3. **Random walks on groups**: Spectral gap (determined by group structure) governs convergence rate

### Challenges Encountered
- Many arXiv IDs for papers on "error propagation on groups" resolve to navigation/robotics papers — the mathematical content is relevant but the application domain is different
- Some arXiv IDs resolved to different papers than expected (PDF download verification important)
- No papers found directly connecting group-theoretic error propagation to proof complexity or formal verification

## Recommendations for Proof Construction

1. **Proof strategy**: Build the framework in layers:
   - (a) Define proof state space as a group/monoid
   - (b) Show error operators form a group under composition
   - (c) Apply Barrau-Bonnabel log-linear property to get linear error dynamics
   - (d) Analyze spectral radius of the resulting linear operator
   - (e) Prove abelian subgroups yield smaller spectral radius

2. **Key prerequisites**:
   - Barrau-Bonnabel Theorem 2 (log-linear property)
   - Gottesman stabilizer formalism (abelian requirement)
   - Standard spectral radius convergence theorem
   - Representation theory of finite/compact groups

3. **Computational tools**: Use SymPy for:
   - Verifying group axioms on proof state spaces
   - Computing commutators and spectral radii for small examples
   - Constructing explicit abelian vs non-abelian error propagation examples

4. **Potential difficulties**:
   - Proof state spaces may be infinite-dimensional (need functional analysis tools)
   - The group affine property may need to be weakened for proof applications
   - The abelian/non-abelian convergence gap may depend on specific group structure
