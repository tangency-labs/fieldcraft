# Resonant modes of a rectangular cavity

A 20 × 10 × 15 mm box of air (εᵣ = μᵣ = 1) meshed with 8 × 4 × 6 blocks of tetrahedra, its six walls
perfect electric conductors, solved for its six lowest modes with second-order Nédélec edge elements
(`solve.cavity`). The closed form (Pozar §6.4) gives f_mnl = (c₀/2)·√((m/a)² + (n/b)² + (l/d)²), one
mode for each (m, n, l) with at most one index zero and two (a TE and a TM mode) when none is.

Refinement: the solver's tests (`crates/solver/tests/cavity.rs`) solve a box of the same proportions on
4 × 2 × 3 cells (largest k² error 0.81 % at order 2) and on 8 × 4 × 6 (the rate h⁴, observed 3.69),
so this mesh is converged well inside the 0.2 % tolerance.
