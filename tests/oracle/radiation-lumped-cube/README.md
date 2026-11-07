# A lumped body cooling by radiation

A steel cube 10 mm on a side (ρ = 7800 kg/m³, c = 460 J/(kg·K), k = 15 W/(m·K)) starts at
T₀ = 1000 K and cools by radiation alone from its six faces (ε = 0.7) to large surroundings, at
300 K in step `Room` and at 0 K in step `Space`, for 600 s.

**Closed form.** A body with a uniform temperature obeys ρ c V dT/dt = −ε σ A (T⁴ − T_sur⁴), which
integrates to

t = ρ c V / (4 ε σ A T_sur³) {ln|(T_sur + T)/(T_sur − T)| − ln|(T_sur + T₀)/(T_sur − T₀)| + 2 [atan(T/T_sur) − atan(T₀/T_sur)]},

and to t = ρ c V / (3 ε σ A) (1/T³ − 1/T₀³) when T_sur = 0. T(600 s) is the root of these, by
bisection to round-off: 440.526745 K to 300 K, 425.863199 K to 0 K. The cube loses most of its
heat in the run.

**Why one hex8 is the lumped body:** with a uniform temperature its eight nodes each carry an
eighth of the capacity and of the radiating area, so the temperature stays uniform and the
conductivity never enters.

**Run:** transient heat transfer, Crank–Nicolson (θ = ½) with the radiation weighted by θ like
the conduction, 320 increments of 1.875 s; each increment solves the radiation by Newton's method.

**Time convergence** (`crates/solver/tests/radiation.rs`, the same cube):

| Increments | 80 | 160 | 320 (this case) |
|---|---|---|---|
| Crank–Nicolson error at 600 s, to 300 K | 0.0533 K | 0.0133 K | 0.0033 K |
| Backward Euler error | | 1.28 K | 0.64 K |

Crank–Nicolson converges at second order, backward Euler at first. The tolerance, 2·10⁻⁵ relative
(about 0.009 K), covers the 320-increment error (7.5·10⁻⁶ of T to 300 K, 8.6·10⁻⁶ to 0 K) with a
margin.

Reference: F. P. Incropera, D. P. DeWitt, *Fundamentals of Heat and Mass Transfer*, §5.3 (the
general lumped capacitance analysis, with radiation to large surroundings).
