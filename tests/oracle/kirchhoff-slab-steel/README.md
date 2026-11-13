# Nonlinear conduction through a steel slab (#838)

A slab 10 × 10 × 200 of 1 × 1 × 8 hex8 elements of the library's EN10025-S355 steel, whose
conductivity follows the EN 1993-1-2 carbon-steel curve k = 54 − 3.33·10⁻² θ W/(m K) (θ in °C).
A steady heat transfer step holds z = 0 at 20 °C and z = 200 at 700 °C; the solver iterates the
temperature-dependent conductivity (Picard) to convergence.

With k(T) the Kirchhoff transform θ = ∫k dT is linear in z, so T(z) is the root of
a T + b T²/2 = θ(z) (Carslaw and Jaeger 1959, §1.7). Linear elements taking k at their mean
temperature reproduce it exactly at the nodes, and their flux is uniform and exact. The profile
sags up to 45 K below the straight line a constant conductivity would give: the temperature
dependence acts.
