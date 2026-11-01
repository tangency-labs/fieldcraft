# Thick-walled cylinder: stress linearisation through the wall (Lamé, ASME VIII-2 Annex 5-A)

The Lamé cylinder of `lame-cylinder` (a quarter of a cylinder a = 10, b = 20, one unit long,
plane strain, internal pressure p = 100, E = 210000, ν = 0.3; 8 × 16 hex20 elements), with a
stress classification line from the bore A = (10, 0, 0.5) to the outside B = (20, 0, 0.5), so
t = b − a = 10 and the line's local y is the hoop direction, z the axis.

With A₀ = p a²/(b² − a²) and B₀ = A₀ b², Lamé gives σ_θ = A₀ + B₀/r², σ_r = A₀ − B₀/r² and
σ_z = 2νA₀. Integrating them over the wall as Annex 5-A defines membrane, (1/t)∫σ dx, and
bending, (6/t²)∫σ (t/2 − x) dx, in closed form:

| quantity | closed form | value |
|---|---|---|
| hoop membrane | p a/(b − a) | 100.00 |
| hoop bending at A | (6B₀/t²)[c (1/a − 1/b) − ln(b/a)], c = (a + b)/2 | 45.482 |
| hoop peak at A, at B | σ_θ − membrane ∓ bending | 21.184, 12.149 |
| radial membrane | −p a/(a + b) | −33.333 |
| axial membrane | 2νA₀ | 20.000 |
| von Mises of the membrane tensor | | 116.24 |

The command samples each element's own extrapolated nodal stresses with its shape functions at 41 equally spaced
stations and integrates them with Simpson's rule.
