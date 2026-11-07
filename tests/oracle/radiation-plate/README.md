# A plate cooled by radiation

A steel plate (k = 15 W/(m·K)) 0.02 m thick and 1 m × 1 m across is heated on one face by a
uniform flux q and loses that heat from the other face to large surroundings; its edges are
insulated. In steady state everything q brings in leaves the cooled face, so the face's energy
balance fixes its temperature, and conduction carries the heat across the plate on a straight
line, T(z) = T_L + q (L − z) / k.

**Three cases on one mesh** (hex8, four elements through the thickness):

| Step | Heat in | Cooled face | Face temperature T_L |
|---|---|---|---|
| `Room` | q = 2000 W/m² | radiation, ε = 0.8, to 293.15 K | (q/εσ + T_amb⁴)^¼ = 476.3179 K |
| `Space` | q = 2000 W/m² | radiation, ε = 0.8, to 0 K | (q/εσ)^¼ = 458.2285 K |
| `Mixed` | q = 1500 W/m² | radiation, ε = 0.9, to 300 K, and convection, h = 10 W/(m²·K), to 295 K | root of h (T − T∞) + εσ (T⁴ − T_amb⁴) = q: 379.9649 K |

σ = 5.670374419·10⁻⁸ W/(m²·K⁴). In `Mixed` radiation carries 650.35 W of the 1500 W, about as much
as convection in still air.

**Measures:** the temperatures of both faces (`solve.probe` at the nodes on z = 0 and z = L) and
the heat each load takes out (the solve's `heat_flows`). The field is linear, which every hex8
reproduces exactly, so the results match to the Newton tolerance; the tolerances are 10⁻⁹
relative for the temperatures and the radiated heat that is all of q, 10⁻⁸ for the split
between radiation and convection.

The model is in SI with temperatures in kelvin (`model.units`): radiation needs absolute
temperatures, so 293.15 K is a room at 20 °C.

Reference: F. P. Incropera, D. P. DeWitt, *Fundamentals of Heat and Mass Transfer*, §1.2.3
(radiation exchange with large surroundings) and §1.3.1 (the surface energy balance).
