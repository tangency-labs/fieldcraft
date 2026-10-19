# Natural convection in a square cavity with a centred conducting body, Ra = 10⁵

Square cavity of side L, the left wall hot and the right wall cold, top and bottom adiabatic,
Pr = 0.71, Ra = g β ΔT L³ / (ν α) = 10⁵. A square body of side L / 2 sits at the centre and
conducts heat with k_s / k_f = 0.2, 1 or 5 (no heat generation; the steady state does not depend
on its heat capacity). The quantity is the average Nusselt number of the hot wall, on the
cavity's side and the fluid's conductivity. With no flow the published values are 0.7069, 1 and
1.4138.

The lattice Boltzmann solver checks this case in `crates/lbm/tests/conjugate.rs`
(`conjugate_cavity_matches_house_et_al`, heavy; the conduction limit runs every time) until its
oracle recipe lands (a follow-up: the oracle runs fluid steps since the de Vahl Davis cases).
