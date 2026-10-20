# Heated square cylinder in a plane channel (pending)

Laminar forced convection from a heated block: a square cylinder held at a constant temperature,
centred between the two walls of a plane channel, in air (Pr = 0.7), in the steady flow regime so
that the run can use the steady stop rule. The quantities are the cylinder's mean Nusselt number,
its front, side and rear faces' Nusselt numbers and its drag coefficient.

Reference: A. K. Dhiman, R. P. Chhabra and V. Eswaran, *Flow and heat transfer across a confined
square cylinder in the steady flow regime: effect of Peclet number*, International Journal of
Heat and Mass Transfer 48 (2005) 4598–4614, doi:10.1016/j.ijheatmasstransfer.2005.04.033
(confirmed on Crossref). Its abstract gives the range: 1 ≤ Re ≤ 45, 0.7 ≤ Pr ≤ 4000 and
blockages 1/8, 1/6 and 1/4, with constant-temperature and constant-flux cylinders.

The case is pending because the paper is not open access: neither its values nor the details of
its setup (how Re is defined, the inflow profile, the channel's length upstream and downstream)
have been read from it, and nothing is cited here that was not. Once they are transcribed (table
and page named in `source`), the case gets its recipe: the cylinder as an obstacle body
(`design.json`) with its faces held hot, adiabatic channel walls, the paper's inflow and a
pressure outlet, and the grid refined until each quantity moves by less than its tolerance.

The unconfined cylinder (Sharma & Eswaran, Numer. Heat Transfer A 45 (2004) 247–269,
doi:10.1080/10407780490278562, also not open access) is the alternative at Re ≈ 100; there the
flow sheds vortices, so its Nusselt number needs a time average rather than the steady stop rule.
