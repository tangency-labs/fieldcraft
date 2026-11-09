# mitc6a-bryan-buckle-n8

128 genuine tri6 elements, 289 nodes. Fixed original tolerances. Coarse 32-triangle mesh failed: harmonic components up to 17.4% and Bryan buckling -3.54%; these are not accepted coarse-mesh validations. Same physical loading integrated on the refined shape functions; no tuning of stiffness, mass or tolerance.

Accepted 128-triangle load factor error -1.08%, with the fixed 2% band.

E=1000, nu=0.3, h=0.01, unit square, unit total reference compression
on x=1 with quadratic edge force weights; ux fixed on x=0, uy at the
origin, w on all edges, rx on x edges, ry on y edges, rz everywhere.
The unit preload is a normalization of the eigenproblem, not an asserted
finite 10% strain equilibrium; the critical multiplier scales it to a
small-strain state. The reference is Ncr=4*pi^2*D, D=E*h^3/[12*(1-nu^2)].
