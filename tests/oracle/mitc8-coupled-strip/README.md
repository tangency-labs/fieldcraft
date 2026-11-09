# Unsymmetric two-layer bending

Same 2x1 quad8 strip mesh, L=2, b=1 and end moment 1e-4 as mitc8-pure-bending. Bottom layer E=1200, h=.04; upper E=3000, h=.06; both nu=0, isotropic. Right-hand normal +z, zero offset. A=228, B=2.16, D=.1756 about the nodal plane. End axial force is zero, so kappa=M/(D-B^2/A), epsilon=-B*kappa/A. Probe at the free-edge midpoint: w=-kappa L^2/2, ux=epsilon L. This is an independent coupled patch, not a symmetric Navier formula used outside its scope.
