# Finite beam on a Winkler foundation

An original 32 by 2 MITC4 composite strip has L=1 m, b=0.1 m,
h=0.01 m, E=1 GPa, nu=0 and k=1e4 N/m^3. Cylindrical bending
constraints suppress transverse rotation and membrane motion. Both ends
have w=0 and free bending rotation. Consistent nodal forces integrate
q=sin(pi*x/L) Pa exactly against the linear shape functions.

With I=b*h^3/12, W=q0*b/(E*I*(pi/L)^4+k*b) and
M(x)=E*I*(pi/L)^2*W*sin(pi*x/L). The finite-length separated solution
satisfies the Winkler/Hetenyi equation and both end moment conditions.
The oracle measures w at midspan and the section moment from top/bottom
ply stresses, averaged over four in-plane quadrature points. It uses the
original 2% discretization band; no reference was fitted to FE output.

Source: NPTEL course 105106049, chapter 11, equations 11.1-11.7:
https://archive.nptel.ac.in/content/storage2/courses/105106049/lecnotes/mainch11.html
The command recipe is the original independently constructed verification
model; no commercial example or result was used.
