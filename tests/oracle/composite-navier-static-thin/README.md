# Navier sinusoidal laminate plate: thin

A 1 × 1 m square of thickness 0.01 m has an 8 × 8 MITC4 mesh, [0/90/90/0] equal plies, E=[100,10,10] GPa, nu=0 and G=[5,3.5,3.5] GPa. Directions are fixed to global x/y. The hard simply supported constraints are the same as in the frequency cases.

q=sin(pi*x)sin(pi*y) Pa acts in +z. The nodal forces in the recipe are original consistent load integrals, with a four-point Gauss rule in each element direction. They are not sampled nodal pressure times area.

For alpha=beta=pi, use the 3 × 3 matrix printed in `../composite-navier-frequency-thin/README.md`; solve K·[X,Y,W]=[0,0,1]. Equivalently, W=1/[2*S*pi²-(S*pi)²*(K11+K22-2*K12)/(K11*K22-K12²)]. This gives W_FSDT=9.498338500230013e-07 m. The independent CLT limit is W_CLT=1/[(D11+4*D66+D22)*pi^4]=9.476291312016311e-07 m. B=D16=D26=0; the scalar separated solution must not be applied to the balanced angle-ply case with nonzero D16/D26.

The fixed 3% acceptance band matches the original FEM kernel regression. For h/L=.01 an additional CLT check tests the thin limit; the thick case uses FSDT and retains its shear effect. These are model-theory checks, not equilibrated interlaminar-stress or three-dimensional validation.

Run `cargo xtask oracle --case composite-navier-static-thin` and freeze only after the reference checks pass with `cargo xtask oracle --bless composite-navier-static-thin`. The command recipe uses SI display units explicitly. All geometry and commands are original; no commercial product or deck was used.
