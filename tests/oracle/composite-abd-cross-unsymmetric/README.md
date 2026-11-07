# Laminate ABD: cross-unsymmetric

Bottom-to-top angles are [0, 90] degrees. This algebraic constitutive case uses h=1 m, equal ply thicknesses, E=[100,10,10] Pa, nu=0, G=[5,4,3] Pa and density 2 kg/m³. The reference matrices are rational numbers obtained by integrating Qbar, z Qbar and z² Qbar analytically, independent of the FieldCraft implementation.

For [0/90/90/0], the outer pair has ∫z²dz=7/96 and the inner pair 1/96, giving D11=710/96 and D22=170/96. For [0/90], B11=-45/4 and B22=45/4. In [45/-45/-45/45], A16=A26=B=0 but **D16=D26=45/32**, so balanced symmetry must not remove bending-twisting coupling.

All nine entries of each A/B/D block, both shear directions and all mass moments are checked. These are constitutive checks, not a shell solve or experimental validation. Source: [A. T. Nettles, Basic Mechanics of Laminated Composite Plates, NASA RP-1351 (1994), §§II–IV](https://ntrs.nasa.gov/citations/19950009349).

Run `cargo xtask oracle --case composite-abd-cross-unsymmetric` and freeze only after the reference checks pass with `cargo xtask oracle --bless composite-abd-cross-unsymmetric`. The command recipe uses SI display units explicitly. All geometry and commands are original; no commercial product or deck was used.
