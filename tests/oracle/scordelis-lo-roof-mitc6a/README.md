# scordelis-lo-roof-mitc6a

Original six-node triangular mesh, obtained by splitting the original quadrilateral geometry and adding true quadratic edge midpoints. No hidden central DOFs. 289 nodes, 128 elements. Same continuum data and tolerances as the corresponding MITC8 reference case. Sinusoidal loads, when used, are integrated against triangular quadratic shapes. Curved roof diagonal midpoints are on the exact cylinder.

Quarter cylinder: R=25, full axial length 50, half-opening 40 degrees,
h=0.25, E=4.32e8, nu=0, volume force [0,0,-360]. The x=25 diaphragm
holds uy, uz, rx; x=0 and y=0 are symmetry planes. Original free-edge
probe node 209 is at x=0, phi=40 degrees. Reference global uz=-0.3024,
fixed 2% tolerance. Shared quadratic-surface nodal directors are used.
This is a curved six-node shell model, not a flat-facet substitution.
