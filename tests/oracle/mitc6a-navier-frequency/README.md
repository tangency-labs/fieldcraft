# mitc6a-navier-frequency

Original six-node triangular mesh, obtained by splitting the original quadrilateral geometry and adding true quadratic edge midpoints. No hidden central DOFs. 81 nodes, 32 elements. Same continuum data and tolerances as the corresponding MITC8 reference case. Sinusoidal loads, when used, are integrated against triangular quadratic shapes. Curved roof diagonal midpoints are on the exact cylinder.

Unit square, total h=0.01, symmetric [0/90/90/0], density 1600.
Engineering constants and support definitions match the original
[MITC8 Navier frequency model](../mitc8-navier-frequency/README.md).
The reference is the same independent 3x3 FSDT Navier eigenproblem;
membrane and drilling DOFs are held to isolate the symmetric flexural
branch. First three errors must satisfy 2%, 5%, 5% separately for each
mass option. This does not certify arbitrary coarse laminate spectra.
