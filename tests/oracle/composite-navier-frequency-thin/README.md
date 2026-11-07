# Navier laminate frequencies: thin

A 1 × 1 m square, h=0.01 m, uses 16 × 16 MITC4 elements and equal [0/90/90/0] plies, with E=[100,10,10] GPa, nu=0, G=[5,3.5,3.5] GPa and rho=1600 kg/m³. The reference frequencies are 40.8226067142, 80.7234365986, 139.101797643 Hz. All three transverse mode branches are compared; the first three modes are (m,n)=(1,1),(1,2),(2,1).

Hard simply supported edges hold w on all edges, rotation about x on x=0,1, and rotation about y on y=0,1. In-plane translations and drilling rotations are held throughout; for this B=0 symmetric laminate these constraints isolate the flexural branch without changing its frequencies.

The references follow a separated Navier mode: w=W sin(alpha x) sin(beta y), with alpha=m*pi and beta=n*pi. Let the rotation amplitudes follow the corresponding cosine/sine functions; changing both rotation signs does not change the eigenvalues. The symmetric stiffness and mass matrices are

```text
K11 = D11 alpha² + D66 beta² + S
K22 = D66 alpha² + D22 beta² + S
K12 = (D12 + D66) alpha beta
K13 = S alpha; K23 = S beta; K33 = S(alpha²+beta²)
M = diag(I2,I2,I0)
D11=(710/96)*10^9*h³; D22=(170/96)*10^9*h³
D12=0; D66=(5/12)*10^9*h³; S=(5/6)*3.5*10^9*h
I0=rho*h; I2=rho*h³/12
```

The smallest positive root of det(K-lambda M)=0 gives f=sqrt(lambda)/(2*pi). Reference arithmetic uses these rational thickness moments, not the FEM laminate implementation. The listed values were computed by 80-digit bounded bisection of the cubic and independently placed in the smallest positive root intervals by Jacobian polynomial root isolation (pi represented by a 50-digit decimal rational; exact isolation applies to that rational input). The numerical reference is approximate to pi, not a claim of exact transcendental arithmetic.

The 2% / 5% acceptance bands are fixed mesh-error allowances; they must not be widened to hide a failed mode. Source: [C. W. Bert and C. C. Siu, NASA CR-112140 (1972), §§2.2–2.6 and 4.1](https://ntrs.nasa.gov/citations/19720022240). This is validation against the chosen equivalent-single-layer theory, not three-dimensional elasticity.

Run `cargo xtask oracle --case composite-navier-frequency-thin` and freeze only after the reference checks pass with `cargo xtask oracle --bless composite-navier-frequency-thin`. The command recipe uses SI display units explicitly. All geometry and commands are original; no commercial product or deck was used.
