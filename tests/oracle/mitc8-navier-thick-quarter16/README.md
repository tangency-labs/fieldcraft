# Thick Navier plate: quarter model

Full plate 1 x 1 m, h=0.1 m, equal [0/90/90/0] plies, E=[100,10,10] GPa,
nu=0, G=[5,3.5,3.5] GPa. Global x/y define the ply reference axes. The model
covers 0<=x,y<=0.5 with 16 x 16 quad8 elements (833 nodes), equivalent by symmetry
to a 32 x 32 full mesh. It is not the earlier 4 x 4 model.

Pressure q=sin(pi*x)sin(pi*y) acts in +z. Consistent four-point Gauss integration
per surface direction gives total quarter load 0.1013211836423377 N, compared with
1/pi^2. Hold w on x=0 or y=0, rx on x=0 or y=0.5, ry on y=0 or x=0.5.
In-plane translations and drilling are held: B=0 makes this compatible with the
separated flexural solution. Probe x=y=0.5, the full plate centre.

The independent FSDT reference is 1.1636747564990154e-9 m, with the fixed 5/6 shear
factor; see ../composite-navier-frequency-thick/README.md for the 3 x 3 separated
Navier matrix and rational thickness moments. This is not a three-dimensional
or equilibrated interlaminar-stress reference.

| Full-mesh resolution | Centre w (m, rounded) | Reference error | Fixed 3% gate |
|---|---|---|---|
| 4 x 4 | 1.3286e-9 | 14.2% | fail |
| 8 x 8 | 1.2652e-9 | 8.72% | fail |
| 16 x 16 | 1.2035e-9 | 3.43% | fail |
| 32 x 32 (quarter model here) | 1.1753e-9 | 0.999% | pass |

The coarser results are not substituted by the finer result or declared passing.
The original BD86 mixed strain interpolation is retained; no shear penalty or
reference tolerance was fitted. This observed thick-plate over-flexibility decreases
with refinement. Only the stated fine mesh is qualified by this oracle; a coarse
thick MITC8 mesh needs its own accuracy check. The independent stiffness rank test
and rigid-motion test check zero modes separately from this displacement accuracy.
