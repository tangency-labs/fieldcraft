# A cantilever of a general beam profile (#892)

Two cantilevers of four beam2 elements each, L = 500 mm along X, clamped at x = 0, of one general
profile given only by its section properties: A = 2000 mm², I_y = 4·10⁶ mm⁴, I_z = 10⁶ mm⁴,
J = 5·10⁴ mm⁴, A_s,y = 800 mm², A_s,z = 1200 mm² (no shape has exactly these numbers; they meet
the bounds every section obeys). The orientation is +Y, so local y is the width (bending across
it takes I_z and A_s,y) and local z the height (I_y and A_s,z). Steel, E = 210 000 MPa, ν = 0.3,
G = E/(2(1 + ν)).

The first cantilever is a Timoshenko beam, the second the same profile with
`"beam_theory": "euler_bernoulli"`. Both carry the tip force (1000, 2000, −3000) N and the tip
torque 5·10⁵ N mm about X. The references are the closed forms for a cantilever under end loads:
stretch N L/(E A), deflections P L³/(3 E I) plus, for the Timoshenko beam only, the shear
deflection P L/(G A_s), and twist T L/(G J). The shear terms are 3.8 % (width) and 9.4 % (height)
of the Timoshenko deflections, so the case also tells the two beam theories apart.

Both elements reproduce end-loaded cantilevers exactly at the nodes, so the mesh needs no
refinement and the tolerance is round-off (10⁻⁹).
