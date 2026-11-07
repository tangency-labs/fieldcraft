# Laminate cantilever strip in pure bending

A 1 × 0.1 × 0.02 m strip uses 8 × 2 MITC4 elements and equal [0/90/90/0] plies. E=[100,10,10] GPa, nu=0, G=[5,3.5,3.5] GPa. The x=0 edge is clamped. The three free-edge nodes take global y moments 0.025, 0.05, 0.025 N m, the consistent edge weights of Mx=1 N m/m.

D11=(710/96)·10^9·h³=59166.6666666667 N m; kappa=1/D11. With the existing shell rotation convention, w(x)=-kappa*x²/2 and sigma_xx=E_x*kappa*z. The reference tip is -8.45070422535211e-6 m. At 90-degree plies the global axial stress is material sigma22, not sigma11. The selected stress fields are from the tip-adjacent element 8, node slot 0, at the four ply midplanes; no stress averaging across plies is used.

Run `cargo xtask oracle --case composite-strip-pure-bending` and freeze only after the reference checks pass with `cargo xtask oracle --bless composite-strip-pure-bending`. The command recipe uses SI display units explicitly. All geometry and commands are original; no commercial product or deck was used.

The recipe also checks all four failure measures against exact rational pure-bending
stresses. With only one axial stress s, maximum stress is |s|/S, Tsai–Hill is
(s/S)² and Tsai–Wu is s(1/St−1/Sc)+s²/(St Sc). Hashin uses the fibre branch
for the outer 0° plies and the matrix branch for the inner 90° plies. The
matrix compression polynomial can be negative at low compression; its raw
value is checked separately and the maximum includes the zero fibre branch.
The references are not obtained from the solved stress output. The predeclared
band is 1e-7 relative plus 1e-12 absolute for structural zeros and rounding.
