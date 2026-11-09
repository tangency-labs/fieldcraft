# mitc6a-navier-harmonic-n12

288 genuine tri6 elements, 625 nodes. Fixed original tolerances. Coarse 32-triangle mesh failed: harmonic components up to 17.4% and Bryan buckling -3.54%; these are not accepted coarse-mesh validations. Same physical loading integrated on the refined shape functions; no tuning of stiffness, mass or tolerance.

The intermediate 128-triangle mesh still failed the 3% band (maximum 4.81%); the 288-triangle mesh gives at most 2.24%. Both real and imaginary components at both frequencies, consistent/lumped and direct/modal12, must pass separately. Probe node 313 is (0.5,0.5,0).
