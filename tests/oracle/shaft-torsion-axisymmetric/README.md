# Solid circular shaft in torsion, axisymmetric with twist (Saint-Venant)

A solid shaft of radius 20 and length 100 as an axisymmetric slice r ∈ [0, 20], 0 ≤ z ≤ 100 (y is the
axis, with the twist carried by the model's Z displacement u_z), one quad4 patch of 4 × 10 elements. The bottom end is
held (all degrees of freedom fixed) and the top end carries a torque T = 10⁶ N·mm spread over the end in
proportion to the radius.

The closed form (Saint-Venant, Timoshenko and Goodier, ch. 10) is u_z = −T r z / (G J), τ_zθ (the yz component) = −T r / J and
τ_rθ = 0, with G = E / (2(1 + ν)) = 80769.23 MPa and J = π R⁴ / 2 = 251327.41 mm⁴. At r = 20, z = 100 the
z displacement is −0.0985245 mm and at r = 20, z = 50 the shear stress τ_zθ is −79.577 MPa.

The tolerances are round-off: the bilinear element holds u_z = −κ r z exactly on a rectangle and its 2 × 2
rule integrates the forces of τ_zθ exactly, so the solve reproduces the closed form to round-off.
