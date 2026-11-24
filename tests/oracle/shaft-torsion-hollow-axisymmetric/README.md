# Hollow circular shaft in torsion, axisymmetric with twist (Saint-Venant)

A hollow shaft with bore radius 10 and outer radius 20, length 100, as an axisymmetric slice r ∈ [10, 20],
0 ≤ z ≤ 100 (y is the axis, with the twist carried by the model's Z displacement u_z), one quad8 patch of 2 × 5 elements.
The bottom end is held and the top end carries a torque T = 10⁶ N·mm spread over the end in proportion to
the radius.

The closed form (Saint-Venant, Timoshenko and Goodier, ch. 10) is u_z = −T r z / (G J), τ_zθ (the yz component) = −T r / J and
τ_rθ = 0, with G = 80769.23 MPa and J = π (b⁴ − a⁴) / 2 = 235619.45 mm⁴. At r = 20, z = 100 the hoop
displacement u_z is −0.1050928 mm and at r = 15, z = 40 the shear stress τ_zθ is −63.662 MPa.

The tolerances are round-off: the eight-node element with straight sides holds u_z = −κ r z exactly and its
3 × 3 rule integrates its forces exactly, so the solve reproduces the closed form to round-off.
