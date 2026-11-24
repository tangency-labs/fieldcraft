# Uniform tension through a distributing coupling

Issue #1024. A 100 × 10 × 10 mm steel bar (E = 210000 MPa, ν = 0.3) of 20 × 2 × 2 hex20 elements,
three symmetry planes, and a 1000 N pull on a reference point at the centre of its end face. The
distributing coupling ties the end face's element faces to the point, so the point's force reaches
the nodes as the consistent loads of a uniform plus linear traction.

Closed form: u = F L/(E A) = 1000 · 100/(210000 · 100) = 0.0047619 mm everywhere on the bar, and
σ = F/A = 10 MPa. The case passes exactly: the nodal forces of a uniform traction are the consistent
loads hex20 integrates exactly, so only round-off remains. Before the fit was taken over the faces,
the same bar was 3.5 % long at the point and reached 43 MPa at the end-face corners.
