# Pressurised pipe loop between two anchors

The loop has legs h = 2000 up from two anchors w = 3000 apart, joined by a run. The pipe has outer radius 50 and
wall 5, E = 210000, ν = 0.3, and internal pressure 10. It is six beam2 elements with rigid corners and no elbow
flexibility.

The closed form, with ε = (1 − 2ν) A/E the pressure's free axial strain:

- ȳ = (h² + w h)/(2h + w) is the height of the elastic centre.
- X (I_c/(E I) + w/(E A) + 2h/(G A_s)) = ε w gives the thrust X, with I_c = 2 (h³/12 + h (h/2 − ȳ)²) + w (h − ȳ)².
- The thrust fixes the four results: the anchor reaction, the corner and crown displacements, and the run's wall force.

The element is exact for this problem (end loads and a uniform initial strain), so the tolerance is round-off.
