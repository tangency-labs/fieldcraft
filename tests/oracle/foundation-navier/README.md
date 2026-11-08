# Navier plate with a Winkler foundation

The geometry, laminate, simply supported conditions and consistent sinusoidal
load are the original FieldCraft `composite-navier-static-thick` model. The
foundation acts normally on the shell reference surface with k=1e8 N/m³.

For the separated FSDT system, adding k W to the transverse equilibrium changes
only K33. Eliminating the two rotations gives W=1/(1/W0+k), where the parent
case has W0=1.1636747564990154e-9 m at unit pressure. The mesh acceptance remains
3%, fixed before the run. This tests composite-shell foundation integration
through public commands. It is not a measured or commercial reference.
