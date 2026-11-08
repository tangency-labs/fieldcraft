# Axial solid on a surface foundation

Unit cube, E=1000 Pa, ν=0, pressure 1 Pa, normal and tangential foundation
100 N/m³ on x=0. Equilibrium and Hooke's law give u_x(x)=1/100+x/1000 m.
Thus u_tip=0.011 m and du_tip/dk=-0.0001 m/(N/m³).
The joint inverse uses noiseless independent analytical observations u(0)=0.008 m
and u(1)=0.0088 m, corresponding to E=1250 Pa and k=125 N/m³. No measured or
commercial results are used. This checks actual commands, not only FEM kernels.
It does not replace finite-foundation beam, Navier plate or BASE acceptance.
