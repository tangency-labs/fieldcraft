# Independent SDOF reference (SI, one-sided PSD per Hz)
The recipe assembles a mass=1 kg and spring=4 N/m; only x is free. Modal damping zeta=.02. Unit base inertial load is -M r; solve.random receives base_direction=[1,0,0]. No analytic transfer is passed into the command.
With x=f/fn, a=sqrt(1-zeta²), D=(1-x²)²+4 zeta² x²:
I0(x)=ln((x²+2ax+1)/(x²-2ax+1))/(8a)+(atan((x+a)/zeta)+atan((x-a)/zeta))/(4zeta).
I2 has the negative log term and the same atan term. Their derivatives are 1/D and x²/D.
Force displacement variance on [l,h] is S fn/(m² wn⁴) [I0(h/fn)-I0(l/fn)]. Base absolute acceleration variance is S fn [I0+4 zeta² I2] between limits.
Full-band exact base variance = S fn*pi/(4zeta)*(1+4zeta²); Miles variance omits (1+4zeta²). White case is the finite band 0..100fn: its exact primitive is the primary oracle, not Miles. The report also compares the same FE RMS to Miles with a predeclared approximation allowance; finite-band exact/Miles ratios are listed below. No fixture is derived from a formula compared to itself.

Finite-band exact RMS=3.5383604678701843, infinite-band exact=3.5383612025908269, Miles=3.5355339059327378, finite/Miles relative difference=0.00079947244536482387.

Public Miles formula: NASA MSFC-STD-3676B, Eq. (6), https://standards.nasa.gov/sites/default/files/standards/MSFC/B/0/msfc-std-3676b.pdf . The exact finite-band integral above is independently derived, not a NASA approximation treated as exact.
