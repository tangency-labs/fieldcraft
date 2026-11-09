# Fatigue (development)

The fatigue kernels and scalar random-stress command path are under development.
The required solver/engine tests, formatting and Clippy checks pass on Linux.
The full oracle retains its 15 existing known failures (750/765 quantities pass).
Release-mode rod examples check coherent stress PSD and moments against an
analytic transfer function, narrow-band damage, and two-reference Miner damage.
This does not establish full-field BASE or published fatigue-example acceptance.

## Conventions

Internal stress is Pa, time seconds, frequency Hz and strain dimensionless.
A range is peak-to-valley; amplitude is half the range. Finite-record rainflow
retains residual half cycles and never closes a record implicitly. Basquin uses
cycles, N = Nref (aref/a)^k. Strain life uses reversals, 2N.
Material fatigue constants must be supplied from applicable measurements; elastic
material constants alone do not define fatigue strength.

The solver library provides Goodman, Gerber and Soderberg corrections and
Coffin–Manson–Basquin inversion with Morrow or Smith–Watson–Topper (SWT).
The chosen policy gives no beneficial credit to compressive mean stress in
Goodman/Gerber/Soderberg/Morrow. SWT requires a positive tensile maximum and
rejects compression-only cycles. This is an explicit scope restriction, not a
prediction of zero damage. Strain-life input must be the local **total** strain
amplitude and the corresponding cycle stress. No elastic-to-plastic conversion
is inferred. Numerical strain-life inversion is bounded to 0.5–1e16 cycles;
out-of-bracket values fail rather than implying an endurance limit.
These library functions do not yet provide a strain-life user command.

## Stored harmonic response to scalar stress PSD

Solve a harmonic step with stress output and a unit spatial loading pattern.
Then extend `solve.random` with a `stress_projection` object:

```json
{
  "step": "Sweep",
  "nodes": [2],
  "input_psd": [[1, 4], [1.02, 4]],
  "stress_projection": {
    "weights": [1, 0, 0, 0, 0, 0],
    "input_model": "common_scalar_unit_pattern",
    "fatigue": {
      "method": "dirlik",
      "curve": {"amplitude_ref": 100000000, "cycles_ref": 1000000,
                "exponent": 3, "endurance_limit": 0},
      "duration_seconds": 60,
      "mean_stress_pa": 0,
      "correction": {"type": "none"},
      "stationary_gaussian": true
    }
  }
}
```

`input_psd` is one-sided, in squared SI excitation units per Hz, on a band
covering every stored frame. The spatial pattern may contain multiple deterministic
loads with prescribed relative phases, but all must be driven by the **same scalar
random input**. The declaration is the user's normalization contract: the program
cannot deduce it from a saved response. Independent or unknown-correlated inputs
are not supported by this interface.

Weights follow xx, yy, zz, xy, xz, yz. The scalar transfer is the complex sum of
weighted components, squared only after summation. A fixed unit plane normal n
uses [nx², ny², nz², 2nxny, 2nxnz, 2nynz]. Neither principal-stress nor von-Mises
spectra can be reconstructed from six component auto-spectra. The legacy summed
`stress_psd` is not used here. Missing stress phasors cause an error.

`projected_stress` returns per-node PSD (Pa²/Hz), fluctuation RMS (Pa), moments
[m0,m1,m2,m4] using Hz, and optional damage. All remain SI regardless of display
units. The piecewise-linear output spectrum is integrated by three-point Gaussian
quadrature per interval; this exactly integrates the moment polynomials, but does
not recover an unresolved resonance. Omit `fatigue` for projection only. The earlier `narrow_band` object remains a narrow-band-only alias; specifying both objects is refused.

The Rayleigh approximation requires a stationary Gaussian scalar linear response,
a narrow spectrum and a single Basquin slope without endurance cutoff. A constant
cycle mean is supplied separately; the PSD does not encode it. The output reports
the irregularity factor and upcrossing rate; it does not certify Gaussianity,
bandwidth, coverage or mesh/modal convergence. Broader spectra need other methods.
Dirlik and Benasciutti–Tovo (2005 four-moment w2 variant) are implemented and exercised by the rod command example. Fatigue contours, full-field BASE validation and transient integration remain pending.

## References

- Stephens et al., *Metal Fatigue in Engineering*, 2nd ed. (2000): stress/strain-life.
- [Altair uniaxial fatigue reference](https://help.altair.com/hwsolvers/os/topics/solvers/os/analysis_fatigue_uniaxial_r.htm): public CMB, Morrow and SWT equations; independently implemented here with the stated domain policy.
- [Dirlik & Benasciutti (2021)](https://doi.org/10.3390/met11091333), equations 1, 2 and 5: spectral moments and Rayleigh damage. Our moments use Hz, so the crossing rate is sqrt(m2/m0) in cycles/s.

## Broad-band and reference-load interfaces (development)

`fatigue.method` selects `narrow_band`, `dirlik`, or `benasciutti_tovo`.
Broad-band methods require zero mean, correction `none`, and exponent 1–20.
They do not replace the distribution of cycle means by a global mean correction.
Invalid mixture coefficients fail; irregularity within 1e-8 of unity uses the
analytical narrow-band limit and reports `narrow_band_limit`.

`solve.fatigue.history` reads current stored linear static steps:

```json
{"channels":[{"step":"UnitX","coefficients":[-1,1,-1]}],
 "times_seconds":[0,1,2],"nodes":[2],"weights":[1,0,0,0,0,0],
 "curve":{"amplitude_ref":100000000,"cycles_ref":1000000,
          "exponent":3,"endurance_limit":0},"correction":{"type":"none"}}
```

Define all reference steps before solving: model edits invalidate stored results.
Channels must share the current mesh and homogeneous supports; nonlinear, contact
and thermal base states are refused. Coefficients multiply the complete reference
pattern, which need not be a one-newton load. Times are SI seconds and match every
channel without interpolation. Six components are summed before projection and
counting; no dynamic inertia is added. Output includes damage, cycle count,
finite representable repeated-record life and ranking of selected nodes only.
Null life is not an endurance certification. Repetition refers to the counted
finite-record spectrum, not periodic rainflow closure. Persistent fatigue contour
fields are still pending.
