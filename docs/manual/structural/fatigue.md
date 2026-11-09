# Fatigue (development)

The fatigue kernels and scalar random-stress command path are under development.
The required solver/engine tests, formatting and Clippy checks pass on Linux.
The full oracle retains its 15 existing known failures (750/765 quantities pass).
Release-mode rod examples check coherent stress PSD and moments against an
analytic transfer function, narrow-band damage, and two-reference Miner damage.
This does not establish full-field BASE acceptance or real material fatigue life.

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
`solve.fatigue.strain_cycles` exposes these laws for explicitly counted local
total-strain cycles; it does not infer plastic strain from elastic stress.

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
a narrow spectrum. Its legacy closed form uses single-slope Basquin without
endurance cutoff; the table path integrates piecewise laws and a cutoff. A constant
cycle mean is supplied separately; the PSD does not encode it. The output reports
the irregularity factor and upcrossing rate; it does not certify Gaussianity,
bandwidth, coverage or mesh/modal convergence. Broader spectra need other methods.
Dirlik and Benasciutti–Tovo (2005 four-moment w2 variant) are implemented and exercised by the rod command example. Native fatigue contours and the exported SI transient adapter are available; full-field BASE validation remains open.

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
finite representable repeated-record life and ranking of the requested nodes.
Null life is not an endurance certification. Repetition refers to the counted
finite-record spectrum, not periodic rainflow closure. For full node coverage and
stored contours, see the output instructions below.


## Exported transient stress history

`solve.fatigue.transient` consumes the actual exported `solve.probe` JSON as
`history`. It requires `spectrum:false`, explicit units `stress:Pa`, `time:s`,
the component order `[xx,yy,zz,xy,xz,yz]`, and the current exporter location
`global nodal average; shell top surface`. Missing/null stresses are errors.
Provide `step`, the exact `nodes` set, a fixed `weights` vector, `curve`,
`correction`, and `max_time_step_seconds`. All sample times must match the frame
axis exactly, increase strictly, and satisfy the declared maximum gap. No peaks
are reconstructed between samples; a small gap is not proof of adequate sampling.

`source_id` retains the caller-verified model/export identity. Labels refer only
to that external source, never to the current session mesh. The file-importing
workflow must verify the original model and raw export hashes before this command;
the command does not authenticate a caller's identity string. The integration
uses the exported interface without merging an unpublished transient solver branch.

Results include per-node rainflow cycles (range, amplitude, mean, count), per-cycle
Miner increments and cumulative damage in counting order, not a chronological
crack-growth history. Residual half cycles are retained; the supplied record is not
closed into a period. Repeated-record life repeats the counted spectrum only.
The initial handoff is a two-node rod with 1001 samples; it is not BASE evidence.
Linux engine tests and refusal cases pass. The original 1001-frame rod export
passes release-mode checks against its analytic step-response stress and one
residual half cycle, with uncorrected and Goodman Miner damage. This verifies
the exported-data path, not a native transient solver or BASE fatigue life.


## BASE demonstration scope

A release-mode face-supported BASE case (314430 nodes, 198174 Tet10 elements)
uses 40 modes, six residual vectors and 71 frequency points over 0–2000 Hz.
Two probes bind to the actual source mesh. Narrow-band, Dirlik and BT results
share the same scalar stress projection; independent moment integration and
narrow-band damage checks pass. These checks validate the command path and
postprocessing of the sampled spectrum, not frequency-grid or modal convergence.
The S–N data are synthetic. The ranked probes are not a full-field hotspot search,
and their reported life is not an engineering prediction for the real component.


## Counted total-strain cycles

```json
{"curve":{"young_pa":200000000000,"strength_pa":1000000000,
          "ductility":0.5,"b":-0.5,"c":-0.5},
 "cycles":[{"total_strain_amplitude":0.00505,"count":2.5}],
 "correction":"none"}
```

Pass this object to `solve.fatigue.strain_cycles`. It reports cycles-to-failure N,
per-cycle damage and cumulative Miner damage (here N=5000, D=0.0005).
For `morrow`, each cycle must also contain `mean_stress_pa`; for `swt`, each must
contain positive `maximum_stress_pa`. Supply stress and total strain for the same
local cycle. No stress-to-plastic-strain conversion, automatic stress/strain cycle
pairing, or life in seconds is inferred. Zero amplitude produces zero damage and
null life, not a validated fatigue limit. The existing 0.5–1e16 cycle bracket applies.

## Delivered scope and remaining limits

Time-domain commands accept the legacy Basquin curve or a segmented S–N table:

```json
{"points":[{"cycles":100,"amplitude":400000000},
           {"cycles":10000,"amplitude":40000000},
           {"cycles":10000000,"amplitude":4000000}],
 "endurance_limit":2000000,"surface_factor":0.8,"size_factor":0.5,
 "extrapolation":"reject"}
```

Amplitudes and endurance limit are fully reversed Pa before reduction. Cycles must
strictly increase and amplitudes decrease (2–1024 points). Interpolate linearly in
log amplitude / log cycles; `reject` refuses amplitudes outside the table, while
`end_slopes` explicitly extrapolates the first/last segment. Amplitudes at or below
the reduced endurance limit cause zero damage. Surface and size factors are
user-supplied strength multipliers in (0,1], not empirical estimates: after mean
correction divide amplitude by their product before evaluating the unreduced
curve. They also reduce the endurance limit; they do not change ultimate/yield
strength in the mean correction. Two points can represent a single Basquin slope.
Spectral tables are supported by narrow-band and Dirlik as detailed below; TB
retains only its single-slope Basquin interface.

`solve.fatigue.history` accepts `nodes:"all"` or explicit labels and saves the field
on its first reference step. `solve.random` with harmonic frames saves projected
fatigue fields on the harmonic step; omit `nodes` for all nodes. Existing bounded
work limits apply (history: 20 million node-time-channel products; spectral:
4 million node-frame samples). A successful subset is never called full coverage.
Both commands report coverage and descending damage rankings with label tie-breaks.

Use `results.contour {"step":"A","field":"fatigue_damage"}` or
`fatigue_life` for native nodal contours. Damage is dimensionless; life is seconds
in storage, converted to display time units. Missing coverage and zero-damage
unbounded life have no colour value. These accumulated fields cannot be animated
as proportional stress. `file.save` / `file.open` preserve fields, exact command
inputs and all reference definition hashes; changing any reference invalidates
the result, including a secondary superposition channel. New result kinds 6/7
carry the extension; old static/harmonic binary layouts 0/5 remain readable.
`results.export.vtu` writes damage/life plus `has_fatigue_damage` and
`has_fatigue_life` masks. Masked zeros are placeholders, not computed values.
Harmonic VTU output contains fatigue fields, not complex phasor vectors.
External transient exports retain their source labels and JSON cycle results;
they are not attached to an unrelated current mesh.

Verification includes the seven individual finite-record cycles of ASTM
E1049-85(2017) §5.4.4.2 / Fig.6, and the original authors' numerical Table 1
in Dirlik & Benasciutti (2021), DOI 10.3390/met11091333: rectangular spectral
band ratios 5 and 20, slopes 3/5/7, published TB/Dirlik damage ratios. The 0.001
absolute comparison tolerance covers printed three-decimal precision. These are
numerical benchmarks, not material validation or an experiment.

Stress histories support fixed projection and proportional-tensor signed von Mises /
maximum-principal reductions below. General critical-plane fatigue remains open.
Persistent model-level fatigue material records are not implemented. The
strain-cycle command consumes already-counted cycles; transient integration
consumes exported SI probe data. The BASE demonstration above remains a
**two-candidate, synthetic-S–N demonstration**, not a full-field hotspot assessment.


### Complete tensor history criteria

`solve.fatigue.history` and `solve.fatigue.transient` accept **either** legacy
`weights` **or** `criterion: {"type":"signed_von_mises"}` /
`criterion: {"type":"max_principal"}`. Explicit projection is
`criterion: {"type":"projection","weights":[1,0,0,0,0,0]}`.
Tensor order is `[xx,yy,zz,xy,xz,yz]`, global coordinates, Pa, with tensor shear
stresses (no engineering-strain factor). Reference-load tensors are summed at
each time before scalarization; imported transient samples must contain all six
components with the documented source identity, location and units. Times are s.

Signed von Mises is sqrt(3 J2) times the sign of the principal stress with largest
absolute value, following [Altair's public equivalent nominal stress convention](https://help.altair.com/hwsolvers/os/topics/solvers/os/analysis_fatigue_uniaxial_r.htm).
Opposite-sign ties within relative 1e-10 (including pure shear) are rejected:
a sign cannot be inferred. Hydrostatic stress gives zero von Mises.
`max_principal` is the **algebraically largest** eigenvalue, not the signed
absolute-largest principal. For uniaxial compression its value is zero; thus
`[-2,2,-2]` Pa uniaxial history becomes `[0,2,0]` Pa under this criterion.
Eigenvalues use the existing Smith (1961), CACM 4(4):168 symmetric-tensor routine.

Both nonlinear reductions are deliberately limited to proportional tensors
sigma(t)=lambda(t) sigma0, including reversals and zero samples. Every nonzero
sample is normalized by its largest absolute component and checked against the
first nonzero reference (absolute normalized residual <=1e-9). A rotating or
non-collinear prestressed tensor history is rejected. This is a nominal uniaxial
S–N approximation requiring a material curve justified for that scalar measure,
not a general multiaxial fatigue or plastic strain model. A justified fixed
projection remains available; the caller must ensure temporal extrema are sampled.

Rainflow acts on the resulting scalar samples. Cycle range=max-min,
amplitude=range/2, mean=(max+min)/2; mean correction acts on each scalar cycle,
then Miner sums count/N. No tensor-mean substitution is made. Residual cycles
retain count=0.5. `solve.random` refuses nonlinear criteria: component PSDs do not
identify signs, phases, joint tensor histories or nonlinear principal-stress PSDs.
Even a cross-spectral matrix does not by itself specify a realized sign history.

Analytic regression: sigma0=[3,0,0,2,0,0] Pa has eigenvalues (4,0,-1),
VM=sqrt(21). Its 45-degree rotated representation [-0.5,3.5,0,1.5,0,0]
gives identical results. For lambda=[-1,1,-1], the two residual half cycles
sum to one cycle. Synthetic N=100/a² gives signed VM D=0.21 (mean zero).
Max principal samples [1,4,1] give a=1.5, mean=2.5; Goodman ultimate=10 Pa
corrects a to 2 Pa and gives D=0.04. Uniaxial sigma0=2 Pa gives D=0.04
(signed VM) or 0.01 (max principal, no correction). These are analytic
implementation checks, not measured material data or BASE validation.

### Segmented S–N in spectral fatigue

`solve.random` → `stress_projection.fatigue.curve` accepts the same S–N table
as time-domain fatigue: `points:[{cycles,amplitude}]` in cycles/Pa,
`endurance_limit` in Pa, `surface_factor`, `size_factor` in (0,1], and
`extrapolation:"end_slopes"`. Tables require 2..64 points and each log-log
slope in 1..20. Node count times point count is capped at 256 per command
to bound quadrature work. Strength factors multiply the S–N stress axis and endurance
limit; the applied PSD is unchanged. A factor below one increases damage.

`narrow_band` integrates the Rayleigh amplitude density with the zero-crossing
rate. `dirlik` integrates its exponential/two-Rayleigh mixture with the peak
rate; in its near-narrow-band limit it uses Rayleigh and the crossing rate.
Both split the integral at every table knot and endurance cutoff. The public
basis is Dirlik & Benasciutti (2021), Metals 11:1333, equations (4), (7), (8):
D = duration × cycle-rate × integral[p(amplitude)/N(amplitude)]. The distribution
is over amplitude, half the rainflow range; no stress history is constructed.

The numerical implementation integrates each power segment after t=u² using
bounded Simpson refinement (relative successive-integral tolerance 2e-9).
The exponential tail is bounded by a cutoff beyond the mode with an additional
128+20 ln(1+lower) margin for supported slopes. Out-of-domain intervals,
nonconvergence and nonrepresentable damage return errors, never a silent zero.
The table's outer slopes must be explicitly extrapolated: Gaussian amplitude
support is unbounded, so `extrapolation:"reject"` is refused by name.

`benasciutti_tovo` with a table is refused by name: the existing exponent-dependent
single-slope correction has no implemented piecewise extension. Dirlik tables
require zero mean and `correction:{type:"none"}`; narrow-band tables retain the
existing constant scalar cycle-mean correction. Legacy single-slope Basquin
behavior is unchanged. Missing correlations, nonlinear invariant PSDs, and
PSD-derived signed time histories remain unsupported. These are approximate
stationary Gaussian scalar stress methods, not general multiaxial laws.

This extends the former single-slope-only interface; it does not claim measured
material properties or full-field BASE validation. Existing BASE evidence remains
two candidate points with synthetic material data.
