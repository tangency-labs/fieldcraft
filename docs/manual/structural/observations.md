# Measured modal observations

`observations.import` reads numeric CSV tables with columns `mode,frequency_hz` and optional
`damping_ratio,sigma_hz`. Header order is arbitrary; labels must be unique positive integers.
Damping is a dimensionless fraction, not a percent. Unknown columns, nonfinite numbers, missing
fields, invalid uncertainties and duplicate labels are refused. This command generates frequency observations. `format:"uff"` reads ASCII datasets 15, 55
and 58; 55 normal-mode frequencies become observations, while `raw_uff` preserves coordinates,
coordinate-system IDs, complex global shape components, reference/response directions and FRFs
in their original file units. Unknown datasets are listed. Binary 58b and unsupported 55/58
variants are refused. Dataset 55 complex eigenvalue analysis types 3/7 are not yet supported.
The supported 55 normal-mode vector can itself be real or complex. 58 supports complex single
and double precision with even or uneven abscissa spacing. `format:"accelerance_csv"` (in both
`observations.read` and `observations.import`) reads a measured accelerance table with the header
`frequency_hz,real_accel_mm_s2_per_N,imag_accel_mm_s2_per_N`, nonnegative strictly increasing
frequencies and finite values; it is returned as measured, never converted to receptance and never
turned into solver observations. Conversion of shapes or FRFs to solver observations remains
pending.

```json
{"command":"observations.import","params":{"path":"modes.csv","step":"Modes",
 "mode_indices":[6,7],"sigma_hz":0.1,
 "sensors":[{"id":"tip-z","xyz_m":[1,0,0],"direction_global":[0,0,-1]}],
 "max_distance_m":0.002}}
```

`mode_indices` explicitly assigns each experimental row to a zero-based solver result index.
The command does not guess which rigid-body modes were returned, verify pairing by MAC, or run
an inverse solve. Copy the returned `observations` array to `solve.inverse.joint`. The original
mode labels, damping, line numbers and file SHA256 remain in the reply. File sigma takes
precedence over the fallback argument; `sigma_sources` records which was used. An assumed
weight must not be reported as measured uncertainty or used to claim calibrated confidence.

Sensor positions are explicit metres in the model's global frame, independent of display units.
Directions are signed global unit vectors. Nearest-node distances are in metres, exact ties use
the smaller label, and an excessive distance refuses the whole command. Mapping preserves an
oblique direction as a projection; it does not substitute the dominant axis. Mapping currently
reports geometry only, not whether a model equation is active or constrained. Input is limited
to 16 MiB, 100000 rows, 4096 bytes per line and 10000000 sensor/node distance evaluations.

Frequency-only observations cannot identify homogeneous Young's modulus and density together:
they constrain their ratio. Fix density and use `solve.inverse.frequency`, or choose independent
parameters and let the joint solver check rank. Importing a table is not experimental validation
of the FE model. Real GARTEUR model validation and MAC comparison require matched geometry,
instrument loading, suspension and measured shape vectors.

UFF layout references: UC-SDRL [15](https://www.ceas3.uc.edu/sdrluff/files/15.asc),
[55](https://www.ceas3.uc.edu/sdrluff/files/55.asc),
[58](https://www.ceas3.uc.edu/sdrluff/files/58.asc). No coordinate-system, unit, Fourier-sign
or modal-scale convention is inferred. Dataset 164 is reported as skipped; its units are not
applied. Non-global dataset-15 coordinates cannot be used as global sensor positions without
a verified transformation. Complex shape metadata must not be silently reduced to real values.
