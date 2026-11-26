# Measure cylinder shedding independently of the lift spectrum

Follow-up to [#1153](https://github.com/tangency-labs/fieldcraft/issues/1153) and
[#796](https://github.com/tangency-labs/fieldcraft/issues/796). The driver currently identifies
shedding with the strongest lift peak. The Re 3900 seed comparison also uses 100 D/U of data,
shorter than the roughly 156 D/U wake modulation reported by Lehmkuhl et al. An independent velocity
signal and a longer coarse-grid record would make the comparison interpretable.
[Lehmkuhl et al., Fig. 4](https://doi.org/10.1063/1.4818641)

## Suggested change

- At the existing field readback, sample cross-stream velocity at the shear-layer probe (x/D, y/D) =
  (0.71, 0.66) and streamwise velocity at the wake-centreline probe (2, 0), as in Fig. 4. Specify
  the cylinder-centred coordinates, spanwise location, interpolation and actual sample times. Also
  accumulate a spanwise variance of velocity at a fixed wake station to test for three-dimensional
  motion. Save scalar histories, with base pressure and rear-surface recirculation length per
  averaging block.
  [Field samples](https://github.com/tangency-labs/fieldcraft/blob/8802ec4a845556afc5af282004a6b0a75fb04b8a/crates/lbm-gpu/tests/bluff_body_les.rs#L312-L328) ·
  [Rear-surface convention](https://github.com/tangency-labs/fieldcraft/blob/8802ec4a845556afc5af282004a6b0a75fb04b8a/crates/lbm-gpu/tests/bluff_body_les.rs#L350-L379)
- Compute each spectrum using its own cadence. Lift is sampled every 30 steps and fields every 90.
  The current direct Fourier transform scales quadratically with sample count. For long records,
  bound the work by processing declared windows with adequate frequency resolution. Report the bin
  spacing and half-window drift.
  [Cadences](https://github.com/tangency-labs/fieldcraft/blob/8802ec4a845556afc5af282004a6b0a75fb04b8a/crates/lbm-gpu/tests/bluff_body_les.rs#L33-L42) ·
  [Peak estimator](https://github.com/tangency-labs/fieldcraft/blob/8802ec4a845556afc5af282004a6b0a75fb04b8a/crates/lbm-gpu/tests/bluff_body_les.rs#L155-L177)

## Acceptance

- [ ] Two synthetic signals at the force and field cadences recover the same known frequency to
  within one bin. A spanwise-uniform control gives zero variance to roundoff.
- [ ] Run the original D/16 no-ramp pair with and without the probe's 1% seed at initialization.
  Compare lift and probe peaks, spanwise variance and blockwise wake statistics over several
  modulation periods if stable. Record early stops. Compare shedding and rear-surface recirculation
  length with the experimental measurements.
  [Parnaudeau et al.](https://doi.org/10.1063/1.2957018)
