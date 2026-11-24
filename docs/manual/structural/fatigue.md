# Fatigue (stress-life)

Stress-life (S-N) fatigue is post-processing of a solved linear static step. The step's stress
field is the response to one reference load; a fatigue assessment says how that load repeats (a
constant amplitude cycle, a stress ratio, or a history of load factors counted by rainflow), and
each material's S-N data says how many cycles it survives. At every node FieldCraft reads a
signed scalar stress, corrects each cycle for its mean stress, reads the life off the S-N curve,
and sums the damage with the Palmgren–Miner rule: the result is the life (cycles, or blocks of a
history), the damage per cycle or block, and the safety factor against the endurance limit, as
fields for contours, probes and path plots, and a summary of where the life is least.

| | |
|---|---|
| **Commands** | `model.material.fatigue` (S-N data), `model.fatigue.set` and `model.fatigue.delete` (the loading), `results.fatigue` (summary); fields `fatigue_log_life`, `fatigue_damage`, `fatigue_safety_factor` in `results.contour`, `solve.probe` and `results.path` |
| **Step procedure** | `static` (linear) |
| **Solver** | post-processing in `crates/fem/src/fatigue.rs` (the method) and `crates/engine/src/fatigue.rs` (the fields), from the step's stored nodal stresses |
| **Method** | Basquin's S-N relation (Basquin 1910; Dowling 2013, ch. 9); Goodman (1899), Gerber (1874), Soderberg (1930) and Smith–Watson–Topper (1970) mean stress corrections (Dowling 2013, ch. 9; Shigley); rainflow counting (Matsuishi & Endo 1968; ASTM E1049-85; Downing & Socie 1982); Palmgren–Miner damage (Miner 1945) |
| **Verified by** | `crates/fem/src/fatigue.rs` tests; `crates/engine/src/cmd/fatigue/tests.rs` |

## Synopsis

```json
[
  {"command": "model.material.fatigue", "params": {"name": "Steel", "fatigue": {
    "curve": {"type": "basquin", "coefficient": 1000, "exponent": -0.1},
    "endurance_limit": 150, "ultimate_strength": 800, "yield_strength": 600}}},
  {"command": "solve.static", "params": {"step": "Load"}},
  {"command": "model.fatigue.set", "params": {"step": "Load",
    "loading": {"type": "stress_ratio", "ratio": -1},
    "measure": "signed_von_mises", "mean_stress": "goodman"}},
  {"command": "results.fatigue", "params": {"step": "Load"}},
  {"command": "results.contour", "params": {"step": "Load", "field": "fatigue_log_life"}}
]
```

The numbers above are an illustration in MPa, not a real material's data. **The S-N data is
yours**: FieldCraft ships no fatigue constants. Neither the data nor the assessment is read by a
solver, so setting or changing them keeps the step's result current: no solve is needed.

## Parameters

### Fatigue data (`model.material.fatigue`, or `fatigue` in `model.material.create`; `null` removes it)

The curve is the stress amplitude of fully reversed cycles (stress ratio R = −1) that smooth
specimens of the material survive for N cycles.

| Field | Value | Default | Meaning |
|---|---|---|---|
| `curve` | `{"type": "basquin", "coefficient": σf′, "exponent": b}` | — | σa = σf′ (2N)^b, N in cycles (2N reversals); σf′ above 0, −1 < b < 0 |
| | `{"type": "table", "points": [{"cycles", "amplitude"}, …]}` | — | 2 to 64 points, cycles rising from 0.5 and amplitudes falling, joined by straight lines on log-log axes |
| `endurance_limit` | stress above 0, below the curve's amplitude at its short-life end (Basquin: σf′, at one reversal; table: the first point) and, for a table, not above its last point (that would cut off part of the curve: refused) | none | the fully reversed amplitude at and below which the life is infinite; the safety factor is against it |
| `ultimate_strength` | stress above 0 | none | S_u: the Goodman and Gerber corrections (SWT needs no strength) |
| `yield_strength` | stress, at most S_u | none | S_y: the Soderberg correction, and the check for yield |

Basquin's relation is stored in reversals, as the strain-life literature writes it (Dowling
2013, ch. 9). The same curve in cycles is S = A N^b with **A = σf′ 2^b**: a curve given as
S = A N^b converts by σf′ = A / 2^b, with the same b.

### The loading (`model.fatigue.set`; one per step, replaced by the next)

| Field | Value | Default | Meaning |
|---|---|---|---|
| `step` | a linear static step | — | its stress field is the response to the reference load |
| `loading` | `{"type": "constant_amplitude", "max", "min"}` | — | one cycle repeated, from `min` to `max` times the step's load (`max` ≥ `min`) |
| | `{"type": "stress_ratio", "ratio": R, "scale"?: 1}` | — | one cycle from `scale` to R · `scale` times the step's load, −1 ≤ R ≤ 1 |
| | `{"type": "history", "factors": […], "counting"?: "repeated"}` | — | one block of 2 to 100 000 load factors, rainflow counted: `repeated` (the block repeats; every cycle closes) or `once` (the history happens once; ranges left at the end count as half cycles) |
| `measure` | `signed_von_mises`, `max_abs_principal` | — | the scalar stress read at each node ([Stress measures](#stress-measures)) |
| `mean_stress` | `none`, `goodman`, `gerber`, `soderberg`, `swt` | — | the mean stress correction ([Mean stress](#mean-stress)) |
| `stress_factor` | above 0, at most 100 | 1 | multiplies the whole stress history, amplitude and mean, before the curve and the strengths are read |

Load factors are at most 10⁹ in magnitude. The method and the correction have no default: they
are the engineer's choice.

**The stress factor** is one multiplier of the whole stress history, amplitude and mean alike
(the step is linear, so a factor on the load is a factor on every stress and on the strength
checks): a stress concentration the mesh does not resolve, or 1/k for a surface finish or size
factor k of the part against the specimens of the curve (Shigley, on the fatigue stress
concentration, surface and size factors). It is one number for both uses, so choose it for the
one you mean: a fatigue notch factor is usually applied to the amplitude and, for a notch that
does not yield, to the mean too, and a surface factor is applied to the strength. A notch the mesh
resolves already has its stress concentration in the solved stress: do not count it twice.

### `results.fatigue`

| Field | Value | Default | Meaning |
|---|---|---|---|
| `step` | a solved step with an assessment | the step shown | |
| `set` | set name | every node | summarise over the set's nodes only, for example to leave out the singular stress at a point support |

## Theory

### Proportional loading

The step is linear, so its stress is proportional to its load. The assessment scales the whole
field by load factors λ(t): every stress component at a node follows the same history, σ(t) =
λ(t) σ_ref, and the principal directions stay fixed. This is **proportional loading**, the case
for which a scalar equivalent stress and uniaxial S-N data apply (Bannantine, Comer & Handrock
1990, ch. 5). Loads that vary out of phase with each other (two load cases with their own
histories) are non-proportional: they are not expressible as one step's factors, and need a
multiaxial (critical-plane) method this page does not provide.

### Stress measures

From the reference stress σ_ref at a node:

- `signed_von_mises`: the von Mises stress σ_vm, with the sign of the principal stress of largest
  magnitude (for ductile metals, whose yield and fatigue follow the distortion energy);
- `max_abs_principal`: the principal stress of largest magnitude, with its sign (for brittle
  materials, and cracks that grow normal to the largest principal stress).

**The sign rule.** A state counts as tensile when its largest tensile principal stress is at
least 99.9 % of the largest compressive one in magnitude (`SIGN_TIE_BAND` = 10⁻³ in
`crates/fem/src/fatigue.rs`), and so does any state whose largest-magnitude principal stress is
tensile. Pure shear sits exactly on the tie, and a nodal stress never lands exactly there: with a
strict comparison, rounding noise (1 µMPa in one component of a shear of 200 MPa) flipped the
sign and moved a Goodman life by a factor of 11.5, and under SWT a node read as compressive from
zero load has σ_max = 0 and no damage. Tension is the conservative side for every correction here.
A state with no tensile principal stress is compressive however close its two largest are. The
rule moves the discontinuity, it does not remove it: a node just outside the band (0.1 %
compressive) still reads the opposite sign, so a contour near a tie can still show a step, and a
shear-dominated region is read with that sign convention in mind. The sign is a
convention that carries tension and compression into the mean stress: under proportional loading
the measure of λ σ_ref is λ times the measure of σ_ref, so the node's stress history is m λ(t)
with m the measure of σ_ref. The von Mises stress itself is never negative; the sign decides only
whether a mean is tensile or compressive.

### Cycles

A constant amplitude cycle between factors λ_max and λ_min has the range Δλ = λ_max − λ_min and
the mean λ̄ = (λ_max + λ_min)/2. A history is counted by **rainflow** (Matsuishi & Endo 1968) on
its peaks and valleys:

- `once`: the three-point procedure of ASTM E1049-85. With X the most recent range and Y the
  one before it, Y is counted when X ≥ Y: as one cycle, or as half a cycle when Y holds the
  history's first point (which then moves on). The ranges left at the end are half cycles.
- `repeated`: Downing and Socie's method (1982) for a history that repeats: the block is
  rotated to start at its value of largest magnitude and closed back to it, so the same three-point
  rule closes every cycle, the one across the join of two blocks included.

Counting looks only at ranges, so the cycles of m λ(t) are those of λ(t) with ranges times |m|
and means times m: the history is counted once, and each node reads its own m. At a node the
cycle's stress amplitude and mean are

```text
σa = K · |m| Δλ / 2,   σm = m λ̄            (K the stress factor)
```

### Mean stress

A tensile mean shortens the life. The curve is read at the equivalent fully reversed amplitude
σ_ar (Dowling 2013, ch. 9; Bannantine et al. 1990, ch. 2):

| `mean_stress` | σ_ar for σm > 0 | Safety factor n against S_e (Shigley) |
|---|---|---|
| `none` | σa | S_e / σa |
| `goodman` | σa / (1 − σm/S_u) | 1/n = σa/S_e + σm/S_u |
| `gerber` | σa / (1 − (σm/S_u)²) | n = ½ (S_u/σm)² (σa/S_e) [−1 + √(1 + (2 σm S_e/(S_u σa))²)] |
| `soderberg` | σa / (1 − σm/S_y) | 1/n = σa/S_e + σm/S_y |
| `swt` | √(σ_max σa), σ_max = σa + σm | n = S_e / √(σ_max σa) |

For Goodman, Gerber and Soderberg, a compressive mean (σm ≤ 0) is given no correction (σ_ar = σa,
n = S_e/σa): the lines are drawn through tensile means, and a benefit from compression is not
taken (Shigley). **SWT** (Smith, Watson and Topper 1970, in the stress-life form of Dowling
2013, ch. 9) is a different kind of rule: the damage parameter is the product of the cycle's
largest stress σ_max = σa + σm and its amplitude, so the curve is read at the fully reversed
amplitude with the same product, √(σ_max σa). It needs no strength, it is also read at a
compressive mean (a mean below zero lowers σ_max and so the parameter), a fully reversed cycle is
read at its own amplitude, and a cycle whose σ_max is not tensile (σm ≤ −σa) does no damage
(σ_ar = 0). For the other three, a mean at or
past the strength of its line is a static failure. The safety factor scales the cycle's amplitude
and mean together onto the line; a history's is that of its worst cycle. Gerber's parabola lies
outside Goodman's line, and Goodman's outside Soderberg's: for the same cycle n_Gerber > n_Goodman >
n_Soderberg.

### Life and damage

Basquin's relation gives the cycles to failure at σ_ar:

```text
N = ½ (σ_ar / σf′)^(1/b)
```

A tabulated curve is the same relation piece by piece: between two points the life follows the
straight line through them on log-log axes, and beyond the end points the end segments are
extended (the summary counts the nodes read there). At and below the endurance limit the life is
infinite. The life never falls below half a cycle, one reversal: an equivalent amplitude at or
above the curve's value there (σf′ for Basquin), or a mean past the strength, is a **static
failure**, and so is a **peak stress at or above the ultimate strength**: the largest tensile
stress of the history, stress factor · m · λ_max (or λ_min for a negative m), is checked against
S_u wherever the material has one, whatever the correction and even for a load with no cycle at
all (a constant load, max = min). A static failure has life 0 (the field draws it at half a cycle,
log₁₀ 0.5), damage at least 2 and a safety factor at most 1 (S_u divided by the peak), and is
counted in `results.fatigue` (`static_failures`) with a warning. A material without an ultimate
strength gets no static check: the reply says `static_check: false` and a warning names it (SWT
with no strengths, for example, has none). Compressive peaks are not checked.

The damage of one block is the Palmgren–Miner sum over its cycles (Miner 1945),

```text
D = Σ nᵢ / Nᵢ,
```

and failure is predicted at D = 1. What the life means depends on the loading:

- a constant amplitude loading: the block is one cycle, D = 1/N, and the life is N cycles;
- a history counted `repeated`: the block repeats with every cycle closed, and the life is 1/D
  blocks (how many times the history can be applied);
- a history counted `once`: the ranges left at the end are half cycles, so D is the damage of one
  pass, but 1/D is not a number of repeats (a block applied again would close those half cycles
  into full ones; for [0, 1, −1, 0] the once damage is about half the repeated damage). The
  reply gives the damage and no life (`life_unit` is `none`, the node entries have a null
  `life`), and `fatigue_log_life` is refused; use `fatigue_damage`, or `repeated` for a life.

A node on the boundary of two materials with fatigue data takes the
worse of them.

## Output

| Field | Meaning | Where it has no value |
|---|---|---|
| `fatigue_log_life` | log₁₀ of the life: cycles for a constant amplitude loading, blocks for a history counted `repeated` (not available for `once`); infinite below the endurance limit (in `results.contour` drawn in the top band and left out of the legend's range; in `results.path` the interpolated value is null and the reply's `infinite_life` is true there, with a matching CSV column, null in both outside the part) | nodes of no element whose material has fatigue data |
| `fatigue_damage` | the Miner damage of one cycle or block (0 below the endurance limit) | as above |
| `fatigue_safety_factor` | the safety factor against the endurance limit | as above, and every node when no material in use has an endurance limit (the field is refused) |

The life is shown on a log scale because it spans orders of magnitude: a linear legend would put
every node but the hot spot in one band. The numbers are pure (no unit). `solve.probe` adds a
`fatigue` entry to each node of the step: `life`, `log10_life`, `damage`, `safety_factor`,
`material`, the governing cycle's `stress_amplitude`, `mean_stress` and `equivalent_amplitude`
(stresses in the display units), and the flags `static_failure`, `above_yield` and
`extrapolated` (an infinite life reads `"infinite"`). `results.fatigue` gives the counted
`cycles` (range, mean and count in load factors, at most 200 listed, `distinct_cycles` and
`cycles_per_block`), the `life_unit`, the node of least life with its coordinates and the same
entry, the node of least safety factor, counts of the nodes assessed, of infinite life, of static
failures, of those above the yield strength and of those read beyond a table, and `warnings`.

Nothing is stored with the result: the fields are computed from the step's saved result and the
saved assessment and data whenever they are shown, so they survive saving and opening a file
without any change to its results part, and follow an edit of the loading or the data at once.

## Verification

| Case | Reference | Result | Tolerance | Test |
|---|---|---|---|---|
| Basquin's curve read both ways, σa = σf′ (2N)^b at N = 1 to 10⁹; the cycle form S = A N^b, A = σf′ 2^b | theory (Basquin 1910) | exact to 10⁻¹² | 10⁻¹² | `the_basquin_curve_is_reproduced_exactly` |
| Tabulated curve: exact at the points; at the geometric mean of two amplitudes the geometric mean of their lives; extended end segments | log-log interpolation | exact to 10⁻¹² | 10⁻¹² | `a_table_interpolates_on_log_log_axes` |
| Mean stress corrections, σa = 150, σm = 200, S_u = 800, S_y = 600, S_e = 300: σ_ar Goodman 200, Gerber 160, Soderberg 225; n Goodman 4/3, Soderberg 1.2, Gerber's (n σa, n σm) on its parabola | the formulas (Dowling ch. 9; Shigley), by hand | exact to 10⁻¹³ | 10⁻¹³ | `mean_stress_corrections_follow_their_formulas` |
| SWT, σa = σm = 100 (σ_max = 200): σ_ar = √(200 · 100); σa = 150, σm = 200: √52 500 and n = S_e/√52 500; σ_max ≤ 0: no damage, n infinite; no strength needed; the same bar in Pa and GPa | the formula (Smith, Watson, Topper 1970), by hand | exact to 10⁻¹⁵ (σ_ar, n); 10⁻⁹ in log₁₀ N on the bar | 10⁻⁹ | `mean_stress_corrections_follow_their_formulas`, `a_uniform_bar_reproduces_the_basquin_curve_and_the_mean_stress_corrections`, `fatigue_data_and_replies_follow_the_display_units` |
| The sign of a signed measure near pure shear: shear 200 with σ_xx = 0, ±10⁻⁶, ±0.05 reads tensile with the same damage (Goodman, SWT; both measures); σ_xx = −2 reads compressive | the sign rule | damage within 1 % | 10⁻² | `the_sign_of_a_signed_measure_is_stable_near_pure_shear` |
| Static overload: a mean of 1000 against S_u = 800 with no cycle gives life 0 and n = 0.8; a peak of 1000 past S_u is a failure under every correction; just under S_u is not; the stress factor scales the mean (σ_ar = 200/0.75) | by hand | exact | 10⁻¹⁵ | `a_static_overload_and_a_peak_at_the_strength_are_failures`, `a_static_overload_is_reported_not_called_safe`, `the_stress_factor_scales_the_mean_too` |
| A history counted once and repeated: [0, 1, −1, 0] gives N(200) blocks repeated; once gives damage ½/N(200) + 1/N(100) and no life | by hand | exact | 10⁻⁹ | `a_history_counted_once_has_damage_and_no_life` |
| Editing an assessed step to a procedure that is not linear static is refused; every accepted command leaves a file that opens | — | refusal; round trip | — | `an_assessed_step_cannot_become_a_step_fatigue_cannot_read`, `an_assessment_round_trips_and_belongs_to_a_static_step` |
| A path through nodes of infinite life marks them | — | `infinite_life` true | — | `a_path_marks_where_the_life_is_infinite` |
| Rainflow of ASTM E1049-85's example history −2, 1, −3, 5, −1, 3, −4, 4, −2: ½ cycle of range 3, 1½ of 4, ½ of 6, 1 of 8, ½ of 9, each mean by hand | the example history as widely reproduced from ASTM E1049-85 (the standard itself was not read; the count redone by hand with the three-point procedure) | exact | exact | `rainflow_counts_the_astm_e1049_example` |
| The same history repeated (Downing and Socie's method, counted by hand from 5): full cycles of 3, 4, 7 and 9 | Downing & Socie 1982 | exact | exact | `rainflow_counts_the_astm_e1049_example` |
| Rainflow invariants: sign and scale of the history, plateaus and non-reversals, rotation of a repeated block, half a cycle per range | theory | exact | exact | `rainflow_is_invariant_and_counts_every_reversal` |
| Miner's sum of a two-block history, 3 cycles at ±1 and 5 at ±0.5 of 400 MPa | hand calculation (Miner 1945) | exact to 10⁻¹² | 10⁻¹² | `miner_sums_a_two_block_history` |
| Uniform bar, 10 hex8, σ = 200 MPa, at R = −1 and R = 0 (none, Goodman, Gerber, Soderberg) and a cycle from −0.5 to 1.5: every node's life Basquin's at the formula's σ_ar | theory | equal to 10⁻⁹ in log₁₀ N at all 44 nodes | log₁₀ N to 10⁻⁹ | `a_uniform_bar_reproduces_the_basquin_curve_and_the_mean_stress_corrections` |
| Uniform bar, S_e = 150 MPa: n = 0.75 at σa = 200; at σa = 100 n = 1.5, infinite life, no damage | theory | exact to 10⁻⁹ | 10⁻⁹ | `the_endurance_limit_gives_infinite_life_and_the_safety_factor` |
| Uniform bar, the two-block history as load factors: damage 3/N(200) + 5/N(100) per block | hand calculation | exact to 10⁻⁹ | 10⁻⁹ | `a_two_block_history_gives_miners_sum` |
| Plate with a hole in tension, plane stress quad8 (the mesh of `kirsch-plate-plane-stress`), S = 100 MPa fully reversed: least life at the top of the hole, at Kt S | Kirsch (1898): Kt = 3, N = ½ (300/1000)^(−10) = 84 675 | Kt 3.0029 (+0.10 %); N 83 867 (−0.95 %), log₁₀ N off by 0.0042 | Kt within 2 % (the oracle case's), so log₁₀ N within 10 log₁₀ 1.02 = 0.086 | `the_life_at_a_hole_follows_its_stress_concentration` |
| Saved and opened: the same fields, the result current; undo of an assessment | — | bit for bit | exact | `fatigue_fields_survive_save_and_open` |
| Hostile input: bad curves, strengths, ratios, factors, histories, names and choices; missing strengths named; refusals change nothing | — | located errors | — | `hostile_assessments_are_refused`, `hostile_fatigue_data_is_refused`, `hostile_fatigue_input_is_refused` |

## Limitations

- **Stress-life only.** No strain-life (ε-N, Coffin–Manson, Morrow, and the strain form of
  Smith–Watson–Topper; the `swt` correction here is its stress-life form) and no
  Neuber or Glinka notch plasticity: where the peak stress yields, the elastic stress the method
  reads is wrong. The summary counts the nodes whose von Mises stress at the largest load factor
  exceeds the yield strength and warns that stress-life does not apply there.
- **Proportional loading only.** One step's stress field times load factors. Several load cases
  with their own histories, rotating principal directions and multiaxial critical-plane methods
  (Findley, Dang Van, Fatemi–Socie) are not provided.
- **Linear damage.** Palmgren–Miner sums damage without regard to the order of the cycles; real
  sequences (high then low, low then high) can fail at sums well away from 1.
- **No welds**, no welded-joint detail categories or hot-spot stresses, no crack growth (Paris),
  no fracture mechanics.
- **Cost and cancelling.** The fatigue fields are recomputed from the result on every contour,
  probe, path and summary call (nothing is cached), over all nodes and distinct cycles; the work
  is capped at 5·10⁸ node-cycle evaluations, and a cancel is honoured every 2048 nodes. A long
  history on a large mesh is slow to redraw.
- **Signed measures near pure shear.** Where shear dominates, the sign of a signed measure is
  decided by the sign rule above, not by the physics: see [Stress measures](#stress-measures).
- **Nodal stresses.** The stress at a node is the average over the elements around it, as the
  step stores it; a coarse mesh at a notch under-reads the peak, and so over-reads the life
  (by |1/b| times the stress error on log₁₀ N: ten times for b = −0.1). Refine there. Shells are
  read at their top surface (the other surface is not assessed, with a warning); beams and
  trusses are not assessed.
- **One curve per material**, at one temperature, environment and surface; no statistical
  scatter or survival probability: the curve's own basis (median, design curve) carries over.
- **The data is the user's.** FieldCraft ships no fatigue constants; the results are as good as
  the curve, its endurance limit and the stress factor entered.
- **Display.** The desktop app's contour menu lists the fatigue fields only through the command
  (`results.contour`); a linear static animation scales them with the load as if they were
  stresses.

## Implementation

| File | Contents |
|---|---|
| `crates/model/src/materials.rs` | `Fatigue`, `SnCurve`, `SnPoint`, their validation |
| `crates/model/src/fatigue.rs` | `FatigueAssessment`, `Loading`, `Counting`, `StressMeasure`, `MeanStress`; `Model::set_fatigue` |
| `crates/model/src/steps.rs` | `StepInputs`: materials without their fatigue data, so results stay current |
| `crates/fem/src/fatigue.rs` | `measure`, `reversals`, `rainflow`, `rainflow_repeated`, `life`, `equivalent_amplitude`, `safety_factor`, `assess` |
| `crates/engine/src/fatigue.rs` | `evaluate` over a result, the nodes' materials, the fields |
| `crates/engine/src/cmd/fatigue.rs` | `model.fatigue.set`, `model.fatigue.delete`, `results.fatigue`; the probe entries |
| `crates/engine/src/cmd/define.rs` | `model.material.fatigue` |

## References

- O. H. Basquin. *The exponential law of endurance tests.* Proc. ASTM 10:625–630, 1910.
- M. A. Miner. *Cumulative damage in fatigue.* J. Appl. Mech. 12:A159–A164, 1945.
- W. Gerber. *Bestimmung der zulässigen Spannungen in Eisenconstructionen.* Z. Bayer. Architekten- u. Ingenieur-Vereins 6:101–110, 1874.
- J. Goodman. *Mechanics Applied to Engineering.* Longmans, Green, London, 1899.
- C. R. Soderberg. *Factor of safety and working stress.* Trans. ASME 52 (APM-52-2):13–28, 1930.
- K. N. Smith, P. Watson, T. H. Topper. *A stress-strain function for the fatigue of metals.* J. Materials 5(4):767–778, 1970.
- T. Matsuishi, T. Endo. *Fatigue of metals subjected to varying stress.* Japan Society of Mechanical Engineers, Fukuoka, 1968.
- S. D. Downing, D. F. Socie. *Simple rainflow counting algorithms.* Int. J. Fatigue 4(1):31–40, 1982.
- ASTM E1049-85 (reapproved 2017). *Standard Practices for Cycle Counting in Fatigue Analysis.* ASTM International.
- N. E. Dowling. *Mechanical Behavior of Materials*, 4th ed. Pearson, 2013 (chapter 9).
- R. I. Stephens, A. Fatemi, R. R. Stephens, H. O. Fuchs. *Metal Fatigue in Engineering*, 2nd ed. Wiley, 2001.
- J. A. Bannantine, J. J. Comer, J. L. Handrock. *Fundamentals of Metal Fatigue Analysis.* Prentice Hall, 1990 (chapters 2 and 5).
- R. G. Budynas, J. K. Nisbett. *Shigley's Mechanical Engineering Design*, 10th ed. McGraw-Hill, 2015 (the chapter on fatigue failure resulting from variable loading; section numbers not checked).
- G. Kirsch. *Die Theorie der Elastizität und die Bedürfnisse der Festigkeitslehre.* Z. VDI 42:797–807, 1898.
