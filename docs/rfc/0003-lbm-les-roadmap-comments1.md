# Comments on RFC 0003 — LBM-LES verification and research priorities

I support resuming the startup investigation and adding the stream-stability warning (options A
and C). I would reserve the D/80 GPU run until D/40 gives us a statistically stationary wake
and an independent measurement of shedding frequency. A uniform D/40 run can answer those
questions. On the GB10, the driver's default 24 D × 12 D × 3.125 D box at D/40 needs
960 × 480 × 125 cells. At 76 bytes per cell its single population buffer needs 4,377,600,000
bytes, exceeding the reported per-buffer limit of 4,294,967,292 bytes. A uniform **GPU** run
in that box therefore depends on the buffer split in #699. A uniform CPU run remains an option
for measuring the wake. Coarse-level stability also needs a multilevel run that exercises the
speeds, relaxation times and interfaces planned for D/80.
[GPU memory limits](../manual/fluid/gpu.md) · [GPU buffer split](https://github.com/tangency-labs/fieldcraft/issues/699)

The relative contributions of resolution and startup remain open. At 16 cells per diameter, changing only the wake
perturbation moved the circular cylinder's mean drag coefficient from 1.437 to 1.113 and its
mean recirculation length from 0.60 to 2.1 diameters. That is already a large change at fixed
cost and resolution. [RFC 0003, §2.1](0003-lbm-les-roadmap.md#21-bluff-bodies) ·
[Startup reproduction](https://github.com/tangency-labs/fieldcraft/issues/1153)

## Decisions before D/80

**Separate the ramp from the wake perturbation.** In the short fixed-σ reproducer, the ramped
inlet starts at zero while the initial stream-function field approaches the free-stream speed
away from the cylinder. Add a quiescent-start option to the benchmark driver and run the pair
with the same ramp, collision and zero perturbation. Record where and when the first speed or
density excursion appears. Keep both ramp runs unseeded. The 2.5-period ramp takes about
20.8 D/U in this box. At time zero the probe's 1% perturbation would add cross-flow to fluid
still at rest. Its spanwise mean cross-stream component is nonzero, and the added velocity
field has nonzero divergence.

For the seed comparison, reproduce the original pair with the inlet at full speed and the
stream-function start with no ramp. Apply the probe's perturbation to one member of the pair
at initialization through `set_equilibrium`. Measure spanwise velocity fluctuations in both wakes. This
tests the proposed three-dimensional transition under the conditions that produced the drag
difference. Keep the perturbation as a named benchmark input pending the late-time comparison.
[Bluff-body driver](../../crates/lbm-gpu/tests/bluff_body_les.rs)

**Require agreement on the wake as well as the drag.** The square-cylinder run gave a mean drag
within 4% of its reference while the strongest lift peak was at 0.057, against a shedding
Strouhal number of 0.132. For the circular cylinder, agree on experimental bands for drag and
shedding frequency before the run, confirm shedding with a downstream velocity probe, and check
the mean recirculation length against measured wake data. Measure from the rear surface
($x_c+D/2$) to the downstream zero crossing of the mean centreline velocity. Parnaudeau et al.
report $L_r/D=1.51$ from particle-image velocimetry (Table II). The driver's peak picker
reports Fourier bins. With a record of 100 D/U, adjacent bins are about 0.01 apart in Strouhal number.
For a tighter frequency tolerance, use a longer record or an estimator with demonstrated
uncertainty. Parnaudeau et al. provide both hot-wire and particle-image measurements at Re 3900.
[RFC 0003, §2.1](0003-lbm-les-roadmap.md#21-bluff-bodies) ·
[Parnaudeau et al. (2008)](https://doi.org/10.1063/1.2957018)

**Keep the channel study at a fixed physical matching height.** The matched wall hardcodes
four rows, so doubling the resolution from 16 to 32 cells per half-height moves the sample
from 0.25δ to 0.125δ. Make the sample row and the corresponding number of GPU mixing-length
rows configurable, then choose grid/row pairs with the same sampling height in units of δ.
[On-node wall](../../crates/lbm/src/bounce.rs) · [GPU setup](../../crates/lbm-gpu/src/native.rs)

The Musker momentum diagnostic excludes on-node walls. The matched Spalding wall already
accounts for each wall node's half-cell body force: its diagonal links carry the requested
tangential stress less ½ρa in lattice units. The short CPU test checks the total wall force
against the requested shear to a relative tolerance of 10⁻¹², and the GPU has parity tests.
In the long GPU run, compare requested and applied shear in f32 with the same half-cell
accounting. Over the averaging window, compare the force on both walls with the total drive
$a(M+\tfrac{1}{2}\sum_w \rho_w)$ in lattice units, where M is the interior fluid mass and each wall node
uses its sample-cell density. Report the change in fluid momentum over that window. All six
completed runs hit their time cap, so stationarity still needs to be demonstrated before
interpreting bulk velocity and the mean profile against Lee and Moser.
[Musker diagnostic](../verification.md#requested-applied-and-carried-wall-stress) ·
[Matched-wall CPU check](../../crates/lbm/tests/wall_matching.rs) ·
[GPU parity](../../crates/lbm-gpu/tests/verification.rs) ·
[Long channel study](../../crates/lbm-gpu/tests/matched_channel.rs) ·
[Lee & Moser (2015)](https://doi.org/10.1017/jfm.2015.268)

**Check the curved-wall force used in LES.** The study claimed in #940 covers halfway
bounce-back with TRT at τ = 3. For the LES question, start with the existing interpolated-wall
TRT test: an exact sphere at radii of 3, 6 and 9 cells at τ = 1, plus a τ study at radius 3.
The cylinder driver uses interpolated links and HRR on a triangulated surface. Extend the
sphere test at fixed solid fraction with a finer radius and controlled sub-cell offsets.
Compare exact and triangulated surfaces at the same placement, then exercise HRR at relevant
low τ where a creeping-flow reference can be resolved. Compare drag against the analytical
sphere-array reference at each step.
[Halfway study](https://github.com/tangency-labs/fieldcraft/issues/940) ·
[Interpolated-wall test](../../crates/lbm/tests/interpolated.rs) ·
[Cylinder geometry and collision](../../crates/lbm-gpu/tests/bluff_body_les.rs)

**Exercise the coarse levels before D/80.** On the multilevel run, record the collision,
characteristic lattice speed, relaxation time and any growing disturbances on each level.
I support warning when a level enters the instability region found in the sampled
uniform-stream analysis. A hard refusal would need stronger evidence about the bound's
applicability to this flow. Use comparable averaging windows, startup and wall treatment
for the D/40 and D/80 refinement comparison.
[Stability analysis](https://github.com/tangency-labs/fieldcraft/issues/1154)

## Papers and the work they motivate

My first two research priorities are the wake probes for option A and a convected-vortex test
for the first two-level CPU implementation in option D. Both can inform the D/80 decision.
For changes to the solver itself, I would prototype consistent initialization for established
flows and measure the channel stress balance before changing its near-wall model.

1. **Wake sampling — Lehmkuhl et al. (2013).** Their Re 3900 cylinder DNS found a slow
   modulation of the recirculation region at a nondimensional frequency of about 0.0064, or a
   period of roughly 156 D/U. The RFC's 100–200 D/U seed comparison averages over less than
   one such period. Add the two velocity probes from their Fig. 4, one in the shear layer and
   one on the wake centreline. Record the signals alongside lift, base pressure and
   recirculation length. Use stable D/16 runs to compare the two starts over several modulation
   periods, and repeat at D/40 if affordable. The RFC budgets about five GB10-hours for
   **150 D/U total** at refined D/80, including startup. At the same throughput, four periods
   (625 D/U) would cost about 21 GB10-hours for sampling alone. If the shorter D/80 record
   cannot support the agreed uncertainty, budget more time or report it as inconclusive.
   The source of FieldCraft's 0.057 peak remains open. Use Parnaudeau et al.'s experimental
   observations for physical validation.
   [Lehmkuhl et al., *Low-frequency unsteadiness in the vortex formation region of a circular cylinder*](https://doi.org/10.1063/1.4818641) ·
   [Parnaudeau et al.](https://doi.org/10.1063/1.2957018)

2. **Refinement noise — Astoul et al. (2020).** A vortex crossing a change in lattice spacing
   can generate spurious sound and vorticity. Add a convected vortex crossing one 2:1 interface
   to the CPU tests in #1107, which already cites this paper. Run a uniform-grid control and
   measure the excess pressure pulse and change in transmitted vorticity in each crossing
   direction. Couette stress continuity checks steady transport, and the sponge's incident
   acoustic pulse checks absorption. The vortex tests what happens when the wake itself crosses
   the interface. [Astoul et al.](https://arxiv.org/abs/2004.11863) ·
   [FieldCraft refinement design](https://github.com/tangency-labs/fieldcraft/issues/1107)

3. **Initialization — Mei, Luo, Lallemand & d’Humières (2006).** For a prescribed
   divergence-free velocity field, their iteration obtains consistent initial pressure and
   non-equilibrium populations. Start with the cheaper quiescent-fluid ramp control above.
   If FieldCraft also needs to start from established flow, prototype the iteration on the
   paper's two-dimensional decaying Taylor–Green vortex, where pressure and velocity decay
   can be checked analytically. Then compare initializers on the cylinder with the same
   velocity field and compatible inlet. Measure the initial pressure pulse and the location
   and time of the first speed excursion. This isolates the contribution of equilibrium
   initialization in a sheared flow while holding inlet compatibility fixed.
   [Mei et al.](https://doi.org/10.1016/j.compfluid.2005.08.008)

4. **Channel stress — Liu, Xu & Huang (2025).** At fixed matching height, plot the mean
   viscous, resolved Reynolds and modeled shear across the channel. In the incompressible,
   statistically stationary, fully developed plane channel, their sum follows
   $\tau_w(1-y/\delta)$ from the wall to the centreline, where δ is the half-height.
   Inspect rows 3–5 in particular. At Re_τ 1000 and 16 cells per δ, the profile errors
   relative to DNS are −5.65%, +0.05% and +5.56% respectively. Row 4 is the sample and the
   first row outside the mixing-length layer. That layer uses κ = 0.384 with a Spalding
   slope of 0.41, giving a log-layer eddy-stress fraction of about 0.877 in the wall-law
   derivation. Adapt the GPU test's row-wise `Budget` accumulator to use the actual layer
   and bulk relaxation rates and non-equilibrium momentum flux. Its current gradient-based
   Smagorinsky estimate would miss the change of closure. If wall momentum transfer closes
   but the stress partition and mean velocity profile are wrong, prototype Liu et al.'s
   near-wall eddy-viscosity correction. Their finite-difference formulation requires a
   discrete LBM derivation.
   [Measured profiles](../verification.md#matched-on-node-wall-on-the-gpu-against-lee--moser-at-re_τ-1000-2000-and-5200-1108) ·
   [Layer rates](../../crates/lbm/src/bounce.rs) · [Budget helper](../../crates/lbm-gpu/tests/verification.rs) ·
   [Liu, Xu & Huang](https://doi.org/10.1016/j.jcp.2025.114029) ·
   [Lee & Moser](https://doi.org/10.1017/jfm.2015.268)

5. **Collision near τ = ½ — Wissocq, Coreixas & Boussuge (2020).** Their two-dimensional
   analysis shows that an under-resolved physical wave can experience excessive,
   direction-dependent damping under regularization. Extend the D3Q19 analysis to measure
   shear-wave decay as well as growth for axial and oblique waves at the coarse levels'
   actual speeds and relaxation times. Excess damping at wake-relevant scales or growth of
   an unstable mode would give the proposed D3Q27 collision a specific requirement to meet.
   [Wissocq et al.](https://arxiv.org/abs/2006.07353) ·
   [FieldCraft stability study](https://github.com/tangency-labs/fieldcraft/issues/1154)

6. **Design of the D/40 and D/80 runs — Suss & Mary (2026).** This LBM study of the Re 3900
   cylinder considers mesh resolution, wall treatment and statistical sampling together.
   I would use its study design when specifying the two runs and set the acceptance bands
   from experimental measurements.
   [Suss & Mary, AIAA 2026-3485](https://doi.org/10.2514/6.2026-3485)

## Suggested issues

Each draft can be filed as a follow-up on its linked tracking issue. Start with the ramp
comparison and wake probes. The interface test depends on the first two-level CPU stepper.
The established-flow initializer is conditional on what the ramp investigation finds.

- [Compare quiescent and stream-function starts in the ramped cylinder test](0003-issue-drafts/01-cylinder-startup.md)
- [Measure cylinder shedding independently of the lift spectrum](0003-issue-drafts/02-cylinder-wake-probes.md)
- [Close the force budget over the matched-wall GPU averaging window](0003-issue-drafts/03-matched-wall-force.md)
- [Specify the matched Spalding wall’s sampling distance in grid rows](0003-issue-drafts/04-matching-height.md)
- [Measure the shear-stress partition across the matched-wall sampling row](0003-issue-drafts/05-channel-stress-profile.md)
- [Measure the pressure pulse from a vortex crossing one CPU refinement interface](0003-issue-drafts/06-refinement-vortex.md)
- [Isolate triangulation error in interpolated-wall Stokes drag](0003-issue-drafts/07-interpolated-sphere.md)
- [Measure D3Q19 HRR shear-wave damping on a moving stream](0003-issue-drafts/08-oblique-shear.md)
- [Test consistent population initialization on a decaying Taylor–Green vortex](0003-issue-drafts/09-consistent-initialization.md)
- [Specify the D/40 evidence and D/80 budget for cylinder acceptance](0003-issue-drafts/10-cylinder-refinement-protocol.md)
