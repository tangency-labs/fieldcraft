# Stress linearisation

Pressure-equipment design checks classify the stress through a wall into membrane, bending
and peak parts along a stress classification line, following the stress-integration method of
ASME BPVC Section VIII, Division 2, Annex 5-A. `results.linearize` does this on a solved linear
static step, headlessly, through MCP and in the desktop control channel:

```json
{"command": "results.linearize", "params": {"step": "Pressure", "from": [10, 0, 0.5], "to": [20, 0, 0.5], "points": 41}}
```

`from` and `to` are the two ends of the line, on the two surfaces of the wall, in undeformed
reference coordinates (display length units). The line is sampled at `points` equally spaced
stations (odd, 5 to 1001, default 41): each station takes the containing solid element's own
nodal stresses, extrapolated from its integration points and interpolated with its shape
functions, so a jump at a material or thickness change is kept rather than averaged away. The
integrals are evaluated with Simpson's rule, which is exact for a cubic variation; the result
converges as points are added.

## What the reply holds

Every tensor is in the line's local frame (`frame`) and lists its six components in the order
`components` gives: xx, yy, zz, xy, xz, yz. The local x axis runs along the line (the
through-thickness direction); z is the global axis least aligned with the line, Z preferred,
made orthogonal to it; y completes the right-handed frame. A radial line in a horizontal plane
therefore reports the hoop stress as yy and the axial stress as zz.

- `membrane`: the mean of each component through the wall, (1/t) ∫₀ᵗ σ dx.
- `bending`, at `from` and at `to`: the linear part with zero mean, (6/t²) ∫₀ᵗ σ (t/2 − x) dx,
  positive at `from` and negative at `to`. Following the annex, only the two normal stresses in
  the section (yy, zz) and the shear between them (yz) carry a bending part; the normal stress
  along the line (xx) and the shears on its planes (xy, xz) are membrane and peak only.
- `peak`, at both ends: the total there less the membrane and the bending at that end.
- `total`, at both ends: the sampled stress at the surfaces.
- `equivalent`: the von Mises stresses of the membrane tensor, of membrane plus bending at
  both ends and of the total at both ends.
- `length`, `samples` (each station's distance from `from` and its stress), `elements`.

Allowable stresses are not judged: whether Pm or Pm + Pb meets an allowable depends on the
design code edition and the material, which is the engineer's call.

## Refusals

The command refuses, by name and without changing anything: a stale result; a frequency,
nonlinear static, explicit, heat or fluid step; axisymmetric and other 2D models (linearisation
of axisymmetric sections, with the annex's radius-dependent weighting, is a follow-up);
coincident ends or a line shorter than 10⁻⁹ of the model's size; an end or station outside the
solid elements or across a gap (the first such point is named); an even number of points.

## Verification

The oracle case `tests/oracle/lame-cylinder-linearisation` linearises Lamé's thick-walled
cylinder (a = 10, b = 20, p = 100, plane strain) from the bore to the outside and compares the
hoop membrane (100), the hoop bending at the bore (45.48), the hoop peaks (21.18 and 12.15),
the radial membrane (−33.33), the axial membrane (20) and the membrane von Mises stress
(116.24) with the closed-form integrals of Lamé's solution.
