# Path plots

After solving a solid step, open **Results > Path Plot**. Add reference-coordinate points by
typing their coordinates or enabling **Pick surface points** and clicking the visible mesh.
The editor keeps SI coordinates internally: changing display units does not round a pick.
Picking a deformed or curved face maps the actual displayed triangle back to the element's
reference geometry. It does not snap to a nearby node.

Choose a nodal field and the number of equally spaced arc-length stations, then **Sample
path**. Sampling runs in a cancellable background job. The chart retains the sampled step,
field and frame or mode independently of the current contour. Hovering it marks the station
in the viewport. Its overlay follows the viewport's reference or deformed configuration;
coordinates and distances in the query remain reference coordinates. **Export CSV** samples
the draft through the same command and writes the result. Cancellation or an error preserves
the preceding successful plot.

## Commands and units

`results.path` works headlessly, through MCP and in the desktop control channel. Programmatic
calls never open a picker. For example, on the built-in solved cantilever:

```json
{"step":"Load","field":"stress_xx","points":[["20 mm","5 mm","20 mm"],["80 mm","5 mm","20 mm"]],"samples":101,"frame":null,"mode":null,"path":"top-fibre.csv"}
```

The points define a polyline with 2–1000 vertices; repeated vertices are allowed, but its
total length must be positive. `samples` is 2–10000, default 101. Bare point coordinates use
the session's display length unit. Explicit unit strings avoid dependence on that setting.
The reply carries `distance`, `xyz`, `values`, containing `elements`, `units`, `field_title`,
`frame`, `mode` and `configuration: "reference"`. Values and coordinates are converted once
at the display boundary. Modal displacement fields retain the existing dimensionless
normalisation rather than being labelled as a physical length.

The optional CSV has exactly `distance,x,y,z,value` as its header and uses the same numeric
values as the reply, CRLF line endings and empty values at gaps. Units are supplied by the
reply and the chart, rather than extra CSV columns. Writes respect the session's file roots.

## Interpolation and limits

Sampling uses the containing solid's shape functions: tet4, tet10, wedge6, pyramid5, hex8 and
hex20, including curved quadratic geometry. The lowest element label wins on shared faces.
This interpolates the existing nodal field; it does not recompute stresses at integration
points or make a discontinuous stress field continuous. Numerical contact fields are
available where their contributing nodal values are defined. Categorical `contact_status`
is refused because fractional status codes would be misleading.

Proven outside points and undefined nodal contributions give `null` and visible chart gaps;
no nearest-node extrapolation or gap bridging is used. An inverse search that cannot establish
containment or exclusion within its work limit returns an error, rather than calling the
point outside. A straight segment between two picks on a concave surface can leave the solid.
Positive local Jacobian checks do not certify that an entire element is globally injective.

Stale results are refused. Omitted/null `frame` selects the final state; stored frames start
at zero. Nonlinear stored frames contain displacements only: stress, contact and temperature
require the final state until complete intermediate states are stored. Heat flux is final-only;
stored heat frames provide temperatures. Frequency results use `mode` from one (default one)
and reject frames. Other procedures reject non-null modes. Fluid cell fields and non-solid
element families are outside this command's scope.

The index is cached for an exact result and model definition. Construction is bounded by
500000 mesh elements, 2000000 result nodes, 96 MiB index storage and 14000000 basis/enclosure
work units. A query separately permits 4000000 basis/enclosure work units, 2000000 BVH visits
and 4096 candidate elements per station. Index construction and sampling check cancellation
and report progress. Higher-order element bounds use conservative Bernstein enclosures,
including outward rounding; a nodal bounding box alone is insufficient for curved elements.

## Reproducible temperature check

Run the original [conductor recipe](recipes/path-temperature.json):

```text
fieldcraft-cli run docs/recipes/path-temperature.json
```

It builds a 20 × 30 × 120 mm insulated rectangular conductor, with an explicitly assumed
constant conductivity of 15 W/(m·K), held at 20 °C and 30 °C at its ends, then solves steady
conduction on the CPU. The 48 hex8 elements reproduce the linear solution:

```text
T(z) = 20 + z / 12 °C  (z in mm)
q_z = -1250 W/m²
```

The centreline path contains 121 stations, so the endpoint temperatures are 20 and 30 °C and
the midpoint is 25 °C. This checks interpolation and unit handling against an analytic field;
the assumed constant is not a material grade or a validated heat-exchanger model. To inspect
it in the app, run the recipe through `engine.script`, show its temperature contour, then
open Path Plot. A path extended beyond either end should show gaps, not extrapolated values.

The original [curved-picking recipe](recipes/path-curved-picking.json) defines one hex20 with
the reference map `x = 10 r`, `y = 10 (s + 0.3 r - 0.4 (1-r²))`, `z = 10 t` in millimetres.
Its positive map Jacobian is constant. It clamps the negative-x face and applies 100 N in -z
at the other face, using explicitly assumed elastic constants. After solving, show its
deformed shape, pick surface points and sample them. The regression checks the genuine
display-triangle picks under all five presets; this is a picking/interpolation check, not a
converged stress benchmark. A two-station query checks both picked endpoints without
assuming their straight connecting chord stays inside the concave solid.
