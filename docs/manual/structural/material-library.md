# Cited material library

The Material menu offers **Browse Material Library** and **Add Library Material**. The same
actions are `model.material.library.list` and `model.material.library.add` in scripts, the CLI and
MCP. Browse all ten entries, or filter `family` to `steel`, `stainless`, `aluminium`,
`titanium` or `fluid`. Adding an entry creates an ordinary editable material; it refuses a taken name.

```json
{"command":"model.material.library.list","params":{"family":"stainless"}}
{"command":"model.material.library.add","params":{"entry":"EN10088-14301","name":"Stainless304","reference_temperature":"20 degC"}}
{"command":"model.material.library.add","params":{"entry":"air","name":"Cooling_air","temperature":"40 degC"}}
```

The report gives every supplied value in the current display units, its basis, source
document/edition/table or clause and, where checked against the document, page, with a publisher or catalogue link where available,
and the product/temperature conditions. `material_si`
is explicitly the material as stored, in SI. `values` are the elastic, density and thermal
properties used; `not_applied` contains strength information. `unavailable` names properties
not supplied by this small catalogue. No network request or external data file is needed.

| Entry id | Material and main qualification |
|---|---|
| `EN10025-S235` | S235 hot-rolled steel, nominal thickness up to 40 mm |
| `EN10025-S355` | S355 hot-rolled steel, nominal thickness up to 40 mm |
| `EN10088-14301` | 1.4301 / 304 stainless, hot-rolled plate up to 75 mm |
| `EN10088-14401` | 1.4401 / 316 stainless, hot-rolled plate up to 75 mm |
| `AA6061-T6` | Wrought 6061-T6; no product-specific strength selected |
| `AA6082-T6` | Extruded 6082-T6; strength row for 5-20 mm specimen cross-section |
| `AA7075-T6` | Wrought 7075-T6; no product-specific strength selected |
| `AA2024-T3` | Bare 2024-T3 sheet, 0.010-0.128 inch, explicitly legacy handbook edition |
| `Ti6Al4V-annealed` | Annealed Ti-6Al-4V; strength row for ASTM B265 sheet/plate, 0.025-1 inch |
| `air` | Dry air, 200-400 K: viscosity, conductivity and sea-level density curves, perfect-gas specific heat |

## Temperature-dependent properties

Some entries give a property as a **curve** of (temperature, value) points instead of a constant;
the report lists them under `curves`, each point in the display units, with its source like a
constant. A curve is linear between its points and held at its end values outside them.

| Entry | Curves | Source |
|---|---|---|
| `EN10025-S235`, `EN10025-S355` | conductivity and specific heat 20-1200 C, through the 5000 J/(kg K) specific heat peak at 735 C | The EN 1993-1-2:2005 (3.4.1) carbon-steel design equations as an open paper reproduces them (Abeysiriwardena and Mahendran, Thin-Walled Structures 180 (2022) 109792, accepted manuscript on QUT ePrints, p. 19, Eqs. (27)-(32)), evaluated densely enough that specific-heat interpolation stays within 0.1 % of its equations |
| `air` | viscosity, conductivity, density, 200-400 K every 10 K | U.S. Standard Atmosphere 1976, eq. (51) and (53) p. 19, eq. (42) at 101325 Pa p. 15 |

Adding such an entry makes a material that keeps the whole curves (saved in the model, so a later
library update never changes its results) and takes each curve's value at `temperature`
(default 20 C) as its constant. The heat transfer solver follows the conductivity, specific
heat and density curves: each element takes its properties at its mean temperature, iterated to
convergence (see [heat transfer](heat-transfer.md)), with a warning where the temperatures leave
a curve. The other solvers and the deck writers use the constants. A `temperature` outside a
curve takes its end value and the reply warns.

Every curve cites a document a reader can open and check. Aluminium has no curves yet: the
EN 1999-1-2 design equations are not reproduced in any open source we could check. Liquid water is
not in the library yet: the published tables we found are NIST Standard Reference Data, which may not be
redistributed without NIST's permission.

The air density is that of the perfect gas at the sea-level pressure; at another pressure, edit
the density. The Standard notes its viscosity equation fails at very high and very low
temperatures, so the curves stop at 200 and 400 K.

The material name defaults to the stable entry id. The stress-free `reference_temperature`
defaults to **20 degrees Celsius** in every display system. This is a modelling choice, not a
measured handbook property. Set it to the assembly's actual stress-free temperature when
calculating thermal strain. The report says whether expansion exists to use that reference.

All five [display presets](../../units.md), including imperial and Fahrenheit overrides, use the
existing public `Quantity` conversion boundary. Stored properties remain SI; changing display
units never changes them. The legacy `unspecified` alias means the mm-tonne-second preset,
and `default` follows the user's setting. Thermal properties support the selected temperature
scale, including Fahrenheit.

## Scope and source choices

These are constant isotropic approximations, with each property's own basis and conditions.
`design` identifies a source's nominal calculation constant, `minimum` a published minimum for
its specified product/condition, and `typical` a representative physical value. They are not
interchangeable. Yield and tensile strengths are informational only: adding an entry never
enables plasticity, damage, erosion, a safety factor or a failure criterion. A supplier's
certificate and the applicable design standard remain necessary for an engineering decision.

The original small curation in [materials.json](../../../crates/model/data/materials.json) selects
individual cited facts; it does not copy a handbook table layout or another product's material
library. Standards are cited by number, edition, clause or table and, where checked, page; the catalogue does
not link to third-party copies. Steel/stainless calculation constants come from the explicitly
cited EN 1993 editions;
6082 constants use EN 1999-1-1 and the product-specific Hydro sheet. Stainless physical values
are from Outokumpu's public Core/Supra Table 7. NASA's public 1972 alloy handbooks identify the
6061/7075 elastic constants and physical-property temperature ranges. Ti-6Al-4V combines NASA's
nominal design tensile modulus with TIMET's cited room-temperature physical values: it is an
isotropic approximation, not a single measured lot or texture-dependent stiffness tensor.

The 2024-T3 entry explicitly cites **MIL-HDBK-5H (1998), Table 3.2.3.0(b1)**, in the section
reproduced with **Change Notice 1 (2001)**. Its link is the [official DLA
catalogue](https://quicksearch.dla.mil/qsDocDetails.aspx?ident_number=53876), which identifies the
superseding 5J (2003); these values are deliberately not relabelled as that edition. Older
handbook constants do not establish current procurement requirements or design allowables.

Missing properties stay absent. For example, 2024 has elasticity and density here, 6061/7075
have no conductivity or specific heat, and 6082 has conductivity but no specific heat. Supply
independently justified missing properties with `model.material.edit` before an analysis that
needs them. Mean expansion coefficients over a stated interval are applied as constants.

## Reproduce a conduction check

Run the original [conduction recipe](../../recipes/material-library-conduction.json):

```sh
fieldcraft-cli run docs/recipes/material-library-conduction.json
```

It adds library 304 stainless and uses its cited constant conductivity, 15 W/(m K) at 20 C,
in a 20 by 30 by 120 mm insulated bar with its ends held at 20 and 30 C. With 36 hex8 cells,
the exact constant-property solution is midpoint 298.15 K, axial heat flux -1250 W/m2 and
0.75 W through either end, with opposite signed reactions. The recipe probes those results.
This is a synthetic implementation check of the material/units/solver path, not experimental
material validation, a transient prediction or a coupled heat exchanger. Real thermal design
must account for property variation, contacts, convection and radiation as appropriate.
