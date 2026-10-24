# Cited material library

The Material menu offers **Browse Material Library** and **Add Library Material**. The same
actions are `model.material.library.list` and `model.material.library` in scripts, the CLI and
MCP. Browse all nine entries, or filter `family` to `steel`, `stainless`, `aluminium` or
`titanium`. Adding an entry creates an ordinary editable material; it refuses a taken name.

```json
{"command":"model.material.library.list","params":{"family":"stainless"}}
{"command":"model.material.library","params":{"entry":"EN10088-14301","name":"Stainless304","reference_temperature":"20 degC"}}
```

The report gives every supplied value in the current display units, its basis, source
document/edition/table or clause/page/link, and the product/temperature conditions. `material_si`
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

The material name defaults to the stable entry id. The stress-free `reference_temperature`
defaults to **20 degrees Celsius** in every display system. This is a modelling choice, not a
measured handbook property. Set it to the assembly's actual stress-free temperature when
calculating thermal strain. The report says whether expansion exists to use that reference.

All five [display presets](units.md), including imperial and Fahrenheit overrides, use the
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

The original small curation in [materials.json](../crates/model/data/materials.json) selects
individual cited facts; it does not copy a handbook table layout or another product's material
library. Steel/stainless calculation constants come from the explicitly cited EN 1993 editions;
6082 constants use EN 1999-1-1 and the product-specific Hydro sheet. Stainless physical values
are from Outokumpu's public Core/Supra Table 7. NASA's public 1972 alloy handbooks identify the
6061/7075 elastic constants and physical-property temperature ranges. Ti-6Al-4V combines NASA's
nominal design tensile modulus with TIMET's cited room-temperature physical values: it is an
isotropic approximation, not a single measured lot or texture-dependent stiffness tensor.

The 2024-T3 entry explicitly cites **MIL-HDBK-5H (1998), Table 3.2.3.0(b1)**, in the section
reproduced with **Change Notice 1 (2001)**. Its link is an accessible reproduction of the
original handbook text, rather than a commercial material database. The [official DLA
catalogue](https://quicksearch.dla.mil/qsDocDetails.aspx?ident_number=53876) identifies the
superseding 5J (2003); these values are deliberately not relabelled as that edition. Older
handbook constants do not establish current procurement requirements or design allowables.

Missing properties stay absent. For example, 2024 has elasticity and density here, 6061/7075
have no conductivity or specific heat, and 6082 has conductivity but no specific heat. Supply
independently justified missing properties with `model.material.edit` before an analysis that
needs them. Mean expansion coefficients over a stated interval are applied as constants;
there is no temperature-dependent material curve or interpolation in this first library.

## Reproduce a conduction check

Run the original [conduction recipe](materials/conduction.json):

```sh
fieldcraft-cli run docs/materials/conduction.json
```

It adds library 304 stainless and uses its cited constant conductivity, 15 W/(m K) at 20 C,
in a 20 by 30 by 120 mm insulated bar with its ends held at 20 and 30 C. With 36 hex8 cells,
the exact constant-property solution is midpoint 298.15 K, axial heat flux -1250 W/m2 and
0.75 W through either end, with opposite signed reactions. The recipe probes those results.
This is a synthetic implementation check of the material/units/solver path, not experimental
material validation, a transient prediction or a coupled heat exchanger. Real thermal design
must account for property variation, contacts, convection and radiation as appropriate.
