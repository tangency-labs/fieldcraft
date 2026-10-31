//! The commands' parameter schemas ([`crate::schema`]), area by area. Each describes the JSON
//! the command's `run` reads, as its `params` doc string does in words.

use fieldcraft_model::units::Quantity;

use crate::schema::Collection as C;
use crate::schema::Dim as D;
use crate::schema::*;

// ---- Shared pieces ----

/// The smallest fan relaxation a form takes (the model takes any in (0, 1]).
const MIN_FAN_RELAXATION: f64 = 1e-6; // unit-factors: a fraction per time step
/// The largest merge distance a form takes: a cap on hostile input, in the display unit it is typed in.
const MAX_MERGE: f64 = 1.0e9; // unit-factors: a cap on a typed number, not a conversion

/// `{"set": …}` or `{"reference_point": …}`.
const TARGET_OPTIONS: &[Field] = &[
    reference("set", "Set", C::Set).help("A named set of nodes, elements or faces."),
    reference("reference_point", "Reference point", C::ReferencePoint).help("A reference point: a single node that stands for a rigid connection."),
];
const fn target(key: &'static str, label: &'static str) -> Field {
    either(key, label, TARGET_OPTIONS)
}

/// The name of an existing definition, for delete commands.
macro_rules! existing {
    ($label:literal, $c:expr) => {
        &[reference("name", $label, $c).req().help(concat!("The ", $label, " to act on, by name."))]
    };
}

pub const SET_DELETE: &[Field] = existing!("Set", C::Set);
pub const MATERIAL_DELETE: &[Field] = existing!("Material", C::Material);
pub const MATERIAL_LIBRARY_LIST: &[Field] =
    &[choice("family", "Family", fieldcraft_model::library::FAMILIES).help("Optional family filter; empty lists all nine cited entries.")];
pub const MATERIAL_LIBRARY: &[Field] = &[
    choice("entry", "Library entry", fieldcraft_model::library::CHOICES).req().help("Cited constant isotropic properties, with product/temper/temperature limits in the report. Missing data stays absent; strengths never enable plasticity or damage."),
    text("name", "Material name").hint("entry id").help("Optional unique name; empty uses the stable library entry id."),
    num("reference_temperature", "Stress-free temperature", D::Temperature).hint("20 degC").help("Your modelling choice: temperature of zero thermal strain. Defaults to 20 degrees Celsius, regardless of display units."),
];
pub const SECTION_DELETE: &[Field] = existing!("Section", C::Section);
pub const REFERENCE_POINT_DELETE: &[Field] = existing!("Reference point", C::ReferencePoint);
pub const COUPLING_DELETE: &[Field] = existing!("Coupling", C::Coupling);
pub const MASS_DELETE: &[Field] = existing!("Point mass", C::Mass);
pub const SPRING_DELETE: &[Field] = existing!("Spring", C::Spring);
pub const BOUNDARY_CONDITION_DELETE: &[Field] = existing!("Boundary condition", C::BoundaryCondition);
pub const LOAD_DELETE: &[Field] = existing!("Load", C::Load);
pub const FLUID_DOMAIN_DELETE: &[Field] = existing!("Fluid domain", C::FluidDomain);
pub const CONTACT_DELETE: &[Field] = existing!("Contact", C::Contact);
pub const TIE_DELETE: &[Field] = existing!("Tie", C::Tie);
pub const STEP_DELETE: &[Field] = existing!("Step", C::Step);
pub const SET_RESOLVE: &[Field] = existing!("Set", C::Set);

// ---- Model definitions ----

pub const SET: &[Field] = &[
    NAME,
    rows(
        "geometry",
        "Geometry",
        &[
            reference("body", "Body", C::Body).help("A body of the open design.").req(),
            text("face", "Face").help("A persistent face name of the body (e.g. F3); the whole body when empty.").hint("persistent face name"),
            text("edge", "Edge").help("A persistent edge name of the body."),
            text("vertex", "Vertex").help("A persistent vertex name of the body."),
        ],
    )
    .help("Bodies of the design, or their named faces, edges or vertices. Kept by name, so the set survives remeshing.")
    .hint("whole bodies, or their faces, edges or vertices"),
    labels("nodes", "Nodes", LabelOf::Node).help("Mesh node labels; Use selection takes the picked nodes.").hint("labels: 1, 2, 10-20"),
    labels("elements", "Elements", LabelOf::Element)
        .help("Mesh element labels; Use selection takes the picked elements.")
        .hint("labels: 1, 2, 10-20"),
    sides("faces", "Element faces", Side::Face)
        .help(
            "Element faces as element:face (local face number), ranges of elements sharing one: 1-20:3, 25:1. Use selection takes the picked faces.",
        )
        .hint("element:face, e.g. 1-20:3, 25:1"),
    sides("edges", "Element edges", Side::Edge)
        .help("Element edges as element:edge (local edge number), ranges of elements sharing one: 1-20:3, 25:1.")
        .hint("element:edge, e.g. 1-20:3"),
];

const ELASTIC: Field = variants(
    "elastic",
    "Elasticity",
    &[
        variant(
            "isotropic",
            "Isotropic",
            &[
                num("young", "Young's modulus", D::Stress)
                    .help("Stiffness in tension: stress over elastic strain. Steel ≈ 210 GPa, aluminium ≈ 70 GPa, rubber a few MPa.")
                    .req(),
                range("poisson", "Poisson's ratio", D::None, -1.0, 0.5)
                    .help("Lateral contraction over axial stretch, between -1 and 0.5. Metals ≈ 0.3, rubber close to 0.5.")
                    .req(),
            ],
        ),
        variant(
            "engineering_constants",
            "Orthotropic",
            &[
                vector("young", "E1, E2, E3", D::Stress, &["E1", "E2", "E3"]).help("Young's moduli along the three material axes.").req(),
                vector("poisson", "ν12, ν13, ν23", D::None, &["ν12", "ν13", "ν23"])
                    .help("Poisson's ratios between the material axes (contraction along the second axis for a stretch along the first).")
                    .req(),
                vector("shear", "G12, G13, G23", D::Stress, &["G12", "G13", "G23"]).help("Shear moduli in the three material planes.").req(),
            ],
        ),
    ],
)
.help("How the material deforms elastically: the same in every direction (isotropic) or per material axis (orthotropic).")
.def(r#"{"type": "isotropic"}"#);

const HARDENING: &[(&str, &str)] = &[("isotropic", "Isotropic"), ("kinematic", "Kinematic")];

const PLASTIC: Field = variants(
    "plastic",
    "Plasticity",
    &[
        variant(
            "bilinear",
            "Bilinear",
            &[
                positive("yield_stress", "Yield stress", D::Stress).help("The initial yield stress: where elastic behaviour ends and plastic flow starts. 6061-T6 aluminium ≈ 276 MPa, mild steel ≈ 250 MPa.").req(),
                num("hardening_modulus", "Hardening modulus", D::Stress).help("Slope of stress against plastic strain after yield. 0 is perfectly plastic; metals are often 1–5 % of Young's modulus (e.g. 2000 MPa for steel).").def("0"),
                choice("hardening", "Hardening", HARDENING).help("Isotropic: the yield surface grows equally in every direction. Kinematic: it moves with the loading, so reversed loading yields earlier (Bauschinger effect; for cyclic loads).").def("\"isotropic\""),
            ],
        ),
        variant(
            "table",
            "Table",
            &[rows(
                "points",
                "Hardening curve",
                &[positive("strain", "Plastic strain", D::None).help("Equivalent plastic strain of the point; the first point is 0.").req(), positive("stress", "Stress", D::Stress).help("Yield stress at that plastic strain.").req()],
            ).help("Points of yield stress against plastic strain, plastic strain rising from 0 at the initial yield.")
            .req()],
        ),
        variant(
            "johnson_cook",
            "Johnson–Cook",
            &[
                positive("a", "A", D::Stress).help("Johnson–Cook initial yield stress at the reference strain rate and room temperature. Steel ≈ 300–800 MPa.").req(),
                num("b", "B", D::Stress).help("Johnson–Cook hardening coefficient: the stress added by plastic strain.").req(),
                num("n", "n", D::None).help("Johnson–Cook hardening exponent, typically 0.1–0.6.").req(),
                num("c", "C", D::None).help("Johnson–Cook strain-rate coefficient, typically 0.01–0.03; 0 ignores the strain rate.").def("0"),
                num("m", "m", D::None).help("Johnson–Cook thermal softening exponent, typically about 1.").def("1"),
                positive("reference_strain_rate", "Reference strain rate", D::StrainRate).help("The plastic strain rate at which A, B and n were measured, per time unit (often 1/s).").hint("1/s"),
                num("room_temperature", "Room temperature", D::Temperature).help("The temperature at which A, B and n hold; no thermal softening below it.").req(),
                num("melt_temperature", "Melt temperature", D::Temperature).help("The temperature at which strength drops to zero. Steel ≈ 1800 K, aluminium ≈ 900 K.").req(),
                positive("specific_heat", "Specific heat", D::SpecificHeat).help("Heat that raises a unit mass by one degree, for heating by plastic work. Steel ≈ 460 J/(kg·K).").hint("isothermal when empty"),
                range("taylor_quinney", "Taylor–Quinney", D::None, 0.0, 1.0).help("Fraction of plastic work turned into heat, 0 to 1; usually 0.9.").def("0.9"),
            ],
        ),
    ],
).help("Permanent deformation once the stress passes the yield stress: bilinear, a table of points, or Johnson–Cook.");

const DAMAGE: Field = variants(
    "damage",
    "Damage",
    &[variant(
        "johnson_cook",
        "Johnson–Cook",
        &[
            num("d1", "D1", D::None).help("Johnson–Cook failure strain constant. Typical metals: −0.8 to 0.2.").req(),
            num("d2", "D2", D::None).help("Johnson–Cook failure strain constant multiplying the stress-triaxiality exponential.").req(),
            num("d3", "D3", D::None).help("Johnson–Cook stress-triaxiality exponent constant, usually negative.").req(),
            num("d4", "D4", D::None).help("Johnson–Cook strain-rate constant; 0 when empty."),
            num("d5", "D5", D::None).help("Johnson–Cook temperature constant; 0 when empty."),
        ],
    )],
)
.help("When the material fails: Johnson–Cook failure strain built from D1–D5.");
const EROSION: Field = group(
    "erosion",
    "Erosion",
    &[positive("plastic_strain", "Plastic strain", D::None)
        .help("Equivalent plastic strain at which an element is removed, e.g. 0.3–1 for ductile metals.")
        .req()],
)
.help("Removes an element from an explicit analysis once it reaches the plastic strain given.");

const HYPERELASTIC: Field = variants(
    "hyperelastic",
    "Hyperelasticity",
    &[
        variant("neo_hookean", "Neo-Hookean", &[]),
        variant(
            "mooney_rivlin",
            "Mooney–Rivlin",
            &[num("c01", "C01", D::Stress).help("Mooney–Rivlin second constant; C10 comes from the isotropic shear modulus.").req()],
        ),
        variant(
            "yeoh",
            "Yeoh",
            &[
                num("c20", "C20", D::Stress).help("Yeoh second-order constant.").req(),
                num("c30", "C30", D::Stress).help("Yeoh third-order constant.").req(),
            ],
        ),
    ],
)
.help("A large-strain rubber model (nonlinear steps), built on the isotropic constants.")
.hint("large strain in nonlinear steps, from the isotropic constants");

/// The material axes 1, 2, 3 (the global X, Y, Z until material orientations land).
const AXES_123: &[&str] = &["1", "2", "3"];

const THERMAL: Field = group(
    "thermal",
    "Thermal",
    &[
        per_axis("conductivity", "Conductivity", D::Conductivity, Some(0.0), AXES_123)
            .help("Thermal conductivity, one value or one per material axis. Steel ≈ 50 W/(m·K)."),
        positive("specific_heat", "Specific heat", D::SpecificHeat).help("Heat that raises a unit mass by one degree. Steel ≈ 460 J/(kg·K)."),
        per_axis("expansion", "Expansion", D::Expansion, None, AXES_123)
            .help("Thermal expansion coefficient, one value or one per material axis. Steel ≈ 1.2e-5, aluminium ≈ 2.3e-5."),
        num("reference_temperature", "Reference temperature", D::Temperature)
            .help("The temperature at which there is no thermal strain.")
            .hint("0 °C"),
    ],
)
.help("Heat conduction, specific heat and thermal expansion, for heat-transfer and thermal-stress steps.");

/// One domain's coefficients of the Tait equation (melt or solid).
const TAIT_DOMAIN: &[Field] = &[
    positive("b1", "b1", D::Of(Quantity::SpecificVolume))
        .help("Specific volume at the transition temperature b5 and zero pressure. Polymers ≈ 0.8e-3–1.3e-3 m³/kg.")
        .req(),
    positive("b2", "b2", D::Of(Quantity::SpecificVolumePerTemperature)).help("Rise of the zero-pressure specific volume per degree above b5.").req(),
    positive("b3", "b3", D::Stress).help("Pressure sensitivity B of the Tait equation at b5; typically 50–300 MPa.").req(),
    positive("b4", "b4", D::Expansion).help("Fall of the pressure sensitivity with temperature: B = b3 exp(−b4 (T − b5)).").req(),
];

const PVT: Field = variants(
    "pvt",
    "PVT",
    &[variant(
        "tait",
        "Two-domain Tait",
        &[
            group("melt", "Melt", TAIT_DOMAIN).help("Coefficients above the transition temperature.").req(),
            group("solid", "Solid", TAIT_DOMAIN).help("Coefficients at and below the transition temperature.").req(),
            num("b5", "b5", D::Temperature).help("Transition temperature at zero pressure: the glass transition or crystallisation temperature (above absolute zero).").req(),
            positive("b6", "b6", D::Of(Quantity::TemperaturePerPressure)).help("Rise of the transition temperature with pressure, T_t = b5 + b6 p.").req(),
            positive("b7", "b7", D::Of(Quantity::SpecificVolume)).help("Extra solid volume at the transition of a semi-crystalline polymer; 0 for an amorphous one.").def("0"),
            positive("b8", "b8", D::Expansion).help("Temperature rate of the transition volume term b7 exp(b8 (T − b5) − b9 p).").def("0"),
            positive("b9", "b9", D::Of(Quantity::InversePressure)).help("Pressure rate of the transition volume term.").def("0"),
        ],
    )],
)
.help("Specific volume against temperature and pressure (the modified two-domain Tait equation of a polymer): its density ρ = 1/v where the state is known.");

const PHASE_CHANGE: Field = group(
    "phase_change",
    "Phase change",
    &[
        positive("latent_heat", "Latent heat", D::Of(Quantity::SpecificEnergy))
            .help("Heat a unit mass releases as it freezes. Water ≈ 334 kJ/kg, aluminium ≈ 397 kJ/kg.")
            .req(),
        num("solidus", "Solidus", D::Temperature).help("Below it the material is solid and a fluid stops flowing.").req(),
        num("liquidus", "Liquidus", D::Temperature).help("Above it the material is liquid; equal to the solidus for a pure substance.").req(),
    ],
)
.help("Melting and freezing: the latent heat released between the liquidus and the solidus. Needs a specific heat.");

pub const MATERIAL: &[Field] = &[
    NAME,
    ELASTIC,
    positive("density", "Density", D::Density)
        .help("Mass per volume; needed for dynamics, frequencies and gravity. Steel ≈ 7850 kg/m³, aluminium ≈ 2700 kg/m³."),
    PLASTIC,
    HYPERELASTIC,
    DAMAGE,
    EROSION,
    THERMAL,
    PVT,
    PHASE_CHANGE,
    variants(
        "fluid",
        "Fluid",
        &[
            variant(
                "newtonian",
                "Newtonian",
                &[positive("viscosity", "Viscosity", D::Viscosity).help("Dynamic viscosity. Water ≈ 1e-3 Pa·s (1e-9 MPa·s), air ≈ 1.8e-5 Pa·s.").req()],
            ),
            variant(
                "power_law",
                "Power law",
                &[
                    positive("consistency", "Consistency K", D::Viscosity)
                        .help("The viscosity at the reference shear rate (K in Pa·sⁿ is the viscosity at 1/s).")
                        .req(),
                    positive("index", "Index n", D::None).help("Flow behaviour index: below 1 shear-thinning, above 1 shear-thickening.").req(),
                    positive("reference_shear_rate", "Reference shear rate", D::StrainRate)
                        .help("The shear rate the consistency is the viscosity at.")
                        .hint("1/s"),
                ],
            ),
            variant(
                "carreau_yasuda",
                "Carreau–Yasuda",
                &[
                    positive("zero_shear_viscosity", "Zero-shear viscosity η₀", D::Viscosity).help("The viscosity at rest.").req(),
                    positive("infinite_shear_viscosity", "Infinite-shear viscosity η∞", D::Viscosity)
                        .help("The viscosity at very high shear rates, at most η₀.")
                        .hint("0"),
                    positive("time_constant", "Time constant λ", D::Time).help("1/λ is the shear rate where shear-thinning sets in.").req(),
                    positive("a", "Yasuda exponent a", D::None).help("How sharp the transition between the plateau and the power law is; 2 is Carreau's model.").hint("2"),
                    positive("index", "Index n", D::None).help("Power-law index of the shear-thinning part.").req(),
                ],
            ),
            variant(
                "cross",
                "Cross",
                &[
                    positive("zero_shear_viscosity", "Zero-shear viscosity η₀", D::Viscosity).help("The viscosity at rest.").req(),
                    positive("time_constant", "Time constant λ", D::Time).help("The viscosity is half η₀ at the shear rate 1/λ.").req(),
                    positive("index", "Index n", D::None).help("Power-law index of the shear-thinning part, between 0 and 1.").req(),
                ],
            ),
            variant(
                "cross_wlf",
                "Cross-WLF",
                &[
                    positive("d1", "D1", D::Viscosity).help("The zero-shear viscosity at the reference temperature D2.").req(),
                    num("d2", "D2", D::Temperature).help("The WLF reference temperature, usually the glass transition.").req(),
                    positive("a1", "A1", D::None).help("The WLF shift's first constant.").req(),
                    positive("a2", "A2", D::TemperatureDifference).help("The WLF shift's second constant; the melt is frozen at D2 − A2.").req(),
                    positive("critical_stress", "Critical stress τ*", D::Stress).help("The shear stress where shear-thinning sets in.").req(),
                    positive("index", "Index n", D::None).help("Power-law index of the shear-thinning part, between 0 and 1.").req(),
                ],
            ),
        ],
    )
    .help("Makes the material a fluid for a lattice Boltzmann flow analysis: a constant viscosity, or one that depends on the shear rate (and for Cross-WLF on the temperature, which needs heat in the fluid).")
    .hint("a fluid needs a density"),
];

pub const MATERIAL_PLASTICITY: &[Field] =
    &[reference("name", "Material", C::Material).help("The material to change.").req(), PLASTIC.req().def(r#"{"type": "bilinear"}"#)];
pub const MATERIAL_PVT: &[Field] =
    &[reference("name", "Material", C::Material).help("The material to change.").req(), PVT.req().def(r#"{"type": "tait"}"#)];
pub const MATERIAL_SPECIFIC_VOLUME: &[Field] = &[
    reference("name", "Material", C::Material).help("The material to evaluate.").req(),
    num("temperature", "Temperature", D::Temperature).help("The temperature of the state.").req(),
    positive("pressure", "Pressure", D::Stress).help("The pressure of the state.").def("0"),
];

pub const MATERIAL_PHASE_CHANGE: &[Field] = &[reference("name", "Material", C::Material).help("The material to change.").req(), PHASE_CHANGE.req()];
pub const MATERIAL_DAMAGE: &[Field] = &[
    reference("name", "Material", C::Material).help("The material to change.").req(),
    DAMAGE.def(r#"{"type": "johnson_cook"}"#).hint("give damage, erosion or both"),
    EROSION,
];

const PROFILE: Field = variants(
    "profile",
    "Profile",
    &[
        variant(
            "rectangle",
            "Rectangle",
            &[
                positive("width", "Width", D::Length).help("Cross-section width, along the section's first axis.").req(),
                positive("height", "Height", D::Length).help("Cross-section height, along the section's second axis.").req(),
            ],
        ),
        variant("circle", "Circle", &[positive("radius", "Radius", D::Length).help("Outer radius of the cross-section.").req()]),
        variant(
            "pipe",
            "Pipe",
            &[
                positive("radius", "Radius", D::Length).help("Outer radius of the pipe.").req(),
                positive("thickness", "Wall", D::Length).help("Wall thickness of the pipe.").req(),
            ],
        ),
        variant(
            "i",
            "I-beam",
            &[
                positive("height", "Height", D::Length).help("Overall height of the I-beam.").req(),
                positive("width", "Width", D::Length).help("Flange width.").req(),
                positive("web", "Web", D::Length).help("Web thickness.").req(),
                positive("flange", "Flange", D::Length).help("Flange thickness.").req(),
            ],
        ),
        variant(
            "box",
            "Box",
            &[
                positive("width", "Width", D::Length).help("Outer width of the box.").req(),
                positive("height", "Height", D::Length).help("Outer height of the box.").req(),
                positive("thickness", "Wall", D::Length).help("Wall thickness of the box.").req(),
            ],
        ),
    ],
)
.help("The beam's cross-section shape.")
.req()
.def(r#"{"type": "rectangle"}"#);

pub const SECTION: &[Field] = &[
    NAME,
    reference("set", "Region", C::Set).help("The set of elements the section applies to.").req(),
    reference("material", "Material", C::Material).help("The material of the section.").req(),
    variants(
        "kind",
        "Kind",
        &[
            variant("solid", "Solid", &[]),
            variant(
                "shell",
                "Shell",
                &[
                    positive("thickness", "Thickness", D::Length).help("Shell thickness.").req(),
                    range("offset", "Offset", D::None, -0.5, 0.5)
                        .help("Where the reference surface sits, as a fraction of the thickness: 0 mid-surface, ±0.5 top or bottom.")
                        .def("0")
                        .hint("fraction of the thickness"),
                    int("points", "Integration points", 1, 15)
                        .help("Points through the thickness. 5 is usual; more for plasticity through the thickness.")
                        .def("5"),
                ],
            ),
            variant(
                "beam",
                "Beam",
                &[
                    PROFILE,
                    xyz("orientation", "Orientation", D::None)
                        .help("A direction that sets the beam section's first axis (not along the beam).")
                        .req()
                        .def("[0, 0, 1]"),
                ],
            ),
            variant("rod", "Rod", &[positive("area", "Area", D::Area).help("Cross-section area of the rod.").req()]),
            variant(
                "plane_stress",
                "Plane stress",
                &[positive("thickness", "Thickness", D::Length).help("Out-of-plane thickness of the plane-stress elements.").req().def("1")],
            ),
            variant(
                "plane_strain",
                "Plane strain",
                &[positive("thickness", "Thickness", D::Length)
                    .help("Out-of-plane thickness of the plane-strain slice; results are per this length.")
                    .hint("1 mm")],
            ),
            variant("axisymmetric", "Axisymmetric", &[]),
        ],
    )
    .help("What the elements are: solids, shells, beams, rods or 2D elements.")
    .req()
    .def(r#"{"type": "solid"}"#),
    group(
        "formulation",
        "Formulation",
        &[
            choice("integration", "Integration", &[("full", "Full"), ("reduced", "Reduced")])
                .help("Full: every integration point. Reduced: fewer points, cheaper and less stiff, but may need hourglass control.")
                .def("\"full\""),
            flag("incompatible_modes", "Incompatible modes")
                .help("Extra bending modes for linear hexahedra, so a coarse hex8 mesh bends without locking.")
                .def("false"),
            flag("hybrid", "Hybrid").help("Pressure as an extra unknown, for nearly incompressible materials such as rubber.").def("false"),
        ],
    )
    .help("How the solid elements are integrated; the defaults suit most models."),
];

pub const REFERENCE_POINT: &[Field] = &[
    NAME,
    xyz("xyz", "Position", D::Length).help("The point's coordinates.").req().def("[0, 0, 0]"),
    int("label", "Node label", 0, u32::MAX as i64).help("The node label written to decks.").hint("0: the highest free label"),
];

const DOFS: Field = multi("dofs", "Degrees of freedom", &[("1", "UX"), ("2", "UY"), ("3", "UZ"), ("4", "URX"), ("5", "URY"), ("6", "URZ")])
    .help("Which of the six DOFs to couple (1–3 translations, 4–6 rotations); all when empty.")
    .hint("all");

pub const COUPLING: &[Field] = &[
    NAME,
    reference("reference_point", "Reference point", C::ReferencePoint).help("The reference point the region follows.").req(),
    reference("set", "Region", C::Set).help("The nodes or faces coupled to the reference point.").req(),
    choice("kind", "Kind", &[("kinematic", "Kinematic"), ("distributing", "Distributing")])
        .help("Kinematic: the region moves rigidly with the point. Distributing: forces spread over the region, which stays free to deform.")
        .req()
        .def("\"kinematic\""),
    DOFS,
];

pub const MASS: &[Field] = &[
    NAME,
    TARGET,
    per_axis("mass", "Mass", D::Mass, Some(0.0), XYZ)
        .help("A lumped mass at each node of the target: one value, or one per direction X, Y, Z.")
        .hint("per node"),
    vector("inertia", "Inertia", D::Inertia, &["I11", "I22", "I33", "I12", "I13", "I23"])
        .help("Rotary inertia about the reference point: three moments and three products.")
        .hint("about the node; needs a reference point"),
];

const SPRING_FIRST: Field = target("first", "First end").help("Where the spring starts: a set or a reference point.").req();
const SPRING_SECOND: Field = target("second", "Second end").help("Where the spring ends; ground when empty.").hint("ground when empty");
const SPRING_DIRECTION: Field = either(
    "direction",
    "Direction",
    &[
        int("dof", "DOF", 1, 6).help("1–3 translation along X, Y, Z; 4–6 rotation about them."),
        xyz("vector", "Vector", D::None).help("A direction vector; its length does not matter."),
        unit("line", "Along the ends").help("Along the line from the first end to the second."),
    ],
)
.help("The direction the spring acts in: one DOF, a vector, or the line between its ends.")
.req()
.def(r#"{"dof": 1}"#);

pub const SPRING: &[Field] = &[
    NAME,
    variants(
        "",
        "Type",
        &[
            variant(
                "spring",
                "Spring",
                &[
                    SPRING_FIRST,
                    SPRING_SECOND,
                    SPRING_DIRECTION,
                    positive("stiffness", "Stiffness", D::SpringStiffness)
                        .help("Force per unit stretch (moment per radian for a rotational spring).")
                        .req(),
                ],
            ),
            variant(
                "dashpot",
                "Dashpot",
                &[
                    SPRING_FIRST,
                    SPRING_SECOND,
                    SPRING_DIRECTION,
                    positive("coefficient", "Coefficient", D::SpringDamping).help("Damping force per unit velocity.").req(),
                ],
            ),
        ],
    )
    .help("What kind of connector: a spring (stiffness) or a dashpot (damping).")
    .req()
    .def("\"spring\""),
];

/// The target of a boundary condition, force, moment or temperature.
const TARGET: Field = target("target", "Target").help("Where it applies: a set or a reference point.").req();

pub const BOUNDARY_CONDITION: &[Field] = &[
    NAME,
    variants(
        "",
        "Type",
        &[
            variant("fixed", "Fixed", &[TARGET]),
            variant(
                "displacement",
                "Displacement",
                &[
                    TARGET,
                    free_vector("values", "Values", D::Motion, &["u1", "u2", "u3", "ur1", "ur2", "ur3"])
                        .help("Displacements u1–u3 along X, Y, Z and rotations ur1–ur3; empty entries stay free.")
                        .req()
                        .def("[0, 0, 0, null, null, null]")
                        .hint("empty: free"),
                ],
            ),
            variant(
                "symmetry",
                "Symmetry",
                &[
                    TARGET,
                    choice("plane", "Plane normal", &[("x", "X"), ("y", "Y"), ("z", "Z")])
                        .help("The axis normal to the symmetry plane: its translation and the two rotations about the other axes are held.")
                        .req()
                        .def("\"x\""),
                ],
            ),
            variant("temperature", "Temperature", &[TARGET, num("value", "Temperature", D::Temperature).help("The prescribed temperature.").req()]),
        ],
    )
    .help("What is held: everything fixed, given displacements, a symmetry plane or a temperature.")
    .req()
    .def("\"fixed\""),
];

const DISTRIBUTION: Field = choice("distribution", "Distribution", &[("total", "Total"), ("per_node", "Per node")])
    .help("Total: the value is split over the target's nodes. Per node: every node gets the full value.")
    .def("\"total\"");
const LOAD_SET: Field = reference("set", "Region", C::Set).help("The faces or elements the load acts on.").req();

pub const LOAD: &[Field] = &[
    NAME,
    variants(
        "",
        "Type",
        &[
            variant(
                "force",
                "Force",
                &[TARGET, xyz("vector", "Force", D::Force).help("The force, X, Y and Z components in the global axes.").req(), DISTRIBUTION],
            ),
            variant(
                "moment",
                "Moment",
                &[TARGET, xyz("vector", "Moment", D::Moment).help("The moment, about the X, Y and Z axes.").req(), DISTRIBUTION],
            ),
            variant(
                "pressure",
                "Pressure",
                &[
                    LOAD_SET,
                    num("magnitude", "Pressure", D::Stress)
                        .help("Pressure normal to the faces; positive pushes into the surface. 1 MPa = 1 N/mm².")
                        .req()
                        .hint("positive pushes in"),
                ],
            ),
            variant(
                "traction",
                "Traction",
                &[LOAD_SET, xyz("vector", "Traction", D::Stress).help("Force per area on the faces, X, Y and Z components.").req()],
            ),
            variant(
                "gravity",
                "Gravity",
                &[
                    xyz("acceleration", "Acceleration", D::Acceleration)
                        .help("Gravity or other acceleration. Earth gravity downward is 9.81 m/s² along −Z.")
                        .req(),
                    reference("set", "Region", C::Set).help("Where gravity acts; the whole model when empty.").hint("the whole model"),
                ],
            ),
            variant(
                "body_force",
                "Body force",
                &[LOAD_SET, xyz("vector", "Force per volume", D::ForcePerVolume).help("A body force per unit volume, X, Y and Z components.").req()],
            ),
            variant(
                "centrifugal",
                "Centrifugal",
                &[
                    num("angular_velocity", "Angular velocity", D::AngularVelocity)
                        .help("Spin speed in radians per time unit. 3000 rpm ≈ 314 rad/s.")
                        .req(),
                    xyz("point", "Point on axis", D::Length).help("A point on the axis of rotation.").req().def("[0, 0, 0]"),
                    xyz("axis", "Axis", D::None).help("Direction of the axis of rotation.").req().def("[0, 0, 1]"),
                    reference("set", "Region", C::Set).help("Where the centrifugal load acts; the whole model when empty.").hint("the whole model"),
                ],
            ),
            variant(
                "surface_flux",
                "Surface heat flux",
                &[LOAD_SET, num("flux", "Flux", D::HeatFlux).help("Heat flow per area into the surface.").req()],
            ),
            variant(
                "convection",
                "Convection",
                &[
                    LOAD_SET,
                    positive("film", "Film coefficient", D::FilmCoefficient)
                        .help("Convective heat transfer coefficient. Still air ≈ 5–25 W/(m²·K), water 500–10000.")
                        .req(),
                    num("sink", "Sink temperature", D::Temperature).help("The temperature of the surrounding fluid.").req(),
                ],
            ),
            variant(
                "heat_source",
                "Heat source",
                &[LOAD_SET, num("power", "Power per volume", D::HeatGeneration).help("Heat generated per unit volume.").req()],
            ),
            variant("temperature", "Temperature", &[TARGET, num("value", "Temperature", D::Temperature).help("The prescribed temperature.").req()]),
        ],
    )
    .help("What kind of load: forces, pressures, gravity, heat and more.")
    .req()
    .def("\"force\""),
];

const FACE_CONDITION: &[Variant] = &[
    variant("wall", "Wall", &[]),
    variant("moving_wall", "Moving wall", &[xyz("velocity", "Velocity", D::Velocity).help("The velocity of the moving wall.").req()]),
    variant("slip", "Slip", &[]),
    variant(
        "velocity",
        "Velocity inlet",
        &[
            xyz("velocity", "Velocity", D::Velocity).help("Inflow velocity.").req(),
            choice("profile", "Profile", &[("uniform", "Uniform"), ("parabolic", "Parabolic")])
                .help("Uniform across the inlet, or parabolic (fully developed flow).")
                .def("\"uniform\""),
        ],
    ),
    variant(
        "pressure",
        "Pressure outlet",
        &[
            num("pressure", "Pressure", D::Stress).help("Pressure at the outlet, relative to the domain's reference.").req().def("0"),
            choice("outflow", "Outflow", &[("uniform", "Uniform"), ("convective", "Convective")])
                .help(
                    "Uniform: the flow crosses the outlet normally. Convective: the flow leaves as it arrives, for outlets a wake or vortices reach.",
                )
                .def("\"uniform\""),
        ],
    ),
    variant("periodic", "Periodic", &[]),
];

const BOX_FACES: &[(&str, &str)] =
    &[("x_min", "X min"), ("x_max", "X max"), ("y_min", "Y min"), ("y_max", "Y max"), ("z_min", "Z min"), ("z_max", "Z max")];
const BOX_FACE: Field = choice("", "Box face", BOX_FACES).help("A face of the domain's box, at the smallest or largest X, Y or Z.").req();
const OBSTACLE_FACE_FIELDS: &[Field] = &[
    reference("body", "Body", C::Body).help("An obstacle body of the design.").req(),
    text("face", "Face").help("A persistent face name of the body (e.g. F3); it also selects the face's split pieces.").req(),
];
const OBSTACLE_FACE: Field = group("face", "Face", OBSTACLE_FACE_FIELDS).help("The obstacle face.").req();
const SURFACE_CONDITION: &[Variant] = &[
    variant("wall", "Wall at rest", &[]),
    variant(
        "moving_wall",
        "Moving wall",
        &[xyz("velocity", "Velocity", D::Velocity).help("The wall's tangential velocity; the body itself does not move.").req()],
    ),
];

const THERMAL_CONDITION: &[Variant] = &[
    variant("adiabatic", "Adiabatic", &[]),
    variant("temperature", "Temperature", &[num("value", "Temperature", D::Temperature).help("The temperature held at the boundary.").req()]),
    variant("flux", "Heat flux", &[num("value", "Heat flux", D::HeatFlux).help("Heat per area into the fluid (negative out of it).").req()]),
];

const fn face(key: &'static str, label: &'static str) -> Field {
    variants(key, label, FACE_CONDITION).hint("wall")
}

const fn thermal_face(key: &'static str, label: &'static str) -> Field {
    variants(key, label, THERMAL_CONDITION).hint("adiabatic")
}

const FLUID_THERMAL: Field = group(
    "thermal",
    "Heat",
    &[
        num("initial_temperature", "Initial temperature", D::Temperature)
            .help("The temperature everywhere at the start, and the reference of the buoyancy.")
            .hint("0 °C"),
        xyz("gravity", "Gravity", D::Acceleration).help("Gravity for the buoyancy; zero for none. Needs the expansion coefficient."),
        positive("expansion", "Expansion β", D::Expansion).help("Volumetric thermal expansion coefficient of the fluid, for the buoyancy."),
        group(
            "faces",
            "Box faces",
            &[
                thermal_face("x_min", "X min").help("Thermal condition on the box face at the smallest X."),
                thermal_face("x_max", "X max").help("Thermal condition on the box face at the largest X."),
                thermal_face("y_min", "Y min").help("Thermal condition on the box face at the smallest Y."),
                thermal_face("y_max", "Y max").help("Thermal condition on the box face at the largest Y."),
                thermal_face("z_min", "Z min").help("Thermal condition on the box face at the smallest Z."),
                thermal_face("z_max", "Z max").help("Thermal condition on the box face at the largest Z."),
            ],
        )
        .help("The thermal condition on each face of the box; a velocity inlet carries its held temperature in."),
        rows(
            "surfaces",
            "Obstacle faces",
            &[
                OBSTACLE_FACE,
                variants("condition", "Condition", THERMAL_CONDITION)
                    .help("The thermal condition on that face.")
                    .req()
                    .def(r#"{"type": "adiabatic"}"#),
            ],
        )
        .help("Thermal conditions on single obstacle faces; every other obstacle face is adiabatic."),
        rows(
            "solids",
            "Conducting solids",
            &[
                reference("body", "Body", C::Body).help("An obstacle body that conducts heat.").req(),
                reference("material", "Material", C::Material)
                    .help("Its material: a density, a conductivity (isotropic, or per axis of the box) and a specific heat.")
                    .req(),
                num("power", "Power", D::Power).help("The heat it generates, spread over its volume (a component's dissipation).").def("0"),
                num("limit", "Temperature limit", D::Temperature)
                    .help("The highest temperature it may reach (a junction or case limit); the solve reports its peak against it.")
                    .hint("none"),
            ],
        )
        .help("Obstacles heat conducts through, exchanging heat with the fluid across their faces (conjugate heat transfer); their faces take no thermal condition."),
        rows(
            "contacts",
            "Thermal contacts",
            &[
                references("bodies", "Bodies", C::Body).help("The two touching conducting solids.").req(),
                positive("resistance", "Resistance", D::Of(fieldcraft_model::units::Quantity::ThermalContactResistance))
                    .help("The interface layer's resistance per area (a paste's or pad's data sheet value); or give its thickness and conductivity."),
                positive("thickness", "Layer thickness", D::Length).help("The interface layer's thickness, with its conductivity."),
                positive("conductivity", "Layer conductivity", D::Conductivity).help("The interface layer's conductivity, with its thickness."),
            ],
        )
        .help("Thermal pastes, pads and gap fillers between touching conducting solids, not resolved in the lattice: each adds its resistance to the faces between the two bodies. Solids not listed touch perfectly."),
    ],
)
.help("Heat carried by the flow: needs the material's conductivity and specific heat. Off for an isothermal flow.")
.hint("isothermal when off");

const AXIS: &[(&str, &str)] = &[("\"x\"", "X"), ("\"y\"", "Y"), ("\"z\"", "Z")];

const RESISTANCE: &[Variant] = &[
    variant(
        "darcy_forchheimer",
        "Darcy–Forchheimer",
        &[
            positive("permeability", "Permeability K", D::Area).help("The viscous resistance: a force per volume of μ/K times the velocity.").req(),
            num("inertial", "Inertial C₂", D::InverseLength)
                .help(
                    "The inertial resistance per length: a force per volume of ρ C₂/2 times the velocity squared. Zero for a purely viscous medium.",
                )
                .def("0"),
        ],
    ),
    variant(
        "loss_coefficient",
        "Loss coefficient",
        &[positive("coefficient", "Loss coefficient k", D::None)
            .help("The pressure drop across the zone is k ρ u²/2 at the velocity u through it (a grille's or filter's data sheet value).")
            .req()],
    ),
];

const POROUS_ZONE: &[Field] = &[
    text("name", "Name").help("The vent's name, in the solve report.").req(),
    xyz("min", "Box min", D::Length).help("The zone's corner with the smallest coordinates.").req(),
    xyz("max", "Box max", D::Length).help("The zone's corner with the largest coordinates; it holds the cells whose centres it contains.").req(),
    choice("axis", "Flow axis", AXIS)
        .help("The axis the flow goes through it along: the thickness of a loss coefficient and the direction of its reported flow.")
        .req(),
    variants("resistance", "Resistance", RESISTANCE)
        .help("How the zone resists the flow (in every direction).")
        .req()
        .def(r#"{"type": "loss_coefficient"}"#),
];

const FAN: &[Field] = &[
    text("name", "Name").help("The fan's name, in the solve report.").req(),
    xyz("min", "Box min", D::Length).help("The fan's corner with the smallest coordinates.").req(),
    xyz("max", "Box max", D::Length)
        .help("The fan's corner with the largest coordinates: usually one or a few cells thick along its direction.")
        .req(),
    choice("direction", "Direction", &[("\"+x\"", "+X"), ("\"-x\"", "−X"), ("\"+y\"", "+Y"), ("\"-y\"", "−Y"), ("\"+z\"", "+Z"), ("\"-z\"", "−Z")])
        .help("The direction the fan blows.")
        .req(),
    table("curve", "Fan curve", D::FanCurve, &["flow", "pressure rise"], fieldcraft_model::fluid_zones::MAX_CURVE_POINTS)
        .help("The pressure–flow curve: volume flow (volume per time) increasing from 0, pressure rise not increasing.")
        .req(),
    range("relaxation", "Relaxation", D::None, MIN_FAN_RELAXATION, 1.0)
        .help("How far the pressure rise moves toward the curve each time step; smaller is steadier.")
        .hint("0.01"),
];

pub const FLUID_DOMAIN: &[Field] = &[
    NAME,
    reference("material", "Fluid", C::Material).help("The fluid material: needs a density and a viscosity.").req(),
    xyz("min", "Box min", D::Length).help("The box corner with the smallest coordinates.").req(),
    xyz("max", "Box max", D::Length).help("The box corner with the largest coordinates.").req(),
    positive("spacing", "Cell size", D::Length)
        .help("Lattice cell edge length. Smaller is more accurate and costs as the cube of the refinement.")
        .req(),
    group(
        "faces",
        "Box faces",
        &[
            face("x_min", "X min").help("Condition on the box face at the smallest X."),
            face("x_max", "X max").help("Condition on the box face at the largest X."),
            face("y_min", "Y min").help("Condition on the box face at the smallest Y."),
            face("y_max", "Y max").help("Condition on the box face at the largest Y."),
            face("z_min", "Z min").help("Condition on the box face at the smallest Z."),
            face("z_max", "Z max").help("Condition on the box face at the largest Z."),
        ],
    )
    .help("The condition on each face of the box."),
    references("obstacles", "Obstacles", C::Body).help("Bodies of the design the flow goes around."),
    rows(
        "surfaces",
        "Obstacle faces",
        &[
            OBSTACLE_FACE,
            variants("condition", "Condition", SURFACE_CONDITION).help("The flow condition on that face.").req().def(r#"{"type": "wall"}"#),
        ],
    )
    .help("Conditions on single obstacle faces; every other obstacle face is a wall at rest.")
    .hint("walls by default"),
    xyz("initial_velocity", "Initial velocity", D::Velocity).help("The fluid's velocity at the start."),
    xyz("acceleration", "Acceleration", D::Acceleration).help("A body acceleration on the fluid, such as gravity."),
    positive("tolerance", "Obstacle tolerance", D::Length)
        .help("How close a cell must be to an obstacle surface to count as on it.")
        .hint("a tenth of the cell size"),
    FLUID_THERMAL,
    rows("porous", "Vents and filters", POROUS_ZONE)
        .help("Porous resistance zones: vents, grilles, filters and perforated plates without resolving every hole."),
    rows("fans", "Fans", FAN).help("Fans: zones that raise the pressure by their pressure–flow curve; the solve finds the operating point."),
];

pub const FLUID_BOUNDARY: &[Field] = &[
    reference("domain", "Domain", C::FluidDomain).help("The fluid domain to change.").req(),
    name_or("at", "Boundary", &BOX_FACE, "Obstacle face", OBSTACLE_FACE_FIELDS)
        .help("A face of the domain's box, or a named face of an obstacle body.")
        .req()
        .def("\"x_min\""),
    variants("condition", "Condition", FACE_CONDITION).help("The flow condition on that boundary.").req().def(r#"{"type": "velocity"}"#),
];

pub const FLUID_THERMAL_BOUNDARY: &[Field] = &[
    reference("domain", "Domain", C::FluidDomain).help("The fluid domain to change; it needs heat settings.").req(),
    name_or("at", "Boundary", &BOX_FACE, "Obstacle face", OBSTACLE_FACE_FIELDS)
        .help("A face of the domain's box, or a named face of an obstacle body.")
        .req()
        .def("\"x_min\""),
    variants("condition", "Condition", THERMAL_CONDITION).help("The thermal condition on that boundary.").req().def(r#"{"type": "temperature"}"#),
];

pub const FLUID_COMPONENT_LIMIT: &[Field] = &[
    reference("domain", "Domain", C::FluidDomain).help("The fluid domain whose conducting solid it is.").req(),
    reference("body", "Body", C::Body).help("A conducting solid of the domain.").req(),
    num("limit", "Temperature limit", D::Temperature)
        .help("The highest temperature it may reach (a junction or case limit); empty for none.")
        .hint("none"),
];

const PROCEDURE: Field = variants(
    "procedure",
    "Procedure",
    &[
        variant("static", "Static, linear", &[]),
        variant(
            "nonlinear_static",
            "Static, nonlinear",
            &[
                flag("nonlinear_geometry", "Nonlinear geometry")
                    .help("Account for large displacements and rotations. Needed for buckling-like and large-deflection problems.")
                    .def("true"),
                range("initial_increment", "Initial increment", D::None, 0.0, 1.0)
                    .help("First increment, as a fraction of the step (0.1 = 10 %).")
                    .def("0.1"),
                range("min_increment", "Min increment", D::None, 0.0, 1.0).help("Smallest increment allowed before the step gives up.").def("1e-5"),
                range("max_increment", "Max increment", D::None, 0.0, 1.0).help("Largest increment, as a fraction of the step.").def("1.0"),
                int("max_increments", "Max increments", 1, 1_000_000).help("The most increments the step may take.").def("100"),
                int("max_iterations", "Iterations per increment", 1, 1000).help("Newton iterations before an increment is cut back.").def("16"),
                group(
                    "tolerances",
                    "Tolerances",
                    &[
                        positive("force", "Force", D::None).help("Residual force (moment), relative to the step's forces (moments).").def("1e-6"),
                        positive("displacement", "Displacement", D::None)
                            .help("Displacement (rotation) correction, relative to the increment's displacement (rotation).")
                            .def("1e-5"),
                        positive("energy", "Energy", D::None).help("Work of the residual, relative to the increment's work.").def("1e-10"),
                    ],
                )
                .help("Convergence tolerances; the defaults suit most models."),
                flag("line_search", "Line search")
                    .help("Scale each Newton correction to reduce the residual; helps with hard convergence.")
                    .def("false"),
                table("load_path", "Load path", D::None, &["fraction", "factor"], fieldcraft_model::MAX_LOAD_PATH_POINTS)
                    .help("How the load factor varies through the step: rows of step fraction (0 to 1) and load factor; a linear ramp when empty.")
                    .hint("a ramp when empty"),
                group(
                    "imperfection",
                    "Imperfection",
                    &[
                        reference("step", "Buckling step", C::Step).help("The solved buckling step whose mode is used.").req(),
                        int("mode", "Mode", 1, 100_000).help("Which buckling mode, from 1.").req().def("1"),
                        num("scale", "Largest component", D::Length)
                            .help("The largest displacement of the scaled mode, e.g. 1 % of the thickness.")
                            .req(),
                    ],
                )
                .help("Seed a buckling mode shape into the geometry, for post-buckling analysis."),
                group(
                    "arc_length",
                    "Arc length",
                    &[
                        choice("constraint", "Constraint", &[("cylindrical", "Cylindrical"), ("spherical", "Spherical")])
                            .help("Cylindrical: displacements only. Spherical: displacements and load factor.")
                            .def("\"cylindrical\""),
                        positive("max_load_factor", "Max load factor", D::None).help("Stop when the load factor reaches this.").def("1"),
                        group(
                            "displacement_limit",
                            "Displacement limit",
                            &[
                                int("node", "Node", 1, u32::MAX as i64).help("The label of the node watched.").req(),
                                int("direction", "Direction", 1, 3).help("1, 2, 3: its translation along X, Y, Z.").req().def("3"),
                                num("displacement", "Displacement", D::DofDisplacement)
                                    .help("The step stops when that translation reaches this size (its sign is ignored).")
                                    .req(),
                            ],
                        )
                        .help("Stop when a node's displacement reaches a value; at the max load factor when empty."),
                    ],
                )
                .help("Follow the load-displacement path through limit points (snap-through), with the load factor as an unknown."),
            ],
        ),
        variant(
            "buckle",
            "Buckling",
            &[
                int("modes", "Modes", 1, 100_000).help("How many buckling modes to find.").req().def("4"),
                num("shift", "Shift", D::None).help("Look for eigenvalues near this value; 0 finds the lowest."),
            ],
        ),
        variant(
            "frequency",
            "Frequency",
            &[
                int("modes", "Modes", 1, 10_000).help("How many natural frequencies to find.").req().def("10"),
                positive("min_hz", "Min frequency", D::Hertz).help("Lowest frequency of interest."),
                positive("max_hz", "Max frequency", D::Hertz).help("Highest frequency of interest."),
            ],
        ),
        variant(
            "explicit_dynamic",
            "Dynamic, explicit",
            &[
                positive("end_time", "End time", D::Time).help("Duration of the step. An impact is typically a few milliseconds.").req(),
                range("dt_scale", "Time step scale", D::None, 0.0, 1.0)
                    .help("Fraction of the stable time step used, 0 to 1; 0.9 is usual.")
                    .def("0.9"),
                int("max_increments", "Max increments", 1, 10_000_000_000).help("The most time increments the step may take.").def("10000000"),
                positive("output_interval", "Field output every", D::Time)
                    .help("Time between saved result frames; ten frames over the step when empty."),
                positive("history_interval", "History every", D::Time).help("Time between history samples."),
                reference("history_nodes", "History nodes", C::Set).help("The nodes whose displacements are recorded over time."),
                rows(
                    "initial_velocities",
                    "Initial velocities",
                    &[
                        reference("set", "Set", C::Set).help("The set given the velocity.").req(),
                        xyz("velocity", "Velocity", D::Velocity).help("Its starting velocity.").req(),
                    ],
                )
                .help("Starting velocities of sets, for impacts and drops."),
                rows(
                    "rigid_walls",
                    "Rigid walls",
                    &[
                        reference("set", "Set", C::Set).help("The nodes the wall stops.").req(),
                        xyz("point", "Point", D::Length).help("A point on the wall.").req(),
                        xyz("normal", "Normal", D::None).help("The wall's normal, pointing to the side the nodes stay on.").req(),
                        FRICTION,
                    ],
                )
                .help("Flat walls the set's nodes cannot pass."),
                group(
                    "hourglass",
                    "Hourglass control",
                    &[
                        choice("type", "Type", &[("viscous", "Viscous"), ("stiffness", "Stiffness")])
                            .help("Viscous: damping, for high-speed events. Stiffness: springs, for slower ones.")
                            .def("\"viscous\""),
                        positive("coefficient", "Coefficient", D::None).help("Strength of the hourglass control, typically 0.03–0.15."),
                    ],
                )
                .help("Suppresses zero-energy modes of reduced-integration elements."),
                group(
                    "bulk_viscosity",
                    "Bulk viscosity",
                    &[
                        positive("linear", "Linear", D::None).help("Linear bulk viscosity coefficient, typically 0.06."),
                        positive("quadratic", "Quadratic", D::None).help("Quadratic bulk viscosity coefficient, typically 1.2."),
                    ],
                )
                .help("Damps shock-wave ringing."),
                flag("nonlinear_geometry", "Nonlinear geometry").help("Account for large displacements and rotations.").def("true"),
                choice("stress_rate", "Stress rate", &[("jaumann", "Jaumann"), ("green_naghdi", "Green–Naghdi")])
                    .help("How stresses rotate with the material: Jaumann or Green–Naghdi."),
                group(
                    "mass_scaling",
                    "Mass scaling",
                    &[
                        positive("target_increment", "Target increment", D::Time).help("The time step mass scaling aims for.").req(),
                        reference("set", "Set", C::Set).help("Where mass is added; the whole model when empty."),
                    ],
                )
                .help("Adds mass to small elements so the time step reaches the target; check the added mass stays small."),
            ],
        ),
        variant(
            "heat_transfer",
            "Heat transfer",
            &[group(
                "transient",
                "Transient",
                &[
                    positive("end_time", "End time", D::Time).help("Duration of the transient.").req(),
                    positive("increment", "Increment", D::Time).help("Time step of the transient.").req(),
                    range("theta", "θ", D::None, 0.0, 1.0).help("Time integration: 0.5 Crank–Nicolson, 1 backward Euler (default, most stable)."),
                    num("initial_temperature", "Initial temperature", D::Temperature).help("Temperature at the start of the transient.").hint("0 °C"),
                ],
            )
            .help("Solve over time instead of for the steady state.")
            .hint("steady state when off")],
        ),
        variant("coupled_temperature_displacement", "Coupled temperature–displacement", &[]),
        variant(
            "fluid",
            "Fluid (lattice Boltzmann)",
            &[
                reference("domain", "Domain", C::FluidDomain).help("The fluid domain to solve.").req(),
                positive("end_time", "End time", D::Time).help("Simulated flow time, or give Steps.").hint("or Steps"),
                int("steps", "Steps", 1, 1_000_000_000_000).help("Lattice time steps, instead of End time."),
                positive("time_step", "Time step", D::Time)
                    .help("The lattice time step; chosen for a Mach number of 0.087 when empty.")
                    .hint("chosen for Mach 0.087"),
                int("frames", "Frames", 1, 100_000).help("How many result frames to keep.").def("10"),
                choice("collision", "Collision", &[("trt", "Two relaxation times"), ("regularized", "Regularized"), ("hrr", "Hybrid regularized")])
                    .help("The lattice Boltzmann collision. Regularized stays stable with τ close to ½ (fast flows on coarse grids, sharp corners), at a little more cost per step. Hybrid regularized blends in the strain rate of finite differences (σ), which damps the shortest waves.")
                    .def("\"trt\""),
                positive("magic", "TRT Λ", D::None).help("The two-relaxation-time magic parameter; 0.1875 is usual.").hint("0.1875"),
                range("sigma", "HRR σ", D::None, 0.0, 1.0)
                    .help("The hybrid collision's share of the populations' own stress, above 0 and at most 1 (1 is Regularized); the rest comes from finite differences of the velocity.")
                    .hint("0.99"),
                group(
                    "steady",
                    "Stop at steady state",
                    &[
                        range("tolerance", "Tolerance", D::None, 0.0, 1.0)
                            .help("The relative change per window: steady when, two windows in a row, every monitored quantity varies by less than this fraction of its largest value.")
                            .def("0.001"),
                        positive("window", "Window", D::Time)
                            .help("Flow time over which the change is measured; each window is sampled 8 times.")
                            .hint("a twentieth of the step"),
                        positive("min_time", "Minimum time", D::Time).help("No decision before this flow time.").hint("0"),
                        positive("divergence", "Divergence limit", D::None)
                            .help("Stop as diverged when a monitored quantity changes by more than this fraction over a window (after the minimum time).")
                            .hint("off"),
                        json("monitors", "Also monitor")
                            .help("Quantities watched besides the velocity (and temperature) residual, as JSON: mass_flow through a box face, force or heat on walls, mean_temperature.")
                            .hint("[{\"type\": \"force\", \"faces\": [\"x_min\"]}]"),
                    ],
                )
                .help("Stop on its own once the flow stops changing; End time or Steps is then the longest it runs.")
                .hint("runs its full length when off"),
            ],
        ),
    ],
)
.help("What the step solves: linear or nonlinear static, buckling, frequencies, explicit dynamics, heat transfer or flow.");

const OUTPUT: Field = multi(
    "output",
    "Output",
    &[
        ("\"displacement\"", "Displacement"),
        ("\"reaction_force\"", "Reaction force"),
        ("\"stress\"", "Stress"),
        ("\"strain\"", "Strain"),
        ("\"temperature\"", "Temperature"),
        ("\"heat_flux\"", "Heat flux"),
    ],
)
.help("Which results to keep; by procedure when none is ticked.")
.hint("by procedure");
const STEP_BCS: Field =
    references("boundary_conditions", "Boundary conditions", C::BoundaryCondition).help("The boundary conditions active in the step.");
const STEP_LOADS: Field = references("loads", "Loads", C::Load).help("The loads active in the step.");

pub const STEP_CREATE: &[Field] = &[
    NAME,
    PROCEDURE.req().def(r#"{"type": "static"}"#),
    OUTPUT,
    STEP_BCS,
    STEP_LOADS,
    reference("before", "Insert before", C::Step).help("The step the new one goes before; last when empty.").hint("last"),
];

pub const STEP_EDIT: &[Field] = &[reference("name", "Step", C::Step).help("The step to change.").req(), PROCEDURE, OUTPUT, STEP_BCS, STEP_LOADS];

pub const STEP_ACTIVATE: &[Field] = &[
    reference("step", "Step", C::Step).help("The step to change.").req(),
    references("boundary_conditions", "Boundary conditions", C::BoundaryCondition).help("The boundary conditions to switch on or off in the step."),
    references("loads", "Loads", C::Load).help("The loads to switch on or off in the step."),
    flag("active", "Active").help("On: they act in the step. Off: they are switched off from it on.").def("true"),
];

pub const DATUM_SUPPRESS: &[Field] = &[
    reference("name", "Datum", C::Datum).help("The datum to keep but leave out of the model.").req(),
    flag("suppressed", "Suppressed").help("On: left out (its partitions stop cutting). Off: back in the model.").def("true"),
];

pub const PARTITION_SUPPRESS: &[Field] = &[
    reference("name", "Partition", C::Partition).help("The partition to keep but leave out: its cells merge back.").req(),
    flag("suppressed", "Suppressed").help("On: the design is cut without it. Off: it cuts again.").def("true"),
];

pub const RENAME: &[Field] = &[
    choice(
        "kind",
        "Kind",
        &[
            ("set", "Set"),
            ("material", "Material"),
            ("section", "Section"),
            ("step", "Step"),
            ("boundary_condition", "Boundary condition"),
            ("load", "Load"),
            ("reference_point", "Reference point"),
            ("coupling", "Coupling"),
            ("contact", "Contact"),
            ("tie", "Tie"),
            ("point_mass", "Point mass"),
            ("spring", "Spring"),
            ("fluid_domain", "Fluid domain"),
            ("datum", "Datum"),
        ],
    )
    .help("What kind of definition to rename.")
    .req(),
    text("from", "Name").help("Its current name.").req(),
    text("to", "New name").help("Every definition that refers to it is changed to the new name.").req(),
];

pub const UNITS: &[Field] = &[
    choice(
        "system",
        "Unit system",
        &[
            ("mm-t-s-N-MPa", "mm, t, s, N, MPa, °C"),
            ("m-kg-s-N-Pa", "m, kg, s, N, Pa, K"),
            ("mm-kg-ms-kN-GPa", "mm, kg, ms, kN, GPa, °C"),
            ("mm-g-ms-N-MPa", "mm, g, ms, N, MPa, °C"),
            ("in-lbf-s-psi", "in, lbf·s²/in, s, lbf, psi, °F"),
            ("default", "My setting"),
        ],
    )
    .help("The units this model is shown and typed in. Every value is stored in SI, so changing them never changes the model: the numbers shown follow.")
    .req()
    .def("\"mm-t-s-N-MPa\""),
    choice("temperature", "Temperature", &[("°C", "°C"), ("K", "K"), ("°F", "°F")])
        .help("The temperature scale, when not the system's. Temperature differences and the per-degree units (expansion, conductivity, specific heat, film coefficient) follow it unless overridden too: with °F, ΔT is in °F and expansion per °F.")
        .hint("the system's"),
    json("overrides", "Overrides").help("Units for single quantities, by name: {\"energy\": \"J\", \"angle\": \"deg\"}.").hint("none"),
];

// ---- Contact ----

const FRICTION: Field = group(
    "friction",
    "Friction",
    &[
        positive("static", "Static μ", D::None)
            .help("Friction coefficient at rest. Steel on steel ≈ 0.15–0.8, rubber on concrete ≈ 0.6–1.")
            .hint("0"),
        positive("kinetic", "Kinetic μ", D::None)
            .help("Friction coefficient while sliding, usually a little lower than the static.")
            .hint("the static μ"),
        positive("decay", "Decay", D::InverseVelocity)
            .help("How fast friction falls from static to kinetic with slip rate.")
            .hint("0; per slip rate"),
    ],
)
.help("Coulomb friction between the surfaces; frictionless when off.")
.hint("frictionless");
const PENALTY: Field = range("penalty_scale", "Penalty scale", D::None, 0.0, 1.0)
    .help("Scales the contact stiffness. Lower softens contact and helps convergence, but lets surfaces interpenetrate more.")
    .def("1");
const PENETRATION: Field = choice("initial_penetration", "Initial penetration", &[("adjust", "Adjust nodes"), ("refuse", "Refuse")])
    .help("Adjust: move secondary nodes onto the main surface at the start. Refuse: report overlaps as errors.")
    .def("\"adjust\"");

pub const CONTACT_PAIR: &[Field] = &[
    NAME,
    reference("main", "Main surface", C::Set).help("The stiffer or coarser surface.").req(),
    reference("secondary", "Secondary surface", C::Set).help("The finer surface; its nodes are kept out of the main surface.").req(),
    FRICTION,
    PENALTY,
    PENETRATION,
];
pub const CONTACT_GENERAL: &[Field] = &[
    NAME,
    reference("surface", "Surface", C::Set).help("The faces that may touch each other; every exterior face when empty.").hint("every exterior face"),
    FRICTION,
    PENALTY,
    PENETRATION,
];
pub const CONTACT_EDIT: &[Field] = &[
    reference("name", "Contact", C::Contact).help("The contact to change.").req(),
    variants(
        "surfaces",
        "Surfaces",
        &[
            variant(
                "pair",
                "Pair",
                &[
                    reference("main", "Main", C::Set).help("The stiffer or coarser surface.").req(),
                    reference("secondary", "Secondary", C::Set).help("The finer surface.").req(),
                ],
            ),
            variant(
                "general",
                "General",
                &[reference("surface", "Surface", C::Set)
                    .help("The faces that may touch; every exterior face when empty.")
                    .hint("every exterior face")],
            ),
        ],
    )
    .help("A pair of surfaces, or general contact among many faces.")
    .req(),
    FRICTION,
    PENALTY,
    PENETRATION,
];
pub const TIE: &[Field] = &[
    NAME,
    reference("main", "Main surface", C::Set).help("The surface the secondary is tied to; usually the coarser mesh.").req(),
    reference("secondary", "Secondary surface", C::Set).help("The surface whose nodes are tied to the main surface.").req(),
    positive("position_tolerance", "Position tolerance", D::Length)
        .help("Secondary nodes farther than this from the main surface are left untied.")
        .hint("half the main element's thickness"),
];
pub const TIE_EDIT: &[Field] = &[
    reference("name", "Tie", C::Tie).help("The tie to change.").req(),
    reference("main", "Main surface", C::Set).help("The surface the secondary is tied to; usually the coarser mesh.").req(),
    reference("secondary", "Secondary surface", C::Set).help("The surface whose nodes are tied to the main surface.").req(),
    positive("position_tolerance", "Position tolerance", D::Length)
        .help("Secondary nodes farther than this from the main surface are left untied.")
        .hint("half the main element's thickness"),
];
pub const CONTACT_FIND: &[Field] = &[
    positive("tolerance", "Tolerance", D::Length).help("Faces closer than this count as touching.").hint("a tenth of the mean face size"),
    flag("create", "Create pairs").help("Also create a contact for each touching pair found.").def("false"),
    FRICTION,
    PENALTY,
    PENETRATION,
];

// ---- Meshing ----

const BODY: Field =
    reference("body", "Body", C::Body).help("The body to mesh; every body of the design when empty.").hint("every body of the design");
const SURFACE: Field =
    int("surface", "Imported surface", 0, 1_000_000).help("An imported surface mesh to mesh instead of a design body.").hint("instead of a body");
const TOLERANCE: Field = positive("tolerance", "Chordal tolerance", D::Length)
    .help("Largest distance between the tessellation and the true curved surface. Smaller follows curves closer, with more elements.")
    .hint("1e-3 of the diagonal");
const WELD: Field = range("weld", "Weld", D::None, 0.0, fieldcraft_mesh::surface::MAX_RELATIVE_WELD)
    .help("Points closer than this fraction of the diagonal are merged when the surface is assembled. For a confirmed tiny seam, choose less than features to preserve; this does not repair inconsistent face orientation.")
    .hint("1e-9 of the diagonal");
const WELD_DISTANCE: Field = positive("weld_distance", "Weld distance", D::Length)
    .help("Absolute distance for merging points when the surface is assembled; the larger of this and Weld times the diagonal is used. For a confirmed tiny seam, choose less than features to preserve. The default kernel tolerance does not cover every boolean retry displacement; welding does not repair inconsistent face orientation.")
    .hint("the surface's own");
const FEATURE_ANGLE: Field = range("feature_angle", "Feature angle", D::Degree, 0.0, 180.0)
    .help("Edges where faces meet at more than this angle are kept as sharp features.")
    .hint("45");

pub const TET: &[Field] = &[
    BODY,
    SURFACE,
    positive("size", "Element size", D::Length)
        .help("Target element edge length. Several elements through the thinnest wall are needed for accurate stresses.")
        .hint("1/20 of the diagonal"),
    choice("order", "Order", &[("linear", "Linear (tet4)"), ("quadratic", "Quadratic (tet10)")])
        .help("Linear (tet4) is fast but stiff in bending; quadratic (tet10) is far more accurate for stress.")
        .def("\"linear\""),
    TOLERANCE,
    flag("quality", "Improve quality").help("Smooth and swap elements after meshing to remove slivers.").def("true"),
    range("quality_target", "Quality target", D::Degree, 5.0, 60.0).help("The smallest dihedral angle quality improvement aims for.").hint("40"),
    range("curvature_angle", "Curvature angle", D::Degree, 0.0, 90.0)
        .help("Angle an element may span on a curved surface; smaller refines curved regions.")
        .hint("22.5"),
    range("grading", "Grading", D::None, 0.05, 2.0).help("How fast element size may grow from small to large; smaller is smoother.").hint("0.3"),
    int("max_passes", "Improvement passes", 0, 200).help("The most quality-improvement passes.").hint("20"),
    int("max_points", "Max points", 1, 20_000_000).help("Cap on the points the mesher creates.").hint("1000000"),
    int("max_steiner", "Max Steiner points", 0, 10_000_000).help("Cap on the points added to recover the surface.").hint("50000"),
    int("max_tets", "Max tetrahedra", 1, 50_000_000).help("Cap on the tetrahedra created.").hint("5000000"),
    WELD,
    WELD_DISTANCE,
    FEATURE_ANGLE,
];

const HEX_SIZE: Field = positive("size", "Element size", D::Length).help("Target element edge length.").hint("1/25 of the diagonal");
const HEX_ORDER: Field =
    int("order", "Order", 1, 2).help("1: linear hex8. 2: quadratic hex20, more accurate in bending.").def("1").hint("1: hex8, 2: hex20");

pub const HEX_GENERATE: &[Field] = &[
    BODY,
    SURFACE,
    HEX_SIZE,
    HEX_ORDER,
    choice(
        "method",
        "Method",
        &[("auto", "Automatic"), ("swept", "Any sweep"), ("extrude", "Extrude"), ("sweep", "Sweep"), ("revolve", "Revolve"), ("mapped", "Mapped")],
    )
    .help(
        "How to build the hexes. Automatic tries, in order, an extrusion, a sweep, a revolution and a mapped block, and \
         uses the first that fits; Any sweep tries the first three.",
    )
    .notes(HEX_METHOD_NOTES)
    .def("\"auto\""),
    TOLERANCE,
    WELD,
    WELD_DISTANCE,
    FEATURE_ANGLE,
];

/// What each structured hex method needs of the body, one line each.
const HEX_METHOD_NOTES: &[(&str, &str)] = &[
    ("auto", "Automatic: the first of Extrude, Sweep, Revolve and Mapped that fits the body."),
    ("swept", "Any sweep: the first of Extrude, Sweep and Revolve that fits the body."),
    ("extrude", "Extrude: the body must be one cross-section carried along a straight direction (levels allowed)."),
    ("sweep", "Sweep: a source face carried to a target face through four-sided linking faces."),
    ("revolve", "Revolve: the body must be one cross-section turned around an axis."),
    ("mapped", "Mapped: the body must have six sides, each with four edges (curved faces allowed)."),
];

pub const HEX_SWEEP: &[Field] = &[
    BODY,
    SURFACE,
    HEX_SIZE,
    HEX_ORDER,
    choice("method", "Method", &[("swept", "Any sweep"), ("extrude", "Extrude"), ("sweep", "Sweep"), ("revolve", "Revolve")])
        .help("How to sweep: any sweep, extrusion, sweep along a path or revolution.")
        .notes(HEX_METHOD_NOTES)
        .def("\"swept\""),
    TOLERANCE,
    WELD,
    WELD_DISTANCE,
    FEATURE_ANGLE,
];

pub const HEX_MAPPED: &[Field] = &[BODY, SURFACE, HEX_SIZE, HEX_ORDER, TOLERANCE, WELD, WELD_DISTANCE, FEATURE_ANGLE];

pub const HEX_BLOCKS: &[Field] = &[
    references("bodies", "Cells", C::Body)
        .help("The cells (bodies) to mesh as blocks; every body of the design when empty.")
        .hint("every body of the design"),
    HEX_SIZE,
    HEX_ORDER,
    TOLERANCE,
    WELD,
    WELD_DISTANCE,
    FEATURE_ANGLE,
];

const CELL_METHODS: &[(&str, &str)] = &[
    ("auto", "Auto (hex, else tet)"),
    ("hex", "Hex"),
    ("swept", "Swept"),
    ("extrude", "Extrude"),
    ("sweep", "Sweep"),
    ("revolve", "Revolve"),
    ("mapped", "Mapped"),
    ("tet", "Tet"),
];

pub const CELLS: &[Field] = &[
    references("bodies", "Cells", C::Body)
        .help("The cells (bodies) to mesh: tick them, or pick them in the viewport (Select: Bodies) and Use selection. Every body of the design when none is ticked.")
        .hint("every body of the design"),
    ref_choices("methods", "Methods", C::Body, CELL_METHODS).help("A meshing method for single cells; the others take Other cells' method."),
    choice(
        "default",
        "Other cells",
        &[("auto", "Auto (hex, else tet)"), ("hex", "Hex"), ("swept", "Swept"), ("mapped", "Mapped"), ("tet", "Tet"), ("blocks", "All hex (conforming blocks)")],
    )
    .help("The method for cells not listed in Methods. All hex meshes every cell as a mapped block, with the division counts of all cells assigned together so shared faces match (mesh.hex.blocks).")
    .notes(&[("blocks", "All hex: every cell must have six sides, each with four edges; leave Methods empty.")])
        .def("\"auto\""),
    HEX_SIZE,
    HEX_ORDER,
    TOLERANCE,
    choice("transition", "Hex–tet joint", &[("pyramid", "Pyramids (conforming)"), ("tie", "Tie")])
        .help("How a tet cell meets the hex cells across its partitions: a layer of pyramids on the hex faces, one conforming mesh (linear elements), or a tie between the two meshes.")
        .def("\"pyramid\""),
];

pub const HEX_PROJECT: &[Field] = &[
    reference("body", "Body", C::Body).help("The body to mesh with a projected block topology.").req(),
    HEX_SIZE,
    choice(
        "template",
        "Template",
        &[("grid", "Block grid (carve)"), ("sphere", "Sphere (7 blocks)"), ("butterfly", "Butterfly (5-block disc, extruded)"), ("blocks", "Blocks (given vertices)")],
    )
    .help("The block topology projected onto the body: a grid of blocks you may carve, the seven-block sphere, the five-block butterfly for round sections, or blocks you list.")
    .def("\"grid\""),
    xyz("center", "Centre", POINT).help("Sphere or butterfly centre.").hint("from the body's bounding box"),
    positive("radius", "Radius", D::Length).help("Sphere or butterfly radius.").hint("from the body"),
    xyz("axis", "Axis", D::None).help("Butterfly axis direction.").hint("0, 0, 1"),
    json("levels", "Levels")
        .numbers(&[("", D::Length)])
        .help("Sphere and butterfly levels; for a butterfly, distances along the axis from the centre, as JSON.")
        .hint("[z0, z1, ...]; from the bounding box"),
    range("inset", "Inset", D::None, 0.05, 0.95).help("Inner square or cube corner at this fraction of the radius.").def("0.55"),
    json("grid", "Grid").numbers(&[("/planes/*", D::Length)]).help("The block grid: planes per axis or a cell count, and blocks to carve away, as JSON.").hint(r#"{"cells": [ni, nj, nk], "carve": [{"i": [a, b], "j": [c, d], "k": [e, f]}]}"#),
    json("blocks", "Blocks").numbers(&[("/vertices", D::Length)]).help("Template Blocks: vertices and hex8 blocks, as JSON.").hint(r#"{"vertices": [[x, y, z], ...], "blocks": [[v1..v8], ...]}"#),
    json("operations", "Operations")
        .numbers(&[("/*/offset", D::Length), ("/*/extrude/offset", D::Length), ("/*/distance", D::Length)])
        .help("O-grids, extrusions, connections, removals and merges applied in order, as JSON.").hint(r#"[{"op": "o_grid", "block": 0, "open": ["+z"]}]"#),
    json("faces", "Face targets").help("Block faces sent to chosen B-rep faces; other outer faces go to the nearest surfaces.").hint(r#"[{"block": 0, "face": "+x", "to": ["F1:+x"]}]"#),
    json("divisions", "Divisions").help("Element counts per block direction.").hint(r#"[{"block": 0, "direction": 1, "count": 8}]"#),
    json("fixed", "Fixed faces")
        .help("Faces of the body another body's mesh already covers: its element faces there are kept as they are and this mesh meets them node to node.")
        .hint(r#"["F3:partition:P:b"]"#),
    choice("smoothing", "Smoothing", &[("elliptic", "Elliptic"), ("elliptic-controlled", "Elliptic, controlled (Thomas–Middlecoff)"), ("laplace", "Laplace"), ("none", "None")])
        .help("How interior nodes are smoothed after projection, with the boundary held on the geometry.")
        .def("\"elliptic\""),
    int("iterations", "Iterations", 1, 100_000).help("Smoothing iterations.").def("100"),
    choice("backend", "Backend", &[("cpu", "CPU"), ("coloured", "Coloured (threads)"), ("gpu", "GPU")])
        .help("Where the smoothing sweeps run: the CPU node by node (the same mesh on every machine), colour by colour on every core, or the GPU FIELDCRAFT_GPU picks (single precision, the same verdicts as coloured, smallest scaled Jacobian within 0.02 of it).")
        .def("\"cpu\""),
    range("min_quality", "Minimum quality", D::None, 0.0, 1.0).help("Smallest scaled Jacobian accepted.").def("0.2"),
    range("max_arc_angle", "Max arc angle", D::Degree, 0.0, 90.0)
        .help("The most a curved face may turn across one element; the divisions grow to keep it. 0 for no limit.")
        .def("30"),
    range("merge", "Merge distance", D::Length, 0.0, MAX_MERGE)
        .help("Join this mesh's nodes with other bodies' meshes within this distance; 0 never.")
        .hint("1e-6 of the diagonal"),
    flag("preview", "Preview only").help("Return the block topology without meshing."),
    TOLERANCE,
];

pub const CLEAR_MESH: &[Field] = &[references("bodies", "Bodies", C::Body)
    .help("The bodies whose meshes go. When none is ticked: the meshes of bodies the open design no longer has because its partitions changed (model.check lists them).")
    .hint("the bodies partitions removed")];

const MIRROR_DATUM: Field = reference("", "Datum plane", C::DatumPlane).help("A datum plane of the model.");
const MIRROR_POINT_NORMAL: &[Field] = &[
    xyz("point", "Point", D::Length).help("A point on the plane.").req(),
    xyz("normal", "Normal", D::None).help("The plane's normal; its length does not matter.").req(),
];

pub const MIRROR: &[Field] = &[
    name_or("plane", "Plane", &MIRROR_DATUM, "Point and normal", MIRROR_POINT_NORMAL)
        .help("The mirror plane: a datum plane, or a point on it and its normal.")
        .req(),
    references("bodies", "Bodies", C::Body).help("The meshed bodies to mirror; every meshed body when empty.").hint("every meshed body"),
    reference_text("as", "Images", C::Body)
        .help("The body each mirrored body becomes; when empty, the design's body at the reflected place, else <body>/mirror.")
        .hint("default"),
];

pub const SYMMETRY: &[Field] = &[
    NAME,
    reference("datum", "Datum plane", C::DatumPlane).help("The symmetry plane.").req(),
    reference("set", "Set", C::Set).help("The faces on the plane; found from the design when empty.").hint("the faces of the design on the plane"),
    references("bodies", "Bodies", C::Body).help("The bodies cut by the symmetry plane; every body when empty.").hint("every body"),
];

pub const HEX_OCTREE: &[Field] = &[
    BODY,
    SURFACE,
    HEX_SIZE,
    TOLERANCE,
    int("max_depth", "Refinement levels", 0, 5)
        .help("The most times the octree may refine near thin walls, narrow faces and curves; the levels the part needs are found before meshing.")
        .hint("4"),
    int("max_elements", "Max elements", 1, 10_000_000).help("Cap on the elements created."),
    int("attempts", "Attempts", 1, 8).help("How many times to retry with adjusted settings.").hint("4"),
    int("max_sweeps", "Optimiser sweeps", 0, 10_000).help("The most smoothing sweeps.").hint("300"),
    WELD,
    WELD_DISTANCE,
    FEATURE_ANGLE,
];

/// A face by body and persistent name, a row of a face list.
const FACE_ROW: &[Field] = &[
    reference("body", "Body", C::Body).help("A body of the open design.").req(),
    text("face", "Face").help("A persistent face name of the body (e.g. F3).").req().hint("persistent face name"),
];
const fn face_rows(key: &'static str, label: &'static str) -> Field {
    rows(key, label, FACE_ROW)
}

pub const HEX_EXTRUDE: &[Field] = &[
    face_rows("source", "Source")
        .help("The section: faces already meshed (by this body or a neighbour, whose nodes are shared), or one planar face meshed here. Pick it, then Use as source.")
        .hint("pick the source face"),
    reference("source_set", "Source set", C::Set).help("Instead of Source: a set of quadrilaterals (a 2D mesh, or hex element faces)."),
    face_rows("to_surface", "Target surface")
        .help("The faces the last layer lies on. Pick them, then Use as target surface.")
        .req()
        .hint("pick the target surface"),
    face_rows("sides", "Side faces").help("The faces the section's sides follow (they need not be straight); straight paths when empty.").hint("straight sides"),
    xyz("direction", "Direction", D::None).help("The direction every node moves along; along the section's normal toward the target when empty.").hint("the section's normal"),
    flag("normal", "Along the normal").help("Each node along the section's normal toward the target (the default when no direction is given)."),
    int("layers", "Layers", 1, 10_000).help("The number of element layers; from Size when empty.").hint("from the size"),
    positive("size", "Size", D::Length).help("Layer thickness, and the element size of a section meshed here.").hint("1/25 of the body"),
    range("bias", "Bias", D::None, 1.0 / fieldcraft_mesh::hex::MAX_BIAS, fieldcraft_mesh::hex::MAX_BIAS).help("The last layer's thickness over the first's: below 1 grades toward the target, above 1 toward the source.").def("1"),
    reference("body", "Body", C::Body).help("The body the elements fill; its mesh is replaced.").hint("the source face's body"),
    TOLERANCE,
];

/// A block's edge directions, from corner 1.
const BLOCK_DIRECTIONS: &[&str] = &["1-2", "1-4", "1-5"];

pub const BLOCK: &[Field] = &[
    text("body", "Body").help("Name of the new block body.").req().def("\"Block\""),
    choice(
        "kind",
        "Element",
        &[
            ("hex8", "hex8"),
            ("hex20", "hex20"),
            ("tet4", "tet4"),
            ("tet10", "tet10"),
            ("wedge6", "wedge6"),
            ("quad4", "quad4 (2D)"),
            ("quad8", "quad8 (2D)"),
            ("tri3", "tri3 (2D)"),
            ("tri6", "tri6 (2D)"),
        ],
    )
    .help("The element type of the block.")
    .req()
    .def("\"hex8\""),
    points("corners", "Corners", D::Length, 8).help("The eight corner points (four for 2D, at z = 0).").req(),
    tuple("divisions", "Divisions", BLOCK_DIRECTIONS, 2, true)
        .help("Elements along each edge direction: corner 1 to 2, 1 to 4 and 1 to 5; a 2D block takes the first two.")
        .req(),
    tuple("bias", "Bias", BLOCK_DIRECTIONS, 2, false)
        .help("Last over first element size along each direction; 1 is uniform. A 2D block takes the first two.")
        .hint("uniform"),
    json("edges", "Curved and named edges").help("Edge curves and names, as JSON."),
    json("faces", "Face names").help("Names of the block's faces, as JSON.").hint("{\"1\": \"bottom\"}"),
    json("vertices", "Corner names").help("Names of the block's corners, as JSON."),
    text("face", "2D face name").help("Name of a 2D block's face.").hint("area"),
    positive("merge", "Merge distance", D::Length).help("Nodes closer than this merge with existing nodes.").hint("1e-6 of the diagonal"),
];

pub const PLANAR: &[Field] = &[
    positive("size", "Element size", D::Length).help("Target element edge length.").req(),
    choice("kind", "Element", &[("quad4", "quad4"), ("quad8", "quad8"), ("tri3", "tri3"), ("tri6", "tri6")])
        .help("The 2D element type.")
        .def("\"quad4\""),
    BODY,
    text("face", "Face").help("A planar face of the body to mesh.").hint("a planar face of the body"),
    TOLERANCE,
    json("polygon", "Polygon")
        .numbers(&[("", D::Length)])
        .help("An outline to mesh instead of a face, as JSON points.")
        .hint("[[x, y], …] instead of a face"),
    json("holes", "Holes").numbers(&[("", D::Length)]).help("Holes in the polygon, as JSON point lists."),
    text("name", "Body name").help("Name of the new 2D body.").hint("Plane"),
    int("max_elements", "Max elements", 1, 4_000_000).help("Cap on the elements created."),
];

pub const QUALITY: &[Field] = &[
    labels("elements", "Elements", LabelOf::Element).help("The elements to measure; all when empty.").hint("all"),
    int("bins", "Histogram bins", 1, 1000).help("Number of bars in the quality histogram.").def("10"),
];

pub const QUALITY_MARK: &[Field] = &[
    labels("elements", "Elements", LabelOf::Element).help("The elements to check; all when empty.").hint("all"),
    group(
        "limits",
        "Limits",
        &[
            num("scaled_jacobian", "Scaled Jacobian at least", D::None)
                .help("Elements with a lower scaled Jacobian are marked (1 is a perfect element).")
                .hint("0.2"),
            num("min_dihedral", "Min dihedral at least", D::Degree).help("Tetrahedra with a smaller dihedral angle are marked.").hint("10"),
            num("max_dihedral", "Max dihedral at most", D::Degree).help("Tetrahedra with a larger dihedral angle are marked.").hint("165"),
            num("aspect_ratio", "Aspect ratio at most", D::None).help("Elements with a larger aspect ratio are marked.").hint("10"),
            num("min_angle", "Min angle at least", D::Degree).help("Shells with a smaller corner angle are marked.").hint("15"),
            num("max_angle", "Max angle at most", D::Degree).help("Shells with a larger corner angle are marked.").hint("165"),
        ],
    )
    .help("Elements worse than a limit are marked; a limit left empty keeps its default."),
];

pub const SURFACE_CHECK: &[Field] = &[
    reference("body", "Body", C::Body).help("The body whose surface to check.").req().hint("or an imported surface"),
    SURFACE,
    TOLERANCE,
    WELD,
    WELD_DISTANCE,
    FEATURE_ANGLE,
];

pub const NODE_ADD: &[Field] = &[
    int("label", "Label", 1, u32::MAX as i64).help("The node label.").req(),
    xyz("xyz", "Position", D::Length).help("The node's coordinates.").req().def("[0, 0, 0]"),
];
pub const LABEL: &[Field] = &[int("label", "Label", 1, u32::MAX as i64).help("The label of the node or element.").req()];
pub const ELEMENT_ADD: &[Field] = &[
    int("label", "Label", 1, u32::MAX as i64).help("The element label.").req(),
    choice(
        "kind",
        "Kind",
        &[
            ("rod2", "rod2"),
            ("beam2", "beam2"),
            ("tri3", "tri3"),
            ("tri6", "tri6"),
            ("quad4", "quad4"),
            ("quad8", "quad8"),
            ("tet4", "tet4"),
            ("tet10", "tet10"),
            ("wedge6", "wedge6"),
            ("pyramid5", "pyramid5"),
            ("hex8", "hex8"),
            ("hex20", "hex20"),
        ],
    )
    .help("The element type.")
    .req()
    .def("\"hex8\""),
    labels("nodes", "Nodes", LabelOf::Node).help("The element's node labels in order.").req().hint("labels in order: 1, 2, 3, 4"),
];

const MERGE_TOLERANCE: Field =
    positive("tolerance", "Tolerance", D::Length).help("Nodes closer than this are one node.").hint("a millionth of the mesh's size");

pub const NODE_MERGE: &[Field] = &[
    MERGE_TOLERANCE,
    reference("set", "Set", C::Set).help("Merge only the nodes of this set (and of its elements and faces).").hint("every node"),
    choice("keep", "Keep", &[("lowest", "Lowest label"), ("highest", "Highest label")]).help("Which node of each group survives.").def("\"lowest\""),
];

pub const CHECK_DUPLICATES: &[Field] = &[MERGE_TOLERANCE];

pub const NODE_RELAX: &[Field] = &[
    xyz("center", "Centre", D::Length).help("The centre of the sphere of influence; click the mesh to place it.").hint("click the mesh"),
    positive("radius", "Radius", D::Length).help("Nodes closer to the centre than this move; Alt + scroll over the view changes it."),
    int("around_element", "Around element", 1, 99_999_999)
        .help("Instead of a centre: a sphere round this element (radius default twice its longest edge).")
        .hint("none"),
    face_rows("faces", "Faces")
        .help("Instead of a sphere: every node on these CAD faces and the elements touching them, nodes sliding on the faces. Shift-click faces in the view, then Use selection.")
        .hint("Shift-click faces"),
    labels("nodes", "Nodes", LabelOf::Node).help("Instead of a sphere: these nodes.").hint("none"),
    labels("elements", "Elements", LabelOf::Element).help("Instead of a sphere: the nodes of these elements.").hint("none"),
    references("bodies", "Bodies", C::Body).help("Only the nodes of these bodies.").hint("every body"),
    choice("mode", "Mode", &[("guarded", "Guarded (local)"), ("strong", "Strong (optimise the region)")])
        .help("Guarded: each move keeps the worst element around the node from getting worse. Strong: move the nodes of the bad elements and morph the mesh round them, nodes on CAD edges held.")
        .def("\"guarded\""),
    range("below", "Below", D::None, 0.0, 1.0).help("Strong: the elements below this scaled Jacobian are moved.").hint("0.7"),
    positive("reach", "Reach", D::Length).help("Strong: nodes this close to a moved node follow it, blended to zero at the reach.").hint("five element lengths"),
    choice("guards", "Keep a morph", &[("valid", "Unless an element inverts"), ("strict", "Only when it improves the elements")])
        .help("Strong: valid keeps every morph that inverts no element; strict keeps one only when the worst element it touches improves and no element gets longer.")
        .hint("valid"),
    int("iterations", "Iterations", 1, 200).help("Guarded: sweeps over the nodes (1 to 100). Strong: morph rounds.").hint("5 (strong: 20)"),
    range("strength", "Strength", D::None, 0.0, 1.0).help("How far each node moves toward its relaxed position, 0 to 1.").def("1"),
    choice("falloff", "Falloff", &[("smooth", "Smooth"), ("linear", "Linear"), ("none", "None")])
        .help("How the strength falls from the centre to the radius, so the edit blends into the mesh.")
        .def("\"smooth\""),
    text("field", "Keep as").help("The relaxation field to keep, to apply again after a remesh.").hint("Relax1, Relax2, …"),
    flag("reapply", "Reapply after remesh").help("Apply the field again by itself after each remesh.").def("false"),
];

pub const NODE_MOVE: &[Field] = &[
    labels("nodes", "Nodes", LabelOf::Node).help("The nodes to grab: drag a node in the view, or box-select nodes and Use selection.").hint("drag a node"),
    face_rows("faces", "Faces").help("Instead of nodes: every node on these CAD faces.").hint("none"),
    xyz("center", "Centre", D::Length).help("Instead of nodes: a sphere's centre; every node inside follows by the falloff.").hint("none"),
    xyz("by", "Move by", D::Length).help("The displacement of the grabbed nodes; dragging in the view sets it."),
    xyz("to", "Move to", D::Length).help("Instead of Move by: where the grabbed nodes' centroid (or the centre) goes.").hint("none"),
    positive("radius", "Radius", D::Length)
        .help("Nodes this close to a grabbed node follow it, blended to zero at the radius; Alt + scroll over the view changes it.")
        .hint("only the grabbed nodes"),
    choice("propagate", "Followers", &[("falloff", "Falloff (by distance)"), ("harmonic", "Harmonic (through the mesh)")])
        .help("Falloff: a node follows by the falloff of its distance. Harmonic: the move spreads through the element edges as a membrane carries it, zero at the radius.")
        .def("\"falloff\""),
    choice("falloff", "Falloff", &[("smooth", "Smooth"), ("linear", "Linear"), ("none", "None")])
        .help("Falloff followers: how the share falls from the grabbed nodes to the radius.")
        .hint("smooth"),
    range("strength", "Strength", D::None, 0.0, 1.0).help("Falloff followers: the largest share of the move a follower takes, 0 to 1.").hint("1"),
    flag("on_geometry", "Stay on the geometry")
        .help("Nodes on CAD faces slide on them, on CAD edges along them; nodes on CAD vertices stay.")
        .def("true"),
    choice("guard", "Guard", &[("no_inversion", "Scale back to no inverted element"), ("none", "None")])
        .help("No inversion: the move is scaled back until no element inverts or falls below the smallest quality.")
        .def("\"no_inversion\""),
    range("min_quality", "Smallest quality", D::None, 0.0, 1.0)
        .help("No inversion: the smallest scaled Jacobian the move may leave an element with (one already below keeps its own).")
        .def("0.05"),
];

pub const NODE_RELAX_APPLY: &[Field] = &[
    text("field", "Field").req().help("The relaxation field kept by Relax Nodes."),
    references("bodies", "Bodies", C::Body).help("Only the nodes of these bodies.").hint("the field's bodies"),
    range("polish", "Polish", D::None, 0.0, 1.0).help("Strength of one guarded smoothing sweep after the field is applied; 0 for none.").def("0.3"),
];

pub const RENUMBER: &[Field] = &[
    choice("what", "Renumber", &[("both", "Nodes and elements"), ("nodes", "Nodes"), ("elements", "Elements")])
        .help("Which labels change.")
        .def("\"both\""),
    choice("by", "Order", &[("compact", "Compact (label order)"), ("bandwidth", "Bandwidth (reverse Cuthill–McKee)"), ("body", "Body by body")])
        .help("The order of the new labels.")
        .def("\"compact\""),
    int("start", "Start", 1, 99_999_999).help("The first new label.").def("1"),
    int("offset", "Offset", -99_999_999, 99_999_999).help("Instead of an order and a start: add this to every label.").hint("none"),
];

pub const SEED_EDGE: &[Field] = &[
    json("edges", "Edges")
        .help("The edges to seed, as JSON; easiest with Use selection.")
        .req()
        .hint("[{\"body\", \"edge\"}, …] — Use selection with Select: Edges"),
    int("count", "Count", 1, 100_000).help("Number of elements along each edge.").hint("or an element size"),
    positive("size", "Element size", D::Length).help("Element length along the edge, instead of a count.").hint("or a count"),
    group(
        "bias",
        "Bias",
        &[
            range("ratio", "Ratio", D::None, 1.0, 1e4).help("Largest element over smallest; 1 is uniform.").req().def("1"),
            choice("toward", "Small toward", &[("start", "Start"), ("end", "End"), ("both", "Both ends")])
                .help("Which end gets the small elements.")
                .req()
                .def("\"end\""),
        ],
    )
    .help("Grade the elements along the edge."),
    choice("constraint", "Constraint", &[("exact", "Exact"), ("min", "At least"), ("free", "Free (a goal)")])
        .help("Exact: the mesher must match the count. At least: it may add more. Free: a goal only.")
        .def("\"exact\""),
];

pub const SEED_CLEAR: &[Field] = &[json("edges", "Edges").help("The edges to clear, as JSON; every seed when empty.").hint("every seed")];

pub const SIZE_SET: &[Field] = &[
    NAME,
    positive("size", "Element size", D::Length).help("The element size in the region.").req(),
    face_rows("faces", "Faces").help("Faces to refine. Pick them, then Use as faces.").hint("or bodies, a box or a sphere"),
    references("bodies", "Bodies", C::Body).help("Whole bodies meshed at this size.").hint("or faces, a box or a sphere"),
    group(
        "box",
        "Box",
        &[
            xyz("min", "Box min", D::Length).help("The box corner with the smallest coordinates.").req(),
            xyz("max", "Box max", D::Length).help("The box corner with the largest coordinates.").req(),
        ],
    )
    .help("An axis-aligned box to refine."),
    group(
        "sphere",
        "Sphere",
        &[
            xyz("center", "Centre", D::Length).help("The sphere's centre; click the model to place it.").req(),
            positive("radius", "Radius", D::Length).help("The sphere's radius.").req(),
        ],
    )
    .help("A sphere to refine."),
];

pub const SIZE_GROWTH: &[Field] =
    &[range("rate", "Growth rate", D::None, 0.05, 2.0).help("The largest change of element size per unit length away from a local size.").def("0.3")];

pub const SIZE_CLEAR: &[Field] =
    &[references("names", "Local sizes", C::SizeControl).help("The size controls to remove; every one when empty.").hint("every one")];

pub const DELAUNAY_POINTS: &[Field] = &[points("points", "Points", D::Length, 20_000).help("The points to triangulate.").req()];

pub const TAG: &[Field] = &[
    group(
        "geometry",
        "Geometry",
        &[
            reference("body", "Body", C::Body).help("A body of the open design.").req(),
            text("face", "Face").help("A persistent face name."),
            text("edge", "Edge").help("A persistent edge name."),
            text("vertex", "Vertex").help("A persistent vertex name."),
        ],
    )
    .help("The geometry the mesh entities are tagged with.")
    .req(),
    labels("elements", "Elements", LabelOf::Element).help("Element labels."),
    sides("faces", "Element faces", Side::Face).help("Element faces as element:face: 1-20:3, 25:1.").hint("element:face"),
    sides("edges", "Element edges", Side::Edge).help("Element edges as element:edge: 1-20:3, 25:1.").hint("element:edge"),
    labels("nodes", "Nodes", LabelOf::Node).help("Node labels."),
];

// ---- Geometry and files ----

pub const OPEN_DESIGN: &[Field] = &[
    path("path", "Design", false).help("The SolidCraft design file to open.").req(),
    choice("model", "Model", &[("keep", "Keep the current model"), ("new", "Start a new model")])
        .help("Keep the current model with the new geometry, or start a new one.")
        .def("\"keep\""),
];
const POINT: Dim = D::Length;
pub const DATUM_PLANE: &[Field] = &[
    NAME,
    variants(
        "by",
        "Definition",
        &[
            variant(
                "offset",
                "Principal plane",
                &[
                    choice("plane", "Plane", &[("yz", "YZ (normal X)"), ("zx", "ZX (normal Y)"), ("xy", "XY (normal Z)")])
                        .help("The principal plane.")
                        .req(),
                    num("offset", "Offset", POINT).help("Distance along the plane's normal.").def("0"),
                ],
            ),
            variant(
                "three_points",
                "Three points",
                &[
                    xyz("a", "Point A", POINT).help("First point.").req(),
                    xyz("b", "Point B", POINT).help("Second point.").req(),
                    xyz("c", "Point C", POINT).help("Third point.").req(),
                ],
            ),
            variant(
                "point_normal",
                "Point and normal",
                &[
                    xyz("point", "Point", POINT).help("A point on the plane.").req(),
                    xyz("normal", "Normal", D::None).help("The plane's normal direction.").req(),
                ],
            ),
            variant(
                "datum_points",
                "Three datum points",
                &[
                    reference("a", "Point A", C::DatumPoint).help("First datum point.").req(),
                    reference("b", "Point B", C::DatumPoint).help("Second datum point.").req(),
                    reference("c", "Point C", C::DatumPoint).help("Third datum point.").req(),
                ],
            ),
            variant(
                "normal_to_axis",
                "Normal to an axis",
                &[
                    reference("axis", "Axis", C::DatumAxis).help("The datum axis the plane is normal to.").req(),
                    reference("point", "Through", C::DatumPoint)
                        .help("A datum point on the plane; the axis's own point when empty.")
                        .hint("the axis's point"),
                ],
            ),
        ],
    )
    .help("How the plane is defined.")
    .req()
    .def(r#"{"type": "offset"}"#),
];
pub const DATUM_AXIS: &[Field] = &[
    NAME,
    variants(
        "by",
        "Definition",
        &[
            variant(
                "principal",
                "Principal axis",
                &[
                    choice("axis", "Axis", &[("x", "X"), ("y", "Y"), ("z", "Z")]).help("The principal direction.").req(),
                    xyz("through", "Through", POINT).help("A point the axis passes through; the origin when empty.").hint("the origin"),
                ],
            ),
            variant(
                "two_points",
                "Two points",
                &[
                    xyz("a", "Point A", POINT).help("First point on the axis.").req(),
                    xyz("b", "Point B", POINT).help("Second point on the axis.").req(),
                ],
            ),
            variant(
                "point_direction",
                "Point and direction",
                &[
                    xyz("point", "Point", POINT).help("A point on the axis.").req(),
                    xyz("direction", "Direction", D::None).help("The axis direction.").req(),
                ],
            ),
            variant(
                "datum_points",
                "Two datum points",
                &[
                    reference("a", "Point A", C::DatumPoint).help("First datum point on the axis.").req(),
                    reference("b", "Point B", C::DatumPoint).help("Second datum point on the axis.").req(),
                ],
            ),
        ],
    )
    .help("How the axis is defined.")
    .req()
    .def(r#"{"type": "principal"}"#),
];
pub const RIGID_WALL: &[Field] = &[
    reference("step", "Step", C::Step)
        .help("The explicit dynamics step the wall acts in; the model's only one when empty.")
        .hint("the explicit step"),
    reference("set", "Set", C::Set).help("The nodes the wall stops: pick them and Use selection.").req(),
    xyz("point", "Point", D::Length).help("A point on the wall.").req().def("[0, 0, 0]"),
    xyz("normal", "Normal", D::None).help("The wall's normal, pointing to the side the nodes stay on.").req().def("[0, 0, 1]"),
    FRICTION,
];
pub const RIGID_WALL_DELETE: &[Field] = &[
    reference("step", "Step", C::Step).help("The explicit dynamics step.").req(),
    int("index", "Wall", 0, MAX_ROWS as i64).help("The wall's place in the step's rigid walls, from 0.").req(),
];
pub const CENTRE: &[Field] = &[
    rows(
        "geometry",
        "Geometry",
        &[
            reference("body", "Body", C::Body).help("A body of the open design.").req(),
            text("face", "Face").help("A persistent face name of the body."),
            text("edge", "Edge").help("A persistent edge name of the body."),
            text("vertex", "Vertex").help("A persistent vertex name of the body."),
        ],
    )
    .help("Faces or edges (a hole's faces or rim give its centre), or one vertex."),
    labels("nodes", "Nodes", LabelOf::Node).help("Mesh node labels: their mean.").hint("labels: 1, 2, 10-20"),
];
pub const DATUM_POINT: &[Field] = &[
    NAME,
    variants(
        "by",
        "Definition",
        &[
            variant("coordinates", "Coordinates", &[xyz("point", "Point", POINT).help("The point's coordinates.").req()]),
            variant(
                "vertex",
                "Vertex",
                &[
                    reference("body", "Body", C::Body).help("The body the vertex belongs to.").req(),
                    text("vertex", "Vertex")
                        .help("The vertex's persistent name: pick it with Select: Vertices and Use selection (geometry.bodies lists every body's vertices).")
                        .req(),
                ],
            ),
            variant(
                "edge",
                "On an edge",
                &[
                    reference("body", "Body", C::Body).help("The body the edge belongs to.").req(),
                    text("edge", "Edge").help("The edge's persistent name.").req(),
                    range("fraction", "Fraction", D::None, 0.0, 1.0).help("How far along the edge from its start: 0.5 is its midpoint.").def("0.5"),
                ],
            ),
            variant(
                "face_centre",
                "Face centre",
                &[
                    reference("body", "Body", C::Body).help("The body the face belongs to.").req(),
                    text("face", "Face").help("The face's persistent name.").req(),
                ],
            ),
            variant(
                "axis_plane",
                "Axis and plane",
                &[
                    reference("axis", "Axis", C::DatumAxis).help("The datum axis.").req(),
                    reference("plane", "Plane", C::DatumPlane).help("The datum plane the axis crosses.").req(),
                ],
            ),
        ],
    )
    .help("How the point is defined.")
    .req()
    .def(r#"{"type": "coordinates"}"#),
];
pub const VIEW_SAVE: &[Field] = &[
    text("name", "Name").help("The view's name; a view of that name is replaced.").req(),
    json("camera", "Camera")
        .help("ui.camera's parameters: yaw and pitch in degrees, distance and target in the display units; the camera last aimed when empty.")
        .hint(r#"{"yaw": -70, "pitch": 18, "fit": true}"#),
    flag("results", "Results display").help("Keep what the results display shows with the camera.").def("true"),
];
pub const SOLVE_SETTINGS: &[Field] = &[
    STEP,
    choice("backend", "Backend", &[("gpu", "GPU"), ("cpu", "CPU"), ("none", "No preference")])
        .help("The backend the step prefers to be solved on, saved with the model; this machine may override it.")
        .req(),
    choice("fallback", "Fallback", &[("cpu", "CPU"), ("none", "None (refuse)")])
        .help("Where a GPU preference falls back on a machine with no suitable GPU.")
        .hint("CPU"),
];
pub const VIEW_NAME: &[Field] = &[text("name", "View").help("The saved view.").req()];
pub const RENDER_SAVE: &[Field] = &[
    text("name", "Name").help("The render's name; a render of that name is replaced.").req(),
    text("view", "View").help("The saved view the clip shows (view.save keeps one).").req(),
    text("step", "Step").help("The step whose frames it plays.").hint("the view's step"),
    json("frames", "Frames").help("\"all\", or [from, to]: the frames drawn, from 0, inclusive.").hint("\"all\""),
    int("fps", "Frames per second", 1, 120).help("The clip's frame rate.").def("12"),
    group(
        "hold",
        "Hold",
        &[
            int("each", "Each picture", 1, 10_000).help("Clip frames each picture shows.").def("1"),
            int("last", "Last picture", 1, 10_000).help("Clip frames the last picture shows: the settled state held.").def("1"),
        ],
    )
    .help("How long each picture shows, in frames of the clip."),
    tuple("size", "Size", &["Width", "Height"], 2, true)
        .help("The picture's width and height in pixels, 16 to 4096; the caption bar adds its height under it.")
        .hint("[1280, 720]"),
    text("caption", "Caption").help("A line of text under the picture."),
    flag("legend", "Legend").help("The contour legend under the picture.").def("true"),
    flag("time_label", "Time").help("The frame's time under the picture.").def("true"),
    choice("format", "Format", &[("png-sequence", "PNG sequence (encode with ffmpeg)")])
        .help("Numbered PNG pictures and the ffmpeg command line that makes a video of them; FieldCraft bundles no video encoder.")
        .def("\"png-sequence\""),
];
pub const RENDER_NAME: &[Field] = &[text("name", "Render").help("The saved render.").req()];
pub const EXPORT_ANIMATION: &[Field] = &[
    text("render", "Render").help("The saved render to draw (results.render.save keeps one).").req(),
    path("path", "File", true)
        .help("The clip to make: the pictures go beside it (clip_00000.png, …), with the ffmpeg command line that encodes them as it.")
        .req(),
];
pub const VIEW_ITEMS: &[Field] = &[rows(
    "items",
    "Items",
    &[
        choice(
            "kind",
            "Kind",
            &[
                ("body", "Body"),
                ("surface", "Surface"),
                ("datum", "Datum"),
                ("boundary_condition", "Boundary condition"),
                ("load", "Load"),
                ("reference_point", "Reference point"),
                ("coupling", "Coupling"),
            ],
        )
        .help("What kind of item it is.")
        .req(),
        text("name", "Name").help("The item's name.").req(),
    ],
)
.help("The viewport items.")
.req()];
pub const DATUM_DELETE: &[Field] = &[
    reference("name", "Datum", C::Datum).req().help("The Datum to act on, by name."),
    flag("cascade", "Delete its partitions").help("Also delete the partitions that cut by it; without it, a datum in use is refused."),
];
pub const PARTITION: &[Field] = &[
    NAME,
    reference("datum", "Datum plane", C::DatumPlane).help("The plane that divides the bodies.").req(),
    references("bodies", "Bodies", C::Body).help("The bodies to divide; every body the plane cuts when empty.").hint("every body the plane divides"),
];
pub const PARTITION_DELETE: &[Field] = existing!("Partition", C::Partition);
pub const FACES_MERGE: &[Field] = &[
    NAME,
    face_rows("faces", "Faces")
        .help("Two or more adjacent faces of one body, meshed as one face. Shift-click faces in the view, then Use selection.")
        .req()
        .hint("Shift-click faces"),
];
pub const FACE_ABSORB: &[Field] = &[
    text("name", "Name").help("The absorption's name (with From check: the prefix of each one's, Name-1, Name-2, …).").hint("Absorb"),
    group("face", "Face", FACE_ROW).help("The tiny or sliver face to absorb; leave empty with From check."),
    text("into", "Into").help("The neighbour face that absorbs it, by persistent name.").hint("the neighbour sharing the longest edges"),
    flag("from_check", "From check").help("Absorb every small and sliver face the last Check Geometry found, each into its best neighbour."),
];
pub const STITCH: &[Field] = &[
    NAME,
    text("surface", "Surface").help("The imported surface to stitch, by name.").req(),
    positive("tolerance", "Tolerance", D::Length)
        .help("Gaps narrower than this close: open points this close become one, and a point this close to another open edge splits it.")
        .req(),
];
pub const STITCH_DELETE: &[Field] = existing!("Stitch", C::Stitch);
pub const FACE_SPLIT: &[Field] = &[
    NAME,
    group("face", "Face", FACE_ROW).help("The face to split.").req(),
    reference("datum", "Datum plane", C::DatumPlane).help("The plane the face is split along.").req(),
];
pub const FACE_SPLIT_DELETE: &[Field] = existing!("Face split", C::FaceSplit);
pub const FACES_MERGE_DELETE: &[Field] = existing!("Face merge", C::FaceMerge);
pub const DEFEATURE: &[Field] = &[
    NAME,
    text("features", "Features")
        .help("The feature ids to take out for meshing (3 or F3), as Propose Defeaturing and the feature graph list them; separate them with commas.")
        .req(),
];
pub const DEFEATURE_DELETE: &[Field] = existing!("Defeature", C::Defeature);
pub const DEFEATURE_PROPOSE: &[Field] = &[
    reference("body", "Body", C::Body).help("The body whose holes, fillets, chamfers and bosses to list.").req(),
    group(
        "below",
        "Propose below",
        &[
            positive("hole_diameter", "Hole diameter", D::Length).help("Propose holes narrower than this."),
            positive("fillet_radius", "Fillet radius", D::Length).help("Propose fillets with a smaller radius."),
            positive("chamfer", "Chamfer", D::Length).help("Propose chamfers smaller than this."),
            positive("boss_height", "Boss height", D::Length).help("Propose bosses lower than this."),
        ],
    )
    .help("Size thresholds: a feature below its threshold is proposed."),
    flag("from_check", "From check").help("Also propose the features with a face the last Check Geometry found."),
    text("name", "Defeature").help("The defeature the list's on/off switches edit.").hint("Defeature"),
];
pub const MESH_PLAN: &[Field] = &[
    reference("body", "Body", C::Body).help("The body of the design to plan (as opened, not a partition cell).").req(),
    text("name", "Name").help("The name the plan is saved under.").hint("Plan-<body>"),
    positive("size", "Element size", D::Length).help("The largest element size; thin walls make it smaller.").hint("1/25 of the diagonal"),
    int("layers", "Layers through walls", 1, 10).help("Elements through the thinnest wall the features give (a shell, a plate).").hint("2"),
    int("max_elements", "Most elements", 1, 5_000_000).help("The size grows past the walls' when more elements would be needed.").hint("300000"),
];
pub const MESH_PLAN_APPLY: &[Field] = &[
    text("plan", "Plan").numbers(crate::units::MESH_PLAN).help("The name of a saved mesh plan (mesh.plan), or a plan object.").req(),
    flag("symmetry", "Use symmetry")
        .help("Mesh one cell of a body that is its own mirror image across planes through its middle, and mirror it.")
        .def("true"),
    flag("corners", "Copy corners")
        .help("With symmetry: when the body's two ends along another axis are alike (an enclosure's corners), mesh one end and mirror it onto the other, the middle meshed to fit both.")
        .def("true"),
    text("name", "Save as").help("The name a plan given as an object is saved under.").hint("Plan-<body>"),
];
pub const FEATURES: &[Field] = &[reference("body", "Body", C::Body).help("Only this body's features.").hint("every body")];
pub const GEOMETRY_CHECK: &[Field] = &[
    reference("body", "Body", C::Body).help("The body to check; every body when empty.").hint("every body"),
    positive("size", "Element size", D::Length)
        .help("The element size you intend to mesh with; features are small against it.")
        .hint("1/20 of the diagonal"),
    range("short_edge", "Short edge", D::None, 0.0, 10.0).help("Edges shorter than this fraction of the size are reported.").def("0.3"),
    range("small_face", "Small face", D::None, 0.0, 10.0).help("Faces with area below this fraction of size² are reported.").def("0.1"),
    range("sliver", "Sliver width", D::None, 0.0, 10.0).help("Faces narrower than this fraction of the size are reported.").def("0.3"),
    range("gap", "Gap", D::None, 0.0, 10.0).help("Parallel faces closer than this fraction of the size, but not touching, are reported.").def("0.3"),
    range("parallel_angle", "Parallel within", D::Degree, 0.0, 60.0).help("Faces this close to parallel count for the gap.").def("10"),
    flag("mark", "Mark findings").help("Mark what is found in the viewport until the design changes.").def("true"),
];
pub const TESSELLATE: &[Field] = &[reference("body", "Body", C::Body).help("The body to tessellate.").req(), TOLERANCE];
pub const IMPORT_MESH: &[Field] = &[path("path", "3MF file", false).help("A 3MF surface mesh with named faces.").req()];
pub const IMPORT_ORPHAN: &[Field] = &[
    path("path", "MSH file", false)
        .help("ASCII Gmsh MSH 4.1. Unitless coordinates use the current length unit. Replaces the mesh; remove existing mesh selections first, or import into a new model.")
        .req(),
    choice("unmatched", "Unmatched surfaces", &[("refuse", "Refuse"), ("shell", "Keep intentional shells")])
        .help("Boundary surfaces select existing volume faces. Keep shells only for intentional surface elements that match no volume face.")
        .def("\"refuse\""),
];
pub const FILE_OPEN: &[Field] = &[path("path", "Model", false).help("The FieldCraft model file to open.").req()];
const SAVE_RESULTS: Field = flag("results", "Keep results")
    .help("Keep every solved step's results in the file, so it opens showing them without solving again; off for a smaller file.")
    .def("true");
const SAVE_GEOMETRY: Field = choice("geometry", "Geometry", &[("refer", "Refer to the files"), ("embed", "Embed a copy")])
    .help("Refer to the design and meshes by path, or embed a copy of them in the file so it opens anywhere: whoever gets the file then gets your geometry too.")
    .hint("as the model's file");
pub const FILE_SAVE: &[Field] =
    &[path("path", "File", true).help("Where to save; the model's own file when empty.").hint("the model's file"), SAVE_RESULTS, SAVE_GEOMETRY];
pub const FILE_SAVE_AS: &[Field] = &[path("path", "File", true).help("The file to save to.").req(), SAVE_RESULTS, SAVE_GEOMETRY];
pub const RECOVERY_LIST: &[Field] =
    &[path("dir", "Folder", false).help("The recovery folder; the default one when empty.").hint("the recovery folder")];
pub const CRASH_REPORTS: &[Field] = &[
    int("limit", "Limit", 1, 1000).help("At most this many reports, newest first.").hint("20"),
    path("dir", "Folder", false).help("The crash folder; the default one when empty.").hint("the crash folder"),
];
pub const CRASH_REPORT: &[Field] = &[
    text("id", "Report").help("The crash report to show.").req(),
    path("dir", "Folder", false).help("The crash folder; the default one when empty.").hint("the crash folder"),
];
pub const CRASH_REPORTS_DELETE: &[Field] = &[
    text("id", "Report").help("The crash report to delete.").hint("or All"),
    flag("all", "All").help("Delete every crash report instead of one."),
    path("dir", "Folder", false).help("The crash folder; the default one when empty.").hint("the crash folder"),
];
pub const BUG_REPORT: &[Field] = &[
    text("title", "Title").help("The issue's title; one is made up when empty."),
    text("observed", "What happened").help("What went wrong; the last error when empty."),
    text("crash", "Crash report").help("A crash report to report (its id)."),
    flag("share_session", "Share my session").help("Attach the commands that built the model (they can name your files).").def("false"),
    flag("include_model", "Include the model").help("Attach the model file too (with Share my session).").def("false"),
];
pub const CRASH_TEST: &[Field] =
    &[choice("mode", "Mode", &[("panic", "Panic"), ("recovered", "Recovered panic"), ("abort", "Abort")]).help("How to crash.").req()];
pub const RECOVER: &[Field] = &[
    text("id", "Autosave").help("The autosave to recover.").req(),
    path("dir", "Folder", false).help("The recovery folder; the default one when empty.").hint("the recovery folder"),
];
pub const RECOVERY_DISCARD: &[Field] = &[
    text("id", "Autosave").help("The autosave to discard."),
    flag("all", "All").help("Discard every autosave."),
    path("dir", "Folder", false).help("The recovery folder; the default one when empty.").hint("the recovery folder"),
];
pub const SAMPLE: &[Field] = &[choice(
    "name",
    "Sample",
    &[
        ("cantilever", "Cantilever"),
        ("impact", "Impact"),
        ("parts", "Parts"),
        ("ibeam", "I-beam cantilever"),
        ("drop-test", "Drop test"),
        ("rubber-seal", "Rubber seal"),
        ("flow-past", "Flow past a bracket"),
        ("heated-flow", "Heated flow past a bracket"),
        ("cooling-heat-sink", "Heat sink in forced air"),
        ("cooling-enclosure", "Sealed enclosure cooling"),
        ("cooling-aio", "AIO cooler radiator"),
        ("heated-bracket", "Heated bracket"),
        ("bracket-hex", "Hex-meshed bracket"),
        ("crank-arm", "Bicycle crank arm"),
        ("crank-mirror", "Mirrored crank arm"),
    ],
)
.help("The sample model to open.")
.req()
.def("\"cantilever\"")];
pub const COUNT: &[Field] = &[int("count", "Count", 1, 1_000_000).help("How many steps.").def("1")];

pub const DECK: &[Field] = &[
    choice("format", "Format", &[("abaqus", "Abaqus input (.inp)"), ("nastran", "Nastran bulk data (.bdf)"), ("lsdyna", "LS-DYNA keyword (.k)")])
        .help("The solver deck format.")
        .req()
        .def("\"abaqus\""),
    path("path", "File", true).help("Where to write the deck.").req(),
    references("steps", "Steps", C::Step).help("The steps to write; all when empty.").hint("all"),
    choice(
        "units",
        "Units",
        &[
            ("mm-t-s-N-MPa", "mm, t, s, N, MPa"),
            ("m-kg-s-N-Pa", "m, kg, s, N, Pa"),
            ("mm-kg-ms-kN-GPa", "mm, kg, ms, kN, GPa"),
            ("mm-g-ms-N-MPa", "mm, g, ms, N, MPa"),
            ("in-lbf-s-psi", "in, lbf·s²/in, s, lbf, psi"),
        ],
    )
    .help("The consistent unit system the deck is written in, stated in its header.")
    .hint("the display system's"),
];

pub const EXPORT_STL: &[Field] = &[
    path("path", "File", true).help("Where to write the STL file.").req(),
    reference("set", "Set", C::Set).help("Write only this set's element faces; the whole outer surface when empty.").hint("the outer surface"),
    choice("format", "Format", &[("binary", "Binary"), ("ascii", "ASCII")]).help("Binary is compact; ASCII is readable text.").def("\"binary\""),
];

// ---- Solving ----

const STEP: Field = reference("step", "Step", C::Step).help("The step to solve; the only one when empty.").hint("the only one");
const SOLVED: Field = reference("step", "Step", C::SolvedStep).help("The solved step; the step shown when empty.").hint("the step shown");
const MAX_DOFS: Field = int("max_dofs", "Max DOFs", 1, 6_000_000).help("Cap on the degrees of freedom solved.").hint("6000000");
const MAX_FACTOR: Field = int("max_factor_entries", "Max factor entries", 1, super::solve::FACTOR_ENTRY_CEILING as i64)
    .help("Cap on factor entries (8 bytes each): 400M by default; raise deliberately up to 2B (16 GB of values alone). Indices, workspace and other solver storage need additional memory.")
    .hint("400000000");
const THREADS: Field = int("threads", "Threads", 1, 256).help("Threads to use; the session's setting when empty.").hint("the session's");

pub const SOLVE_STATIC: &[Field] = &[STEP, MAX_DOFS, MAX_FACTOR, THREADS];
pub const SOLVE_NONLINEAR: &[Field] = &[
    STEP,
    MAX_DOFS,
    MAX_FACTOR,
    int("max_newton_iterations", "Max Newton iterations", 1, 100_000).help("Cap on the Newton iterations over the whole step.").hint("100000"),
    THREADS,
];
pub const SOLVE_FREQUENCY: &[Field] = &[
    STEP,
    choice("mass", "Mass matrix", &[("consistent", "Consistent"), ("lumped", "Lumped")])
        .help("Consistent is more accurate; lumped is cheaper.")
        .def("\"consistent\""),
    num("shift", "Shift", D::Eigenvalue).help("Look for frequencies near this eigenvalue (ω²); 0 finds the lowest.").hint("0"),
    MAX_DOFS,
    MAX_FACTOR,
    int("max_iterations", "Max Lanczos vectors", 1, 3000).help("Cap on the Lanczos vectors.").hint("3000"),
    THREADS,
];
pub const SOLVE_BUCKLE: &[Field] =
    &[STEP, MAX_DOFS, MAX_FACTOR, int("max_iterations", "Max Lanczos vectors", 1, 3000).help("Cap on the Lanczos vectors.").hint("3000"), THREADS];
pub const SOLVE_EXPLICIT: &[Field] = &[
    STEP,
    MAX_DOFS,
    int("max_increments", "Max increments", 1, 1_000_000_000).help("Cap on the time increments.").hint("1000000000"),
    int("max_output_values", "Max output values", 1, 50_000_000).help("Cap on the result values kept.").hint("50000000"),
    int("max_contact_faces", "Max contact faces", 1, 10_000_000).help("Cap on the faces in contact search.").hint("10000000"),
    THREADS,
];
pub const SOLVE_HISTORY: &[Field] = &[
    SOLVED.req(),
    int("every", "Keep one in", 1, 1_000_000_000).help("Keep one history sample in this many."),
    references("sets", "Sets", C::Set)
        .help("Nonlinear static steps: the sets whose mean displacement and reaction sum to give per frame; the boundary-condition sets when none."),
    int("node", "Node", 1, 4_294_967_295).help("Nonlinear static steps: a node whose displacement to give per frame."),
];
pub const SOLVE_FLUID: &[Field] = &[
    STEP,
    choice("backend", "Backend", &[("cpu", "CPU"), ("gpu", "GPU")])
        .help("The CPU (double precision, the same result on every machine) or the GPU the build and FIELDCRAFT_GPU pick (single precision).")
        .hint("the step's preference, else CPU"),
    int("max_cells", "Max cells", 1, fieldcraft_lbm::Limits::MAX_ADDRESSABLE_CELLS as i64)
        .help("Cap on the lattice cells: 16M by default; may be raised up to what this machine's memory holds (the estimate below).")
        .hint("16000000"),
    int("max_steps", "Max steps", 1, 10_000_000).help("Cap on the lattice time steps.").hint("10000000"),
    int("max_cell_updates", "Max cell updates", 1, i64::MAX)
        .help("Cap on cells × steps: 1e12 by default; may be raised up to what this machine's memory holds.")
        .hint("1000000000000"),
    int("max_frame_bytes", "Max frame bytes", 1, i64::MAX)
        .help("Cap on the memory for result frames (about 16 bytes per cell per frame, 20 with heat): 2 GB by default; may be raised up to this machine's bound (the estimate below).")
        .hint("2000000000"),
    THREADS,
    path("checkpoint", "CPU restart file", true)
        .help("Optional native CPU checkpoint destination. Saves complete restart state, distinct from the display results in a model file; kept after success or cancellation.")
        .hint("disabled unless a path is supplied"),
    int("checkpoint_every", "Save every steps", 1, 10_000_000)
        .help("Required with a CPU restart file, including on resume: save every N lattice steps and at a sound final stop. Choose for the run size and storage cost; cadence may change on resume.")
        .hint("required with a restart file"),
    flag("resume", "Resume CPU restart")
        .help("Restore the explicit restart file using the same model, geometry and CPU settings. Thread count and checkpoint cadence may change. GPU and browser restart are not supported yet."),
    flag("estimate", "Estimate only").help("Return the run's memory estimate without solving."),
];
pub const SOLVE_PROBE: &[Field] = &[
    SOLVED.req(),
    labels("nodes", "Nodes", LabelOf::Node).help("Node labels to read.").hint("labels: 1, 2, 10-20"),
    reference("set", "Set", C::Set).help("A set whose nodes to read."),
    points("at", "Points", D::Length, 10_000).help("Points to read at the nearest nodes."),
    int("mode", "Mode", 1, 1_000_000).help("The mode to read, for frequency and buckling steps."),
    int("frame", "Frame", 0, 1_000_000_000).help("The frame to read, for explicit and transient steps."),
    flag("elements", "Element results").help("Also read element results (stresses) around the nodes."),
    points("points", "Fluid points", D::Length, 10_000).help("Points in the fluid to read."),
    json("forces", "Fluid walls").help("Walls whose fluid forces to sum, as JSON.").hint("[\"x_min\", \"F2:top\"]"),
    json("heat", "Fluid walls (heat)").help("Walls whose heat histories to read, as JSON.").hint("[\"y_min\", \"F2:top\"]"),
    group(
        "line",
        "Fluid line",
        &[
            xyz("from", "From", D::Length).help("The line's first end."),
            xyz("to", "To", D::Length).help("The line's last end."),
            int("count", "Count", 1, crate::flow::line::MAX_POINTS as i64)
                .help("Points evenly along the line, both ends included.")
                .hint("one per cell spacing"),
            choice(
                "field",
                "Field",
                &[
                    ("velocity_x", "Velocity X"),
                    ("velocity_y", "Velocity Y"),
                    ("velocity_z", "Velocity Z"),
                    ("speed", "Speed"),
                    ("pressure", "Pressure"),
                    ("temperature", "Temperature"),
                    ("wall_heat_flux", "Wall heat flux (into a box wall)"),
                ],
            )
            .help("The field to read along the line."),
            choice("face", "Wall", BOX_FACES).help("The box wall a wall-heat-flux line lies on."),
        ],
    )
    .help("Read a fluid field along a line, with its largest and smallest value and where they are."),
];

// ---- Results ----

pub const DEFORMED: &[Field] = &[
    SOLVED,
    flag("show", "Show").help("Show the deformed shape.").def("true"),
    number_or_auto("scale", "Scale", D::None, Some(0.0)).help("Displacement magnification: auto, or a factor (1 is true scale).").def("\"auto\""),
];
pub const CONTOUR: &[Field] = &[
    SOLVED,
    choice(
        "field",
        "Field",
        &[
            ("von_mises", "Von Mises stress"),
            ("displacement", "Displacement"),
            ("displacement_x", "Displacement X"),
            ("displacement_y", "Displacement Y"),
            ("displacement_z", "Displacement Z"),
            ("velocity", "Velocity magnitude"),
            ("pressure", "Pressure (fluid)"),
            ("temperature", "Temperature"),
            ("wall_heat_flux", "Wall heat flux (into wall)"),
            ("max_principal", "Max principal"),
            ("mid_principal", "Mid principal"),
            ("min_principal", "Min principal"),
            ("stress_xx", "Stress XX"),
            ("stress_yy", "Stress YY"),
            ("stress_zz", "Stress ZZ"),
            ("stress_xy", "Stress XY"),
            ("stress_xz", "Stress XZ"),
            ("stress_yz", "Stress YZ"),
            ("contact_pressure", "Contact pressure"),
            ("contact_gap", "Contact gap"),
            ("contact_status", "Contact status"),
            ("contact_shear", "Contact shear"),
            ("contact_slip", "Contact slip"),
            ("none", "None"),
        ],
    )
    .help("The result shown as colours.")
    .req()
    .def("\"von_mises\""),
    int("bands", "Bands", 2, 24).help("Number of colour bands.").def("12"),
    group(
        "range",
        "Legend range",
        &[
            num("min", "Min", D::None).help("Values below it take the lowest colour.").hint("the field's minimum"),
            num("max", "Max", D::None).help("Values above it are drawn in the over colour.").hint("the field's maximum"),
        ],
    )
    .help("Clamp the legend, e.g. below a singular peak; empty follows the field."),
];
pub const MODES: &[Field] = &[SOLVED, int("mode", "Mode", 1, 1_000_000).help("The mode shape shown.").def("1")];
pub const ANIMATE: &[Field] = &[flag("on", "Animate").help("Play the deformation or the frames.").def("true")];
pub const FRAME: &[Field] = &[
    SOLVED,
    int("frame", "Frame", 0, 1_000_000_000).help("The frame shown."),
    positive("time", "Time", D::None).help("Show the frame nearest this time (a nonlinear static step: its load factor); the last when empty. A time is in the display system's unit.").hint("the last frame"),
];
const FLUID_FRAME: Field = int("frame", "Frame", 0, 1_000_000_000).help("The frame shown, from 0.").hint("the frame shown");
pub const SLICE: &[Field] = &[
    SOLVED,
    choice("field", "Field", &[("velocity", "Velocity magnitude"), ("pressure", "Pressure"), ("temperature", "Temperature (heated flow)")])
        .help("The field coloured on the slice and the obstacles.")
        .hint("the field shown"),
    choice("axis", "Normal", &[("x", "X"), ("y", "Y"), ("z", "Z"), ("none", "None (hide)")])
        .help("The axis the slice is normal to.")
        .hint("the slice shown"),
    number_or_auto("position", "Position", D::Length, None)
        .help("Where the slice crosses its axis; auto: the domain's middle.")
        .hint("the position shown"),
    FLUID_FRAME,
];
pub const FLOW_SURFACE: &[Field] = &[
    SOLVED,
    flag("show", "Show").help("Colour the obstacles by the selected field.").def("true"),
    choice(
        "field",
        "Field",
        &[("velocity", "Velocity"), ("pressure", "Pressure"), ("temperature", "Temperature"), ("wall_heat_flux", "Wall heat flux (into wall)")],
    )
    .help("Colour the obstacle faces by the adjacent fluid field, or signed wall heat received per area.")
    .hint("the field shown"),
    FLUID_FRAME,
];
pub const STREAMLINES: &[Field] = &[
    SOLVED,
    choice("field", "Field", &[("velocity", "Velocity"), ("pressure", "Pressure (plain lines)"), ("temperature", "Temperature")])
        .help("Colour streamlines by speed or temperature; pressure leaves the lines in one colour.")
        .hint("the field shown"),
    flag("show", "Show").help("Draw the streamlines.").def("true"),
    int("count", "Count", 1, 2000).help("Seeds on the inlet face, or points along the seed line.").def("64"),
    int("max_steps", "Max steps", 1, 100_000).help("Most steps of half a cell each streamline takes.").def("4000"),
    points("seeds", "Seeds", D::Length, 2000).help("The points streamlines start from; the inlet face when empty."),
    group(
        "line",
        "Seed line",
        &[xyz("from", "From", D::Length).help("The seed line's first end."), xyz("to", "To", D::Length).help("The seed line's last end.")],
    )
    .help("Start `count` streamlines evenly along a line instead of the inlet face."),
    FLUID_FRAME,
];
pub const EXPORT_CSV: &[Field] = &[
    choice("what", "Export", &[("reactions", "Reaction sums"), ("probe", "Probed nodes")]).help("What to write.").req().def("\"reactions\""),
    path("path", "File", true).help("The CSV file to write.").req(),
    SOLVED,
    references("sets", "Sets", C::Set).help("Reactions: the sets to sum; the step's BC sets when empty.").hint("the step's BC sets"),
    xyz("about", "About point", D::Length).help("Reactions: the point moments are taken about."),
    labels("nodes", "Nodes", LabelOf::Node).help("Probe: the node labels."),
    reference("set", "Set", C::Set).help("Probe: a set whose nodes to read."),
    points("at", "Points", D::Length, 10_000).help("Probe: points to read at the nearest nodes."),
    int("mode", "Mode", 1, 1_000_000).help("Probe: the mode to read, for frequency and buckling steps."),
    int("frame", "Frame", 0, 1_000_000_000).help("Probe: the frame to read, for explicit and transient steps."),
];
pub const RESULTS_XY: &[Field] = &[
    SOLVED,
    json("y", "Quantity")
        .help("What to plot: a nodal field at a node, an explicit energy or history, a history node's quantity, or \"frequency\".")
        .hint(r#"{"field": "displacement_z", "node": 21}"#)
        .req(),
    choice("x", "Against", &[("auto", "The step's own axis"), ("load_factor", "Load factor"), ("time", "Time"), ("mode", "Mode")])
        .help("The x axis: the load factor (nonlinear static), time (explicit, heat) or mode number (frequency).")
        .def("\"auto\""),
    int("every", "Keep one in", 1, 1_000_000_000).help("Keep one point in this many; by default enough to fit 100000 points."),
    path("path", "CSV file", true).help("Also write the series as a two-column CSV."),
];
pub const RESULTS_PATH: &[Field] = &[
    SOLVED,
    text("field", "Nodal field").req().help("A results.contour field id, such as stress_xx or temperature.").def("\"displacement\""),
    points("points", "Polyline points", D::Length, 1000).req().help("Two or more points in undeformed reference coordinates."),
    int("samples", "Samples", 2, 10_000).help("Equally spaced arc-length stations along the complete polyline, including both ends.").def("101"),
    int("frame", "Frame", 0, 1_000_000_000).help("A stored field frame; empty means final state."),
    int("mode", "Mode", 1, 1_000_000).help("Frequency results only; empty means mode 1."),
    path("path", "CSV file", true).help("Optional distance,x,y,z,value CSV in display units; gaps have empty values."),
];
pub const EXPORT_VTU: &[Field] = &[
    path("path", "File", true).help("A .vtu file for one state, or a .pvd collection for every frame or mode.").req(),
    SOLVED,
    int("mode", "Mode", 1, 1_000_000).help("Frequency steps: the mode to write; every mode in a .pvd when empty.").hint("all (.pvd) or 1"),
    choice("frames", "Frames", &[("all", "Every frame"), ("last", "The last frame")])
        .help("A .pvd series: every frame or only the last.")
        .def("\"all\""),
    choice("encoding", "Encoding", &[("ascii", "Text"), ("binary", "Binary (base64)")])
        .help("Text is readable; binary is smaller and exact.")
        .def("\"ascii\""),
];
pub const INITIAL_FIELD_MAP: &[Field] = &[
    NAME,
    choice("quantity", "Quantity", &[("temperature", "Temperature"), ("stress", "Stress")])
        .help("Map a temperature (where thermal strains are measured from) or a residual stress.")
        .req(),
    text("step", "Source step")
        .help("The solved step to map from: a heat transfer or coupled step's temperatures, or a static step's stresses.")
        .req(),
    path("file", "Source model", false)
        .help("The .fieldcraft model holding the step and its saved results, on its own mesh; this model when empty.")
        .hint("this model"),
    int("frame", "Frame", 0, 1_000_000_000).help("Heat transfer steps: the frame to map; the last when empty.").hint("the last"),
    reference("set", "Target set", C::Set).help("The nodes to map onto; every node of a solid element when empty.").hint("every solid node"),
    positive("tolerance", "Tolerance", D::Length)
        .help("How far outside the source mesh a node may lie and still take the value at the closest point.")
        .hint("a quarter of the source's element size"),
];
pub const INITIAL_FIELD_DELETE: &[Field] = &[text("name", "Initial field").req().help("The initial field to remove, by name.")];

pub const REACTIONS: &[Field] = &[
    SOLVED,
    references("sets", "Sets", C::Set)
        .help("The sets whose reaction forces to sum; the step's boundary-condition sets when empty.")
        .hint("the step's BC sets"),
    xyz("about", "About point", D::Length).help("The point moments are taken about."),
    int("frame", "Frame", 0, 1_000_000_000)
        .help("Nonlinear static steps: the frame (converged increment) whose reactions to sum; the final state when empty."),
];

#[cfg(test)]
mod tests {
    use fieldcraft_model::DisplaySystem;
    use fieldcraft_model::units::Preset;
    use serde_json::{Value, json};

    use super::*;
    use crate::{EngineError, Session};

    /// The form's starting parameters (and those with the required names and references filled
    /// in from the model) are accepted, or refused only for a required field still empty; never
    /// for a value the form put there itself. Errors other than parameter errors (no design, no
    /// results) come after the parameters were read.
    #[test]
    fn form_defaults_are_accepted_or_need_only_required_fields() {
        let (mut checked, mut failures) = (0, Vec::new());
        for c in crate::commands() {
            let Some(form) = c.form else { continue };
            if c.id.starts_with("file.") && c.id != "file.save_as" || c.id == "edit.undo" || c.id == "edit.redo" {
                continue; // they replace or write the session's model; their forms are paths and counts
            }
            for fill_refs in [false, true] {
                let limits = crate::Limits {
                    max_steiner: 10,
                    max_tets: 100,
                    max_hex_elements: 100,
                    max_import_nodes: 10,
                    max_import_elements: 100,
                    max_import_members: 100,
                    threads: 1,
                    memory: 0,
                };
                let mut s = Session { limits, ..Session::default() };
                crate::cmd::define::tests::defined_bar(&mut s);
                s.execute("solve.run", &json!({})).unwrap();
                let mut p = prefill(form, &json!({}));
                if fill_refs {
                    fill_required(&s, form, &mut p);
                }
                let missing: Vec<String> = problems(form, &p).into_iter().filter(|p| p.kind == ProblemKind::Required).map(|p| p.path).collect();
                match s.execute(c.id, &p) {
                    Err(EngineError::Param(e)) if e.kind == ProblemKind::Required && missing.contains(&e.path) => {}
                    Err(EngineError::Param(e)) => {
                        failures.push(format!("{}: {p} refused: {}, but only {missing:?} were missing", c.id, e.for_agents()));
                    }
                    // State, not parameters: the meshers' default is the open design's bodies, the
                    // bar has no explicit step; activation and probes need one of several fields.
                    Err(EngineError::Params(e)) if STATE.iter().any(|m| e.contains(m)) => {}
                    Err(EngineError::Params(e)) if !missing.iter().any(|k| e.contains(k.rsplit('/').next().unwrap_or(k))) => {
                        failures.push(format!("{}: {p} refused for `{e}`, but only {missing:?} were missing", c.id));
                    }
                    Err(EngineError::Params(_)) => {}
                    _ => checked += 1,
                }
            }
        }
        assert!(failures.is_empty(), "{failures:#?}");
        assert!(checked > 20, "{checked}");
    }

    /// Every parameter the engine refuses is a field of the command's form: the problem's path
    /// leads to a field, so a form can show it beside that field, by its label. Each field of
    /// each form is given a value of the wrong kind in turn, over the form's starting values.
    #[test]
    fn every_param_error_names_a_field_of_its_form() {
        let hostile = [json!("x"), json!(-1), json!(0), json!(1e300), json!([]), json!({}), json!(null), json!([1, "x"])];
        let (mut seen, mut failures) = (0, Vec::new());
        for c in crate::commands() {
            let Some(form) = c.form else { continue };
            if c.id.starts_with("file.") || c.id.starts_with("edit.") || form.iter().any(|f| matches!(f.kind, Kind::Path { .. })) {
                continue; // they replace the session's model or write files
            }
            let base = prefill(form, &json!({}));
            let mut tries = vec![json!({}), base.clone()];
            for f in form.iter().filter(|f| !f.key.is_empty()) {
                for v in &hostile {
                    let mut p = base.clone();
                    if let Some(o) = p.as_object_mut() {
                        o.insert(f.key.into(), v.clone());
                    }
                    tries.push(p);
                }
            }
            for p in tries {
                let limits = crate::Limits {
                    max_steiner: 10,
                    max_tets: 100,
                    max_hex_elements: 100,
                    max_import_nodes: 10,
                    max_import_elements: 100,
                    max_import_members: 100,
                    threads: 1,
                    memory: 0,
                };
                let mut s = Session { limits, ..Session::default() };
                crate::cmd::define::tests::defined_bar(&mut s);
                if let Err(EngineError::Param(e)) = s.execute(c.id, &p) {
                    seen += 1;
                    if e.field.is_none() || field_at(form, &e.path).is_none() {
                        failures.push(format!("{}: {p} refused `{}`, not a field of the form ({})", c.id, e.path, e.for_agents()));
                    }
                }
            }
        }
        failures.sort();
        failures.dedup();
        assert!(failures.is_empty(), "{failures:#?}");
        assert!(seen > 50, "{seen}");
    }

    /// A required field left empty says what to do by its kind: a name, a path, else "required".
    #[test]
    fn required_fields_are_phrased_by_kind() {
        let form = |fields: &'static [Field]| problems(fields, &json!({})).iter().map(Problem::for_form).collect::<Vec<_>>();
        assert_eq!(form(RECOVER), ["Autosave: required"]);
        assert_eq!(form(FILE_SAVE_AS), ["File: enter a path"]);
        assert!(form(MATERIAL).contains(&"Name: enter a name".to_string()), "{:?}", form(MATERIAL));
    }

    /// The commands whose form is the JSON text alone (no schema yet). None: keep it that way.
    #[test]
    fn every_command_has_a_form() {
        let without: Vec<&str> = crate::commands().iter().filter(|c| c.form.is_none()).map(|c| c.id).collect();
        assert_eq!(without, Vec::<&str>::new());
        // A command that takes nothing has an empty form, and says so.
        for c in crate::commands() {
            assert_eq!(c.form.is_some_and(<[Field]>::is_empty), c.params.trim_start().starts_with("{}"), "{}", c.id);
        }
    }

    const STATE: &[&str] = &[
        "open a design first",
        "update checks are off",
        "`boundary_conditions` or `loads`",
        "give one of",
        "is not an explicit",
        "is neither an explicit dynamic nor a nonlinear static step",
        "has no explicit dynamics step",
        "has a result",
        "no crash report",
        "FIELDCRAFT_CRASH_TEST",
        "no relaxation field",
        "no saved view",
        "no render",
    ];

    /// Fill top-level required names (a fresh one) and references (the first existing item).
    fn fill_required(s: &Session, form: &[Field], p: &mut Value) {
        let Some(obj) = p.as_object_mut() else { return };
        for f in form {
            if !f.required || obj.contains_key(f.key) {
                continue;
            }
            match f.kind {
                Kind::Text => {
                    obj.insert(f.key.into(), json!("Fresh"));
                }
                Kind::Ref(c) => {
                    if let Some(n) = s.names(c).first() {
                        obj.insert(f.key.into(), json!(n));
                    }
                }
                _ => {}
            }
        }
    }

    /// Every field of every form says what it means, and every number with a dimension shows
    /// its unit, or its dimension when the model has no unit system.
    #[test]
    fn every_field_has_help_and_a_unit() {
        fn walk(id: &str, fields: &[Field], missing: &mut Vec<String>) {
            for f in fields {
                if f.help.trim().is_empty() {
                    missing.push(format!("{id}: {} ({}) has no help", f.label, f.key));
                }
                if f.dim().is_some_and(|d| d != D::None) {
                    for preset in Preset::ALL {
                        // Angles are radians: a unit of their own, not a ratio's.
                        if f.unit_label(&DisplaySystem::new(preset)).is_none_or(|u| u.is_empty()) {
                            missing.push(format!("{id}: {} shows no unit in {preset:?}", f.label));
                        }
                    }
                }
                match f.kind {
                    Kind::Group(g) | Kind::Rows { fields: g, .. } | Kind::Either(g) => walk(id, g, missing),
                    Kind::NameOr { name, fields, .. } => {
                        walk(id, std::slice::from_ref(name), missing);
                        walk(id, fields, missing);
                    }
                    Kind::Variants { variants, .. } => variants.iter().for_each(|v| walk(id, v.fields, missing)),
                    _ => {}
                }
            }
        }
        let mut missing = Vec::new();
        for c in crate::commands() {
            walk(c.id, c.form.unwrap_or_default(), &mut missing);
        }
        assert!(missing.is_empty(), "{missing:#?}");
        let stress = positive("s", "S", D::Stress);
        assert_eq!(stress.unit_label(&DisplaySystem::new(Preset::MKgS)).as_deref(), Some("Pa"));
        assert_eq!(stress.unit_label(&DisplaySystem::default()).as_deref(), Some("MPa"));
        assert!(stress.help_text(&DisplaySystem::default()).contains("Unit: MPa."));
    }

    /// Every default is JSON; every variant field's default names a variant.
    #[test]
    fn defaults_parse() {
        fn walk(fields: &[Field]) {
            for f in fields {
                if let Some(d) = f.default {
                    assert!(serde_json::from_str::<serde_json::Value>(d).is_ok(), "{}: {d}", f.key);
                }
                match f.kind {
                    Kind::Group(g) | Kind::Rows { fields: g, .. } | Kind::Either(g) => walk(g),
                    Kind::Variants { variants, .. } => variants.iter().for_each(|v| walk(v.fields)),
                    _ => {}
                }
            }
        }
        for c in crate::commands() {
            if let Some(form) = c.form {
                walk(form);
            }
        }
    }
}
