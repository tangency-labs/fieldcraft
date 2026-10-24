//! Abaqus input (`.inp`) writer.
//!
//! The deck is written from the public Abaqus keyword and element references, in our own words;
//! nothing here comes from a deck produced by the product.
//!
//! # Layout
//!
//! The deck is *flat*: nodes, elements, sets, materials, sections and couplings sit at model level,
//! with no parts, assembly or instances. FieldCraft has one mesh with model-wide labels and no
//! instance transforms, so a part/assembly wrapper would only prefix every set name with an
//! instance name (`Part-1-1.root`) and bring the 80-character limit on the joined name closer. When
//! the model gains several parts or instances the writer switches to the part/assembly form.
//!
//! In order: heading, nodes (mesh nodes, then reference points), elements grouped by element type,
//! sets (a node set, an element set and an element-face surface per model set, each written when
//! the set resolves to some), materials, sections, couplings, contact, the initial state's
//! boundary conditions, then the steps.
//!
//! # Contact
//!
//! Each contact is a `*SURFACE INTERACTION` with its friction (`*FRICTION`, or
//! `*FRICTION, EXPONENTIAL DECAY` with the static and kinetic coefficients and the decay
//! coefficient, the same law). A pair is `*CONTACT PAIR, TYPE=SURFACE TO SURFACE` with the
//! secondary surface first, and `ADJUST=0.0` when initially overclosed nodes are adjusted. General
//! contact is the model's one `*CONTACT` definition: `*CONTACT INCLUSIONS, ALL EXTERIOR` for every
//! exterior face, a surface with a blank partner for self contact, and a
//! `*CONTACT PROPERTY ASSIGNMENT` line per contact. FieldCraft's penalty scale and its refusal of
//! initial overclosure have no Abaqus counterpart and are noted in comments (Abaqus chooses its
//! own penalty and resolves overclosure in the first increment).
//!
//! # Ties
//!
//! Each tie is `*TIE, NAME=…, ADJUST=NO` with the secondary surface first and its position
//! tolerance when it has one: FieldCraft ties the nodes where they are (no adjustment).
//!
//! # Steps and the activity table
//!
//! The model lists, for every step, exactly the boundary conditions and loads active in it; Abaqus
//! carries them from step to step. The writer bridges the two:
//! - a boundary condition active in every step that holds its degrees of freedom at zero is written
//!   once, in the initial state (model data);
//! - a general (static) step whose active boundary conditions equal the previous general step's
//!   writes none; one that only adds writes the additions (`OP=MOD`); otherwise it writes
//!   `*BOUNDARY, OP=NEW` and restates every active one;
//! - loads in a general step are restated, keyword by keyword, with `OP=NEW` whenever the active set
//!   changes, so a load not listed is removed;
//! - a frequency step is a perturbation of the last general step: its boundary conditions are
//!   written against that base state and do not carry on to later steps; loads have no effect on
//!   natural frequencies and are not written in it;
//! - every step restates its output request (`*OUTPUT, FIELD, OP=NEW`);
//! - a nonlinear static step is `*STEP, NLGEOM=YES` (or `NO` without nonlinear geometry) with
//!   `INC=` its increment cap, `*STATIC` with its initial, minimum and maximum increments over a
//!   step period of 1, and `*CONTROLS, PARAMETERS=LINE SEARCH` when it searches; its iteration cap
//!   and tolerances are left to Abaqus's controls, whose residual norms are defined differently.
//!   Abaqus starts each general step from the end of the previous one, while every model step
//!   starts from the undeformed state; for elastic materials the end state depends on the final
//!   loads alone, so the two agree. A plastic material's state depends on the path, so a deck
//!   with nonlinear steps, a plastic material and more than one general step is refused, as is
//!   a nonlinear step's load path (an amplitude on every load and boundary condition, not
//!   written yet).
//!
//! # Loads
//!
//! - Forces, moments and tractions are written as concentrated loads from the model's resolved
//!   nodal values (the same distribution every writer and the solver use), summed per node and
//!   degree of freedom so loads sharing a node add up.
//! - Pressure is a surface load on the set's element-face surface, summed per set.
//! - Gravity and body forces are element loads on element sets (the section sets when gravity acts
//!   on the whole model); when two active loads of one kind share elements they are summed and
//!   written per element.
//!
//! # Numbers and labels
//!
//! Real numbers are written in the shortest form that reads back to the same `f64` (exponent marker
//! `E`). Lines stay within 256 characters: element connectivity puts at most 15 nodes after the
//! label on the first line and 16 per continuation line; sets list 16 labels per line. Labels are
//! the model's names (a letter, then letters, digits, `_` or `-`, at most 64 characters), which
//! Abaqus reads unquoted. Abaqus ignores case in names, so two names that differ only in case are
//! reported as [`Error::Unsupported`], as is everything else the format cannot express.
//!
//! # Element types
//!
//! | FieldCraft | Solid section | Shell section | Beam | Rod |
//! |---|---|---|---|---|
//! | tet4 | C3D4 (H) | | | |
//! | tet10 | C3D10 (H) | | | |
//! | wedge6 | C3D6 (H) | | | |
//! | pyramid5 | C3D5 (H) | | | |
//! | hex8 | C3D8, C3D8R, C3D8I (H) | | | |
//! | hex20 | C3D20, C3D20R (H) | | | |
//! | tri3, quad4 | | S3 / S3R, S4 / S4R | | |
//! | quad8 | | S8R (reduced only) | | |
//! | beam2, rod2 | | | B31 | T3D2 |
//! | rod3 | | | | T3D3 |
//!
//! Trusses in a 2D model are T2D2 and T2D3.
//!
//! Point masses, springs and dashpots are elements numbered after the mesh's, one per resolved
//! copy, in an element set per definition (`<name>-mass`, `<name>-inertia`, `<name>-spring`, made
//! unique against the model's sets): `MASS` with `*MASS` (one mass for every direction),
//! `ROTARYI` with `*ROTARY INERTIA` (I11, I22, I33, I12, I13, I23), and `SPRING1`/`DASHPOT1` to
//! the ground and `SPRING2`/`DASHPOT2` between two nodes along a global degree of freedom, or
//! `SPRINGA`/`DASHPOTA` along the line of the nodes, with `*SPRING`/`*DASHPOT`. A mass that
//! differs by direction and a spring along a skew fixed vector are reported (they need a local
//! orientation, not written yet).
//!
//! | FieldCraft | Plane stress | Plane strain | Axisymmetric |
//! |---|---|---|---|
//! | tri3 | CPS3 | CPE3 | CAX3 |
//! | tri6 | CPS6 | CPE6 | CAX6 |
//! | quad4 | CPS4, CPS4R, CPS4I | CPE4, CPE4R, CPE4I | CAX4, CAX4R, CAX4I |
//! | quad8 | CPS8, CPS8R | CPE8, CPE8R | CAX8, CAX8R |
//!
//! 2D models: Abaqus's axisymmetric elements, like the model, take r = x and the axis along Y,
//! and their concentrated loads are totals over the circumference; a 2D section's thickness is
//! the `*SOLID SECTION` data line. Surfaces of 2D elements are their edges, `S1`… in the edge
//! order of the model; body forces of axisymmetric elements are `BR`, `BZ`.
//!
//! `H` marks the hybrid variant. A four-node tetrahedron integrates at one point either way, so
//! reduced integration writes C3D4. Not expressible: reduced integration of ten-node tetrahedra and
//! six-node wedges, full integration of eight-node shells, six-node triangular shells, orthotropic
//! materials outside solid sections (shell and beam material axes need orientations), and
//! distributing couplings that leave a translation free.
//!
//! # Materials
//!
//! Isotropic elasticity is `*ELASTIC, TYPE=ISOTROPIC` (or, for a neo-Hookean material in a deck
//! with nonlinear steps, `*HYPERELASTIC, NEO HOOKE` with C₁₀ = μ/2 and D₁ = 2/κ, the same strain
//! energy, and likewise `MOONEY-RIVLIN` (C₁₀, C₀₁, D₁) and `YEOH` (C₁₀, C₂₀, C₃₀, D₁, the higher
//! volumetric terms zero); a linear static step is then refused, since Abaqus's general static step would apply
//! the hyperelastic law), engineering constants
//! `TYPE=ENGINEERING CONSTANTS`, density `*DENSITY`. Von Mises plasticity is `*PLASTIC`: a
//! hardening table as its (yield stress, plastic strain) lines; bilinear isotropic hardening as
//! two lines, the second at a plastic strain of [`BILINEAR_REACH`] (Abaqus holds the last yield
//! stress beyond a table's end, so the line is written out further than any analysis strains);
//! bilinear kinematic hardening as `HARDENING=KINEMATIC` with its yield stress at plastic strains
//! 0 and 1 (Abaqus's linear kinematic model extends that line); Johnson–Cook as
//! `HARDENING=JOHNSON COOK` (A, B, n, m, melting and transition temperature, the transition
//! temperature being the model's room temperature) with `*RATE DEPENDENT, TYPE=JOHNSON COOK`
//! (C, ε̇₀) when C is not zero, and adiabatic heating as `*SPECIFIC HEAT` and
//! `*INELASTIC HEAT FRACTION`. Damage and erosion are failure models of Abaqus/Explicit, and this
//! writer writes Abaqus/Standard procedures: they are reported, not written.
//!
//! Node order: the model stores elements in Abaqus order, so connectivity passes through; the
//! model's face numbers are Abaqus's S1–S6 for solids, and SPOS/SNEG for shell faces 1/2.

use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;

use fieldcraft_model::{
    BeamProfile, BoundaryCondition, ContactSurfaces, CouplingKind, Elastic, Element, ElementKind, Hardening, Hyperelastic, InitialPenetration,
    Integration, Load, Model, OutputVariable, Plastic, Procedure, Resolved, Section, SectionKind, Step, Target,
};

use fieldcraft_model::units::Preset;

use crate::{Error, Result, unsupported};

pub mod read;
pub use read::{Read, Unknown, read, read_in};

const FORMAT: &str = "Abaqus";

/// The plastic strain a bilinear isotropic hardening line is written out to (1000%).
pub const BILINEAR_REACH: f64 = 10.0;

/// Longest line Abaqus reads.
pub const MAX_LINE: usize = 256;

/// Labels per set data line.
const SET_LABELS_PER_LINE: usize = 16;

/// Nodes on an element's first data line (after its label), and on each continuation line.
const FIRST_LINE_NODES: usize = 15;
const CONTINUATION_NODES: usize = 16;

/// Write the model as an Abaqus input deck. Fails when the model does not resolve (its check
/// finds an error) or has something the format cannot express ([`Error::Unsupported`], every such
/// item named).
pub fn write(model: &Model, units: Preset) -> Result<String> {
    let expressed = crate::express(model, units)?;
    write_resolved(&expressed, &expressed.resolve()?, units)
}

/// [`write`] for a model already expressed in `units` ([`crate::express`]) and resolved.
pub fn write_resolved(model: &Model, r: &Resolved, units: Preset) -> Result<String> {
    crate::refuse_fluid(FORMAT, model)?;
    let by_type = elements_by_type(model, r)?;
    check_expressible(model)?;
    let discrete = discrete_groups(model, r)?;
    let mut d = Deck::default();
    d.heading(units);
    d.nodes(model)?;
    d.elements(&by_type);
    d.discrete_elements(&discrete);
    d.sets(model, r)?;
    d.materials(model)?;
    d.sections(model)?;
    d.discrete_properties(&discrete)?;
    d.couplings(model, r);
    d.contacts(model, r)?;
    d.ties(model)?;
    d.equations(model, r)?;
    d.steps(model, r)?;
    Ok(d.out)
}

/// The Abaqus element type of an element of `kind` in `section`, or why there is none.
pub fn element_type(kind: ElementKind, section: &Section) -> std::result::Result<String, String> {
    use ElementKind as K;
    use SectionKind as S;
    let f = section.formulation;
    let reduced = f.integration == Integration::Reduced;
    if (f.hybrid && section.kind.dimension() != 3) || (f.incompatible_modes && !section.kind.is_continuum()) {
        return Err(format!("hybrid or incompatible-mode {} sections", section.kind.name()));
    }
    if let Some(prefix) = match section.kind {
        S::PlaneStress { .. } => Some("CPS"),
        S::PlaneStrain { .. } => Some("CPE"),
        S::Axisymmetric {} => Some("CAX"),
        _ => None,
    } {
        let suffix = match (kind, reduced, f.incompatible_modes) {
            (K::Quad4, false, true) => "4I",
            (_, true, true) => return Err("incompatible modes with reduced integration".into()),
            (_, _, true) => return Err(format!("incompatible modes on {} elements (quad4 only)", kind.name())),
            (K::Tri3, ..) => "3",
            (K::Tri6, true, _) => return Err("reduced integration of tri6 elements (Abaqus has no reduced-integration variant)".into()),
            (K::Tri6, ..) => "6",
            (K::Quad4, r, _) => ["4", "4R"][usize::from(r)],
            (K::Quad8, r, _) => ["8", "8R"][usize::from(r)],
            _ => return Err(format!("{} elements in a {} section", kind.name(), section.kind.name())),
        };
        return Ok(format!("{prefix}{suffix}"));
    }
    let plain = match (&section.kind, kind, reduced, f.incompatible_modes) {
        (S::Solid {}, K::Hex8, false, true) => "C3D8I",
        (S::Solid {}, _, true, true) => return Err("incompatible modes with reduced integration".into()),
        (S::Solid {}, _, _, true) => return Err(format!("incompatible modes on {} elements (hex8 only)", kind.name())),
        (S::Solid {}, K::Tet4, ..) => "C3D4",
        (S::Solid {}, K::Tet10 | K::Wedge6 | K::Pyramid5, true, _) => {
            return Err(format!("reduced integration of {} elements (Abaqus has no reduced-integration variant)", kind.name()));
        }
        (S::Solid {}, K::Tet10, ..) => "C3D10",
        (S::Solid {}, K::Wedge6, ..) => "C3D6",
        (S::Solid {}, K::Pyramid5, ..) => "C3D5",
        (S::Solid {}, K::Hex8, r, _) => ["C3D8", "C3D8R"][usize::from(r)],
        (S::Solid {}, K::Hex20, r, _) => ["C3D20", "C3D20R"][usize::from(r)],
        (S::Shell { .. }, K::Tri3, r, _) => ["S3", "S3R"][usize::from(r)],
        (S::Shell { .. }, K::Quad4, r, _) => ["S4", "S4R"][usize::from(r)],
        (S::Shell { .. }, K::Quad8, true, _) => "S8R",
        (S::Shell { .. }, K::Quad8, false, _) => {
            return Err("full integration of quad8 shells (Abaqus's eight-node shell, S8R, uses reduced integration)".into());
        }
        (S::Shell { .. }, K::Tri6, ..) => return Err("tri6 shells (Abaqus has no general-purpose six-node triangular shell)".into()),
        (S::Beam { .. }, K::Beam2, ..) => "B31",
        (S::Rod { .. }, K::Rod2, ..) => "T3D2",
        (S::Rod { .. }, K::Rod3, ..) => "T3D3",
        (s, ..) => return Err(format!("{} elements in a {} section", kind.name(), s.name())),
    };
    Ok(if f.hybrid { format!("{plain}H") } else { plain.to_string() })
}

/// Whether the model has a nonlinear static step: hyperelastic materials are then written as
/// `*HYPERELASTIC` (in a deck of linear steps their elastic constants are the whole of it).
fn has_nonlinear_step(model: &Model) -> bool {
    model.steps().iter().any(|s| matches!(s.procedure, Procedure::NonlinearStatic(_)))
}

/// Every element, grouped by element type (from its section) and sorted by label.
fn elements_by_type<'m>(model: &'m Model, r: &Resolved) -> Result<BTreeMap<String, Vec<(u32, &'m Element)>>> {
    let mut by_type: BTreeMap<String, Vec<(u32, &Element)>> = BTreeMap::new();
    let mut bad = BTreeSet::new();
    let mut count = 0;
    for (name, section) in model.sections() {
        let mut by_kind: BTreeMap<ElementKind, Vec<(u32, &Element)>> = BTreeMap::new();
        for &e in r.sets.get(&section.set).map(|s| s.elements.as_slice()).unwrap_or_default() {
            if let Some(element) = model.element(e) {
                by_kind.entry(element.kind).or_default().push((e, element));
            }
        }
        for (kind, elements) in by_kind {
            // Trusses of a 2D model are the planar ones.
            let planar = |t: String| if r.is_planar() { t.replacen("T3D", "T2D", 1) } else { t };
            match element_type(kind, section).map(planar) {
                Ok(t) => {
                    count += elements.len();
                    by_type.entry(t).or_default().extend(elements);
                }
                Err(why) => {
                    bad.insert(format!("section `{name}`: {why}"));
                }
            }
        }
    }
    unsupported(FORMAT, bad)?;
    // A resolved model puts every element in exactly one section; this guards the writer.
    if count != model.elements().len() {
        return Err(Error::Unsupported { format: FORMAT, what: "elements without exactly one section".into() });
    }
    for elements in by_type.values_mut() {
        elements.sort_unstable_by_key(|(e, _)| *e);
    }
    Ok(by_type)
}

/// Everything else the format cannot express: case-only name clashes, orthotropic materials
/// outside solids, distributing couplings with a free translation.
fn check_expressible(model: &Model) -> Result<()> {
    let mut bad = BTreeSet::new();
    let namespaces: [(&str, Vec<&str>); 4] = [
        ("sets", model.sets().keys().map(String::as_str).collect()),
        ("materials", model.materials().keys().map(String::as_str).collect()),
        ("couplings", model.couplings().keys().map(String::as_str).collect()),
        ("steps", model.steps().iter().map(|s| s.name.as_str()).collect()),
    ];
    for (what, names) in namespaces {
        let mut seen: BTreeMap<String, &str> = BTreeMap::new();
        for n in names {
            if let Some(other) = seen.insert(n.to_ascii_uppercase(), n) {
                bad.insert(format!("{what} `{other}` and `{n}` differ only in case, and Abaqus ignores case in names"));
            }
        }
    }
    for (name, section) in model.sections() {
        let ortho = matches!(model.material(&section.material).and_then(|m| m.elastic), Some(Elastic::EngineeringConstants { .. }));
        if ortho && !section.kind.is_continuum() {
            bad.insert(format!(
                "section `{name}`: orthotropic material `{}` in a {} section (its material axes need an orientation, not supported yet)",
                section.material,
                section.kind.name()
            ));
        }
    }
    for step in model.steps().iter().filter(|s| matches!(s.procedure, Procedure::ExplicitDynamic(_))) {
        bad.insert(format!("step `{}`: explicit dynamic steps (not written to Abaqus input yet)", step.name));
    }
    if has_nonlinear_step(model)
        && let Some(material) = model.hyperelastic_in_use().into_iter().next()
    {
        for step in model.steps().iter().filter(|s| matches!(s.procedure, Procedure::Static {})) {
            bad.insert(format!(
                "step `{}`: a linear static step beside nonlinear steps of hyperelastic material `{material}` (Abaqus's general static step applies the hyperelastic law; write the linear step to its own deck)",
                step.name
            ));
        }
    }
    for step in model.steps() {
        if let Procedure::NonlinearStatic(n) = &step.procedure
            && n.load_path.is_some()
        {
            bad.insert(format!(
                "step `{}`: a load path (it would be an amplitude on every load and boundary condition of the step; not written yet)",
                step.name
            ));
        }
    }
    let general = model.steps().iter().filter(|s| matches!(s.procedure, Procedure::Static {} | Procedure::NonlinearStatic(_))).count();
    let plastic_in_use = model.sections().values().find(|s| model.material(&s.material).is_some_and(|m| m.plastic.is_some()));
    if has_nonlinear_step(model)
        && general > 1
        && let Some(section) = plastic_in_use
    {
        bad.insert(format!(
            "plastic material `{}` in a deck of {general} general steps (Abaqus carries the plastic state from one general step to the next, while every model step starts from the undeformed state; write each step to its own deck)",
            section.material
        ));
    }
    for (name, m) in model.materials() {
        if m.damage.is_some() || m.erosion.is_some() {
            bad.insert(format!(
                "material `{name}`: damage and erosion (Abaqus/Explicit failure models; this writer writes Abaqus/Standard procedures so far)"
            ));
        }
    }
    bad.extend(crate::thermal_items(model));
    for (name, c) in model.couplings() {
        if c.kind == CouplingKind::Distributing && !(1..=3).all(|d| c.dofs.contains(d)) {
            bad.insert(format!("coupling `{name}`: an Abaqus distributing coupling always couples translations 1-3"));
        }
    }
    for name in inclined_symmetries(model) {
        if !model.steps().iter().all(|s| s.boundary_conditions.contains(name)) {
            bad.insert(format!(
                "boundary condition `{name}`: a symmetry plane not normal to X, Y or Z active in some steps only (it is written as `*EQUATION`, which is model data and holds in every step)"
            ));
        }
    }
    unsupported(FORMAT, bad)
}

/// The symmetry conditions on inclined planes some step uses: linear equations rather than held
/// degrees of freedom.
fn inclined_symmetries(model: &Model) -> BTreeSet<&String> {
    model
        .boundary_conditions()
        .iter()
        .filter(|(n, bc)| {
            matches!(bc, BoundaryCondition::Symmetry { plane, .. } if plane.axis().is_none()) && !model.boundary_condition_users(n).is_empty()
        })
        .map(|(n, _)| n)
        .collect()
}

/// The boundary conditions of `names` that `*BOUNDARY` writes (inclined symmetry planes are
/// `*EQUATION`s).
fn held_conditions(model: &Model, names: &BTreeSet<String>) -> BTreeSet<String> {
    let inclined = inclined_symmetries(model);
    names.iter().filter(|n| !inclined.contains(n)).cloned().collect()
}

/// The property of a group of discrete elements.
enum DiscreteProperty {
    Mass(f64),
    Inertia([f64; 6]),
    /// `*SPRING` or `*DASHPOT`: the degrees of freedom line (empty for the axial form) and the
    /// stiffness or damping coefficient.
    Spring {
        dashpot: bool,
        dofs: String,
        value: f64,
    },
}

/// The discrete elements of one definition: element type, element set, members and property.
struct DiscreteGroup {
    element_type: &'static str,
    elset: String,
    members: Vec<(u32, Vec<u32>)>,
    property: DiscreteProperty,
}

/// Point masses, springs and dashpots as element groups (one per definition and element type,
/// labels counting up after the mesh's), or everything Abaqus cannot take.
fn discrete_groups<'r>(model: &Model, r: &'r Resolved) -> Result<Vec<DiscreteGroup>> {
    let mut next = model.elements().keys().next_back().map_or(1, |l| l.saturating_add(1));
    let mut label = || -> Result<u32> {
        let l = next;
        if l > fieldcraft_model::MAX_LABEL {
            return Err(Error::Limit(format!(
                "no element label left below {} for point masses, springs and dashpots",
                fieldcraft_model::MAX_LABEL as u64 + 1
            )));
        }
        next += 1;
        Ok(l)
    };
    let mut taken: BTreeSet<String> = model.sets().keys().map(|n| n.to_ascii_uppercase()).collect();
    let mut groups: Vec<DiscreteGroup> = Vec::new();
    // The (definition, element type) the last group holds.
    let mut current: Option<(&str, &'static str)> = None;
    let mut push = |name: &'r str, suffix: &str, element_type: &'static str, property: DiscreteProperty, member: (u32, Vec<u32>)| {
        let key = (name, element_type);
        if current.as_ref() == Some(&key)
            && let Some(g) = groups.last_mut()
        {
            g.members.push(member);
            return;
        }
        let elset = crate::unique_name(&format!("{name}-{suffix}"), &mut taken);
        groups.push(DiscreteGroup { element_type, elset, members: vec![member], property });
        current = Some(key);
    };
    let mut bad = BTreeSet::new();
    for m in &r.masses {
        if m.mass.iter().any(|c| *c != 0.0) {
            if m.mass[0] != m.mass[1] || m.mass[0] != m.mass[2] {
                bad.insert(format!("point mass `{}`: a mass that differs by direction (needs a local orientation, not written yet)", m.name));
            }
            push(&m.name, "mass", "MASS", DiscreteProperty::Mass(m.mass[0]), (label()?, vec![m.node]));
        }
    }
    for m in &r.masses {
        if m.inertia.iter().any(|c| *c != 0.0) {
            push(&m.name, "inertia", "ROTARYI", DiscreteProperty::Inertia(m.inertia), (label()?, vec![m.node]));
        }
    }
    for s in &r.springs {
        let what = s.kind();
        let line = matches!(s.direction, fieldcraft_model::ResolvedDirection::Line(_));
        let (element_type, dofs, nodes) = match (s.second, crate::axis_dof(s.direction), line) {
            (Some(b), _, true) => (if s.dashpot { "DASHPOTA" } else { "SPRINGA" }, String::new(), vec![s.first, b]),
            (Some(b), Some(d), false) => (if s.dashpot { "DASHPOT2" } else { "SPRING2" }, format!("{d}, {d}"), vec![s.first, b]),
            (None, Some(d), _) => (if s.dashpot { "DASHPOT1" } else { "SPRING1" }, d.to_string(), vec![s.first]),
            (_, None, _) => {
                bad.insert(format!("{what} `{}`: a direction along no global axis (needs a local orientation, not written yet)", s.name));
                continue;
            }
        };
        push(&s.name, what, element_type, DiscreteProperty::Spring { dashpot: s.dashpot, dofs, value: s.value }, (label()?, nodes));
    }
    unsupported(FORMAT, bad)?;
    Ok(groups)
}

/// Append a real number in the shortest form that reads back to the same value: `0.1`,
/// `210000.0`, `7.85E-9`, `1.0E300`.
fn push_number(out: &mut String, v: f64) -> Result<()> {
    if !v.is_finite() {
        return Err(Error::Unsupported { format: FORMAT, what: format!("the value {v} (not a finite number)") });
    }
    // Rust's `Debug` for f64 is the shortest round-trip form, switching to an exponent for very
    // large and small magnitudes; Abaqus wants the marker `E`, and a mantissa with a point reads
    // the same everywhere.
    let start = out.len();
    let _ = write!(out, "{v:?}");
    if let Some(at) = out.get(start..).and_then(|s| s.find('e')).map(|i| start + i) {
        let point = out.get(start..at).is_some_and(|m| m.contains('.'));
        out.replace_range(at..at + 1, if point { "E" } else { ".0E" });
    }
    Ok(())
}

/// [`push_number`] into a new string.
pub fn number(v: f64) -> Result<String> {
    let mut s = String::new();
    push_number(&mut s, v)?;
    Ok(s)
}

/// The surface face identifier of face `face` (1-based) of an element of `kind`.
pub fn face_id(kind: ElementKind, face: u8) -> Option<&'static str> {
    const SOLID: [&str; 6] = ["S1", "S2", "S3", "S4", "S5", "S6"];
    const SHELL: [&str; 2] = ["SPOS", "SNEG"];
    let i = usize::from(face).checked_sub(1)?;
    if i >= kind.faces().len() {
        return None;
    }
    match kind.dimension() {
        3 => SOLID.get(i).copied(),
        2 => SHELL.get(i).copied(),
        _ => None,
    }
}

/// The surface identifier of edge `edge` (1-based) of a 2D element: `S1`… in the model's edge
/// order (which is Abaqus's face order for plane and axisymmetric elements).
pub fn edge_id(kind: ElementKind, edge: u8) -> Option<&'static str> {
    const EDGES: [&str; 4] = ["S1", "S2", "S3", "S4"];
    let i = usize::from(edge).checked_sub(1)?;
    if kind.dimension() != 2 || i >= kind.edges().len() {
        return None;
    }
    EDGES.get(i).copied()
}

/// The node surface each coupling ties its reference point to: a fresh name per coupling that
/// clashes with no set (sets with faces own surfaces of their name) and no other coupling surface.
fn coupling_surfaces(model: &Model, r: &Resolved) -> BTreeMap<String, String> {
    let mut taken: BTreeSet<String> =
        r.sets.iter().filter(|(n, s)| !s.faces.is_empty() || !r.planar_edges(n).is_empty()).map(|(n, _)| n.to_ascii_uppercase()).collect();
    let mut out = BTreeMap::new();
    for name in model.couplings().keys() {
        out.insert(name.clone(), crate::unique_name(&format!("{name}-nodes"), &mut taken));
    }
    out
}

/// Consecutive runs of degrees of freedom: `[1, 2, 3, 5]` → `[(1, 3), (5, 5)]`.
fn dof_runs(dofs: &[u8]) -> Vec<(u8, u8)> {
    let mut runs: Vec<(u8, u8)> = Vec::new();
    for &d in dofs {
        match runs.last_mut() {
            Some((_, last)) if *last + 1 == d => *last = d,
            _ => runs.push((d, d)),
        }
    }
    runs
}

/// Whether a boundary condition holds every degree of freedom it constrains at zero.
fn holds_at_zero(bc: &BoundaryCondition) -> bool {
    (1..=6).all(|d| bc.value(d).is_none_or(|v| v == 0.0))
}

/// Where an element load acts: an element set, or one element.
#[derive(Clone, Copy)]
enum ElementTarget<'a> {
    Set(&'a str),
    Element(u32),
}

/// Gravity (an acceleration) or a body force (a force per volume).
#[derive(Clone, Copy, PartialEq)]
enum BodyKind {
    Gravity,
    Force,
}

/// A gravity or body force on some elements, before it is written.
struct ElementLoad<'a> {
    /// The element sets that hold exactly its elements.
    elsets: Vec<&'a str>,
    elements: &'a [u32],
    vector: [f64; 3],
}

impl<'a> ElementLoad<'a> {
    fn on_set(r: &'a Resolved, set: &'a str, vector: [f64; 3]) -> ElementLoad<'a> {
        ElementLoad { elsets: vec![set], elements: r.sets.get(set).map(|s| s.elements.as_slice()).unwrap_or_default(), vector }
    }
}

/// Whether no element is in two of the loads.
fn disjoint(loads: &[ElementLoad]) -> bool {
    if loads.len() < 2 {
        return true;
    }
    let mut all: Vec<u32> = loads.iter().flat_map(|l| l.elements.iter().copied()).collect();
    let n = all.len();
    all.sort_unstable();
    all.dedup();
    all.len() == n
}

/// Element loads of one kind: by element set when the loads share no element, otherwise summed
/// per element. Zero vectors are left out.
fn element_loads<'a>(loads: &[ElementLoad<'a>], kind: BodyKind, out: &mut Vec<(ElementTarget<'a>, BodyKind, [f64; 3])>) {
    if disjoint(loads) {
        for l in loads {
            out.extend(l.elsets.iter().map(|s| (ElementTarget::Set(s), kind, l.vector)));
        }
    } else {
        let mut summed: BTreeMap<u32, [f64; 3]> = BTreeMap::new();
        for l in loads {
            for &e in l.elements {
                for (c, add) in summed.entry(e).or_default().iter_mut().zip(l.vector) {
                    *c += add;
                }
            }
        }
        out.extend(summed.into_iter().map(|(e, v)| (ElementTarget::Element(e), kind, v)));
    }
    out.retain(|(_, _, v)| v.iter().any(|c| *c != 0.0));
}

/// A centrifugal load: the element set, ω², a point on the axis and its direction.
type Centrifugal<'a> = (&'a str, f64, [f64; 3], [f64; 3]);

/// The loads of one general step, as the deck writes them.
#[derive(Default)]
struct StepLoads<'a> {
    /// (node, degree of freedom) → summed value.
    concentrated: BTreeMap<(u32, u8), f64>,
    /// Surface (set name) → summed pressure.
    pressure: BTreeMap<&'a str, f64>,
    element: Vec<(ElementTarget<'a>, BodyKind, [f64; 3])>,
    centrifugal: Vec<Centrifugal<'a>>,
    /// Body forces of axisymmetric elements are radial and axial (`BR`, `BZ`).
    axisymmetric: bool,
}

/// Which load keywords a step wrote entries for, so the next step knows what to remove.
#[derive(Default, Clone, Copy)]
struct Written {
    concentrated: bool,
    pressure: bool,
    element: bool,
}

impl<'a> StepLoads<'a> {
    /// `everything` is every element and the section sets that hold them (gravity on the whole
    /// model).
    fn new(model: &'a Model, r: &'a Resolved, step: &'a Step, everything: &'a (Vec<&'a str>, Vec<u32>)) -> StepLoads<'a> {
        let axisymmetric = r.extent() == Some(fieldcraft_model::Extent::Revolved);
        let mut s = StepLoads { concentrated: crate::concentrated_loads(model, r, &step.loads), axisymmetric, ..StepLoads::default() };
        let (mut gravity, mut body) = (Vec::new(), Vec::new());
        for name in &step.loads {
            match model.load(name) {
                Some(Load::Pressure { set, magnitude }) => *s.pressure.entry(set.as_str()).or_default() += magnitude,
                Some(Load::Gravity { acceleration, set: Some(set) }) => gravity.push(ElementLoad::on_set(r, set, *acceleration)),
                Some(Load::Gravity { acceleration, set: None }) => {
                    gravity.push(ElementLoad { elsets: everything.0.clone(), elements: &everything.1, vector: *acceleration })
                }
                Some(Load::BodyForce { set, vector }) => body.push(ElementLoad::on_set(r, set, *vector)),
                Some(Load::Centrifugal { angular_velocity, point, axis, set }) => {
                    let w2 = angular_velocity * angular_velocity;
                    match set {
                        Some(set) => s.centrifugal.push((set, w2, *point, *axis)),
                        None => s.centrifugal.extend(everything.0.iter().map(|set| (*set, w2, *point, *axis))),
                    }
                }
                _ => {}
            }
        }
        s.pressure.retain(|_, p| *p != 0.0);
        element_loads(&gravity, BodyKind::Gravity, &mut s.element);
        element_loads(&body, BodyKind::Force, &mut s.element);
        s
    }

    fn written(&self) -> Written {
        Written {
            concentrated: !self.concentrated.is_empty(),
            pressure: !self.pressure.is_empty(),
            element: !self.element.is_empty() || !self.centrifugal.is_empty(),
        }
    }
}

#[derive(Default)]
struct Deck {
    out: String,
}

impl Deck {
    fn line(&mut self, s: &str) {
        self.out.push_str(s);
        self.out.push('\n');
    }

    fn comment(&mut self, s: &str) {
        self.out.push_str("** ");
        self.line(s);
    }

    /// `, value` for each value, then the end of the line.
    fn numbers_line(&mut self, head: &str, values: &[f64]) -> Result<()> {
        self.out.push_str(head);
        for (i, v) in values.iter().enumerate() {
            if i > 0 || !head.is_empty() {
                self.out.push_str(", ");
            }
            push_number(&mut self.out, *v)?;
        }
        self.out.push('\n');
        Ok(())
    }

    /// Labels, [`SET_LABELS_PER_LINE`] per line.
    fn labels(&mut self, labels: &[u32]) {
        for chunk in labels.chunks(SET_LABELS_PER_LINE) {
            for (i, l) in chunk.iter().enumerate() {
                let _ = write!(self.out, "{}{l}", if i > 0 { ", " } else { "" });
            }
            self.out.push('\n');
        }
    }

    fn heading(&mut self, units: Preset) {
        self.line("*HEADING");
        self.line(&format!("FieldCraft model, units {}", units.label()));
        self.comment(&format!("Written by FieldCraft in one consistent unit system: {}.", crate::units_in_words(units)));
    }

    fn nodes(&mut self, model: &Model) -> Result<()> {
        self.comment("Nodes");
        self.line("*NODE");
        for (label, node) in model.nodes() {
            self.numbers_line(&label.to_string(), &node.xyz)?;
        }
        if !model.reference_points().is_empty() {
            self.comment("Reference points");
            let mut points: Vec<_> = model.reference_points().iter().collect();
            points.sort_by_key(|(_, p)| p.label);
            self.line("*NODE");
            for (name, p) in points {
                self.comment(&format!("reference point {name}"));
                self.numbers_line(&p.label.to_string(), &p.xyz)?;
            }
        }
        Ok(())
    }

    fn elements(&mut self, by_type: &BTreeMap<String, Vec<(u32, &Element)>>) {
        self.comment("Elements");
        for (t, elements) in by_type {
            self.line(&format!("*ELEMENT, TYPE={t}"));
            for (e, element) in elements {
                let _ = write!(self.out, "{e}");
                for (i, n) in element.nodes.iter().enumerate() {
                    // A trailing comma continues the data line.
                    let continues = i >= FIRST_LINE_NODES && (i - FIRST_LINE_NODES).is_multiple_of(CONTINUATION_NODES);
                    let _ = write!(self.out, "{}{n}", if continues { ",\n" } else { ", " });
                }
                self.out.push('\n');
            }
        }
    }

    fn discrete_elements(&mut self, groups: &[DiscreteGroup]) {
        if groups.is_empty() {
            return;
        }
        self.comment("Point masses, springs and dashpots");
        for g in groups {
            self.line(&format!("*ELEMENT, TYPE={}, ELSET={}", g.element_type, g.elset));
            for (label, nodes) in &g.members {
                let nodes: Vec<String> = nodes.iter().map(u32::to_string).collect();
                self.line(&format!("{label}, {}", nodes.join(", ")));
            }
        }
    }

    fn discrete_properties(&mut self, groups: &[DiscreteGroup]) -> Result<()> {
        for g in groups {
            match &g.property {
                DiscreteProperty::Mass(m) => {
                    self.line(&format!("*MASS, ELSET={}", g.elset));
                    self.numbers_line("", &[*m])?;
                }
                DiscreteProperty::Inertia(i) => {
                    self.line(&format!("*ROTARY INERTIA, ELSET={}", g.elset));
                    self.numbers_line("", i)?;
                }
                DiscreteProperty::Spring { dashpot, dofs, value } => {
                    self.line(&format!("{}, ELSET={}", if *dashpot { "*DASHPOT" } else { "*SPRING" }, g.elset));
                    self.line(if dofs.is_empty() { "," } else { dofs });
                    self.numbers_line("", &[*value])?;
                }
            }
        }
        Ok(())
    }

    fn sets(&mut self, model: &Model, r: &Resolved) -> Result<()> {
        if r.sets.is_empty() {
            return Ok(());
        }
        self.comment("Sets");
        for (name, s) in &r.sets {
            if s.is_empty() {
                self.comment(&format!("set {name} resolves to nothing on this mesh"));
                continue;
            }
            if !s.nodes.is_empty() {
                self.line(&format!("*NSET, NSET={name}"));
                self.labels(&s.nodes);
            }
            if !s.elements.is_empty() {
                self.line(&format!("*ELSET, ELSET={name}"));
                self.labels(&s.elements);
            }
            let edges = r.planar_edges(name);
            if !s.faces.is_empty() || !edges.is_empty() {
                self.line(&format!("*SURFACE, TYPE=ELEMENT, NAME={name}"));
                for f in &s.faces {
                    let id = model.element(f.element).and_then(|e| face_id(e.kind, f.face)).ok_or_else(|| Error::Unsupported {
                        format: FORMAT,
                        what: format!("set `{name}`: element {} has no face {}", f.element, f.face),
                    })?;
                    let _ = writeln!(self.out, "{}, {id}", f.element);
                }
                for g in &edges {
                    let id = model.element(g.element).and_then(|e| edge_id(e.kind, g.edge)).ok_or_else(|| Error::Unsupported {
                        format: FORMAT,
                        what: format!("set `{name}`: element {} has no edge {}", g.element, g.edge),
                    })?;
                    let _ = writeln!(self.out, "{}, {id}", g.element);
                }
            }
        }
        Ok(())
    }

    fn materials(&mut self, model: &Model) -> Result<()> {
        if model.materials().is_empty() {
            return Ok(());
        }
        self.comment("Materials");
        for (name, m) in model.materials() {
            self.line(&format!("*MATERIAL, NAME={name}"));
            match m.elastic {
                Some(Elastic::Isotropic { young, poisson })
                    if has_nonlinear_step(model)
                        && let Some(h) = m.hyperelastic =>
                {
                    // The polynomial constants of the option and U = (J − 1)²/D1 with D1 = 2/κ.
                    let (shear, bulk) = Hyperelastic::moduli(young, poisson);
                    let [c10, c01, c20, c30] = h.constants(shear);
                    let (keyword, values) = match h {
                        Hyperelastic::NeoHookean {} => ("NEO HOOKE", vec![c10, 2.0 / bulk]),
                        Hyperelastic::MooneyRivlin { .. } => ("MOONEY-RIVLIN", vec![c10, c01, 2.0 / bulk]),
                        Hyperelastic::Yeoh { .. } => ("YEOH", vec![c10, c20, c30, 2.0 / bulk]),
                    };
                    self.line(&format!("*HYPERELASTIC, {keyword}"));
                    self.numbers_line("", &values)?;
                }
                Some(Elastic::Isotropic { young, poisson }) => {
                    self.line("*ELASTIC, TYPE=ISOTROPIC");
                    self.numbers_line("", &[young, poisson])?;
                }
                Some(Elastic::EngineeringConstants { young, poisson, shear }) => {
                    // E1, E2, E3, ν12, ν13, ν23, G12, G13 on the first line, G23 on the second.
                    self.line("*ELASTIC, TYPE=ENGINEERING CONSTANTS");
                    let [g12, g13, g23] = shear;
                    self.numbers_line("", &[young[0], young[1], young[2], poisson[0], poisson[1], poisson[2], g12, g13])?;
                    self.numbers_line("", &[g23])?;
                }
                None => {}
            }
            if let Some(rho) = m.density {
                self.line("*DENSITY");
                self.numbers_line("", &[rho])?;
            }
            if let Some(p) = &m.plastic {
                self.plastic(p)?;
            }
        }
        Ok(())
    }

    /// Von Mises plasticity (see the module documentation for the mapping).
    fn plastic(&mut self, p: &Plastic) -> Result<()> {
        match p {
            Plastic::Bilinear { yield_stress, hardening_modulus, hardening } => {
                if *hardening_modulus == 0.0 {
                    self.line("*PLASTIC");
                    self.numbers_line("", &[*yield_stress, 0.0])?;
                } else if *hardening == Hardening::Kinematic {
                    self.line("*PLASTIC, HARDENING=KINEMATIC");
                    self.numbers_line("", &[*yield_stress, 0.0])?;
                    self.numbers_line("", &[yield_stress + hardening_modulus, 1.0])?;
                } else {
                    self.line("*PLASTIC");
                    self.numbers_line("", &[*yield_stress, 0.0])?;
                    self.numbers_line("", &[yield_stress + hardening_modulus * BILINEAR_REACH, BILINEAR_REACH])?;
                }
            }
            Plastic::Table { points } => {
                self.line("*PLASTIC");
                for point in points {
                    self.numbers_line("", &[point.stress, point.strain])?;
                }
            }
            Plastic::JohnsonCook { a, b, n, c, m, reference_strain_rate, room_temperature, melt_temperature, specific_heat, taylor_quinney } => {
                self.line("*PLASTIC, HARDENING=JOHNSON COOK");
                self.numbers_line("", &[*a, *b, *n, *m, *melt_temperature, *room_temperature])?;
                if *c != 0.0 {
                    self.line("*RATE DEPENDENT, TYPE=JOHNSON COOK");
                    self.numbers_line("", &[*c, *reference_strain_rate])?;
                }
                if let Some(cp) = specific_heat {
                    self.line("*SPECIFIC HEAT");
                    self.numbers_line("", &[*cp])?;
                    self.line("*INELASTIC HEAT FRACTION");
                    self.numbers_line("", &[*taylor_quinney])?;
                }
            }
        }
        Ok(())
    }

    fn sections(&mut self, model: &Model) -> Result<()> {
        if model.sections().is_empty() {
            return Ok(());
        }
        self.comment("Sections");
        for (name, s) in model.sections() {
            self.comment(&format!("section {name}"));
            let (set, material) = (&s.set, &s.material);
            match s.kind {
                SectionKind::Solid {} | SectionKind::Axisymmetric {} => {
                    self.line(&format!("*SOLID SECTION, ELSET={set}, MATERIAL={material}"));
                    self.line(",");
                }
                SectionKind::PlaneStress { thickness } | SectionKind::PlaneStrain { thickness } => {
                    self.line(&format!("*SOLID SECTION, ELSET={set}, MATERIAL={material}"));
                    self.numbers_line("", &[thickness])?;
                }
                SectionKind::Rod { area } => {
                    self.line(&format!("*SOLID SECTION, ELSET={set}, MATERIAL={material}"));
                    self.numbers_line("", &[area])?;
                }
                SectionKind::Shell { thickness, offset, points } => {
                    let offset = if offset == 0.0 { String::new() } else { format!(", OFFSET={}", number(offset)?) };
                    self.line(&format!("*SHELL SECTION, ELSET={set}, MATERIAL={material}{offset}"));
                    push_number(&mut self.out, thickness)?;
                    let _ = writeln!(self.out, ", {points}");
                }
                SectionKind::Beam { profile, orientation } => {
                    // Abaqus profiles: RECT a, b; CIRC r; PIPE r, t; I l, h, b1, b2, t1, t2, t3
                    // (l: from the section origin to the bottom; centred here); BOX a, b, t1–t4.
                    let (shape, dims) = match profile {
                        BeamProfile::Rectangle { width, height } => ("RECT", vec![width, height]),
                        BeamProfile::Circle { radius } => ("CIRC", vec![radius]),
                        BeamProfile::Pipe { radius, thickness } => ("PIPE", vec![radius, thickness]),
                        BeamProfile::I { height, width, web, flange } => ("I", vec![height / 2.0, height, width, width, flange, flange, web]),
                        BeamProfile::Box { width, height, thickness } => ("BOX", vec![width, height, thickness, thickness, thickness, thickness]),
                    };
                    self.line(&format!("*BEAM SECTION, SECTION={shape}, ELSET={set}, MATERIAL={material}"));
                    self.numbers_line("", &dims)?;
                    self.numbers_line("", &orientation)?;
                }
            }
        }
        Ok(())
    }

    fn couplings(&mut self, model: &Model, r: &Resolved) {
        if model.couplings().is_empty() {
            return;
        }
        let surfaces = coupling_surfaces(model, r);
        self.comment("Couplings: every node of the set weighted equally");
        for (name, c) in model.couplings() {
            let (Some(surface), Some(node)) = (surfaces.get(name), r.reference_point_nodes.get(&c.reference_point)) else { continue };
            self.line(&format!("*SURFACE, TYPE=NODE, NAME={surface}"));
            self.line(&format!("{}, 1.0", c.set));
            self.line(&format!("*COUPLING, CONSTRAINT NAME={name}, REF NODE={node}, SURFACE={surface}"));
            self.line(match c.kind {
                CouplingKind::Kinematic => "*KINEMATIC",
                CouplingKind::Distributing => "*DISTRIBUTING, WEIGHTING METHOD=UNIFORM",
            });
            for (a, b) in dof_runs(&c.dofs.list()) {
                self.line(&format!("{a}, {b}"));
            }
        }
    }

    /// One surface interaction (the friction) per contact, then its pair, and every general
    /// contact in the model's one `*CONTACT` definition.
    fn contacts(&mut self, model: &Model, r: &Resolved) -> Result<()> {
        if model.contacts().is_empty() {
            return Ok(());
        }
        let mut taken: BTreeSet<String> = model.sets().keys().map(|n| n.to_ascii_uppercase()).collect();
        taken.extend(coupling_surfaces(model, r).into_values().map(|n| n.to_ascii_uppercase()));
        self.comment("Contact: a surface interaction (Coulomb friction) per contact, contact pairs, general contact");
        let mut general = Vec::new();
        for (name, c) in model.contacts() {
            let interaction = crate::unique_name(&format!("{name}-interaction"), &mut taken);
            self.line(&format!("*SURFACE INTERACTION, NAME={interaction}"));
            let f = c.friction;
            if f.decay > 0.0 {
                self.line("*FRICTION, EXPONENTIAL DECAY");
                self.numbers_line("", &[f.static_coefficient, f.kinetic_coefficient(), f.decay])?;
            } else {
                self.line("*FRICTION");
                self.numbers_line("", &[f.static_coefficient])?;
            }
            if c.penalty_scale != 1.0 {
                self.comment(&format!(
                    "contact `{name}`: penalty scale {} of FieldCraft's explicit penalty stiffness; Abaqus chooses its own",
                    c.penalty_scale
                ));
            }
            match &c.surfaces {
                ContactSurfaces::Pair { main, secondary } => {
                    // ADJUST=0.0 moves the secondary nodes that start inside the main surface
                    // onto it, strain-free; Abaqus has no refusal of initial overclosure.
                    let adjust = match c.initial_penetration {
                        InitialPenetration::Adjust => ", ADJUST=0.0",
                        InitialPenetration::Refuse => {
                            self.comment(&format!(
                                "contact `{name}`: FieldCraft refuses initially overclosed nodes; Abaqus resolves them in the first increment"
                            ));
                            ""
                        }
                    };
                    self.line(&format!("*CONTACT PAIR, INTERACTION={interaction}, TYPE=SURFACE TO SURFACE{adjust}"));
                    self.line(&format!("{secondary}, {main}"));
                }
                ContactSurfaces::General { surface } => general.push((surface.clone(), interaction)),
            }
        }
        if general.is_empty() {
            return Ok(());
        }
        // General contact: every surface against itself (a blank second surface), or every
        // exterior face; each with its interaction.
        self.line("*CONTACT");
        if general.iter().any(|(s, _)| s.is_none()) {
            self.line("*CONTACT INCLUSIONS, ALL EXTERIOR");
        }
        let named: Vec<&str> = general.iter().filter_map(|(s, _)| s.as_deref()).collect();
        if !named.is_empty() {
            self.line("*CONTACT INCLUSIONS");
            for s in &named {
                self.line(&format!("{s},"));
            }
        }
        self.line("*CONTACT PROPERTY ASSIGNMENT");
        for (surface, interaction) in &general {
            match surface {
                Some(s) => self.line(&format!("{s}, {s}, {interaction}")),
                None => self.line(&format!(", , {interaction}")),
            }
        }
        Ok(())
    }

    /// Ties: `*TIE` per tie, secondary surface first.
    fn ties(&mut self, model: &Model) -> Result<()> {
        if model.ties().is_empty() {
            return Ok(());
        }
        self.comment("Ties: each secondary surface bonded to its main surface");
        for (name, t) in model.ties() {
            let tolerance = t.position_tolerance.map_or(String::new(), |v| format!(", POSITION TOLERANCE={v:?}"));
            self.line(&format!("*TIE, NAME={name}, ADJUST=NO{tolerance}"));
            self.line(&format!("{}, {}", t.secondary, t.main));
        }
        Ok(())
    }

    fn boundary_lines<'n>(&mut self, model: &Model, r: &Resolved, names: impl IntoIterator<Item = &'n String>) -> Result<()> {
        for name in names {
            let Some(bc) = model.boundary_condition(name) else { continue };
            let target = match bc.target() {
                Target::Set(s) => s.clone(),
                Target::ReferencePoint(p) => match r.reference_point_nodes.get(p) {
                    Some(n) => n.to_string(),
                    None => continue,
                },
            };
            self.comment(&format!("boundary condition {name}"));
            match bc {
                BoundaryCondition::Fixed { .. } => self.line(&format!("{target}, ENCASTRE")),
                BoundaryCondition::Symmetry { plane, .. } => {
                    let s = match plane.axis() {
                        Some(0) => "XSYMM",
                        Some(1) => "YSYMM",
                        Some(2) => "ZSYMM",
                        // An inclined plane is an `*EQUATION` (model data), never here.
                        _ => continue,
                    };
                    self.line(&format!("{target}, {s}"));
                }
                // Refused up front (check_expressible).
                BoundaryCondition::Temperature { .. } => {}
                BoundaryCondition::Displacement { values, .. } => {
                    for (dof, v) in (1u8..).zip(values) {
                        if let Some(v) = v {
                            self.numbers_line(&format!("{target}, {dof}, {dof}"), &[*v])?;
                        }
                    }
                }
            }
        }
        Ok(())
    }

    /// The inclined symmetry planes: one `*EQUATION` per node and held direction (n · u = 0, and
    /// the rotations about the plane's axes at shell and beam nodes), whose first term Abaqus
    /// eliminates; it is chosen among the degrees of freedom no boundary condition, coupling
    /// or tie constrains ([`crate::symmetry_equations`]).
    fn equations(&mut self, model: &Model, r: &Resolved) -> Result<()> {
        let inclined = inclined_symmetries(model);
        if inclined.is_empty() {
            return Ok(());
        }
        let steps = model.steps();
        // Per node and degree of freedom: the value it is held at in each step that holds it.
        let mut held: BTreeMap<(u32, u8), Vec<f64>> = BTreeMap::new();
        for step in steps {
            for name in held_conditions(model, &step.boundary_conditions) {
                let Some(bc) = model.boundary_condition(&name) else { continue };
                for &node in r.target_nodes(bc.target()) {
                    for dof in bc.dofs().list() {
                        let v = if matches!(step.procedure, Procedure::Frequency { .. } | Procedure::Buckle { .. }) {
                            0.0
                        } else {
                            bc.value(dof).unwrap_or(0.0)
                        };
                        held.entry((node, dof)).or_default().push(v);
                    }
                }
            }
        }
        let mut dependent: BTreeSet<(u32, u8)> = BTreeSet::new();
        for c in model.couplings().values().filter(|c| c.kind == CouplingKind::Kinematic) {
            for &node in r.sets.get(&c.set).map(|s| s.nodes.as_slice()).unwrap_or_default() {
                dependent.extend(c.dofs.list().into_iter().map(|d| (node, d)));
            }
        }
        for t in model.ties().values() {
            for &node in r.sets.get(&t.secondary).map(|s| s.nodes.as_slice()).unwrap_or_default() {
                dependent.extend((1..=6).map(|d| (node, d)));
            }
        }
        let state = |node: u32, dof: u8| {
            if dependent.contains(&(node, dof)) {
                return crate::Held::Sometimes;
            }
            match held.get(&(node, dof)) {
                None => crate::Held::Free,
                Some(v) if v.len() == steps.len() && v.iter().all(|x| *x == v[0]) => crate::Held::Always(v[0]),
                Some(_) => crate::Held::Sometimes,
            }
        };
        let mut problems = BTreeSet::new();
        let equations = crate::symmetry_equations(model, r, inclined, &model.rotating_nodes(), state, &mut problems);
        unsupported(FORMAT, problems)?;
        self.comment("Symmetry planes not normal to X, Y or Z: per node, no motion along the normal (and no rotation about the plane's axes)");
        for e in &equations {
            self.line("*EQUATION");
            self.line(&e.terms.len().to_string());
            let mut fields = Vec::new();
            for (dof, c) in &e.terms {
                fields.push(format!("{}, {dof}", e.node));
                fields.push(number(*c)?);
            }
            self.line(&fields.join(", "));
        }
        Ok(())
    }

    /// Boundary conditions of a step against the state it starts from.
    fn step_boundaries(&mut self, model: &Model, r: &Resolved, base: &BTreeSet<String>, active: &BTreeSet<String>) -> Result<()> {
        let (base, active) = (&held_conditions(model, base), &held_conditions(model, active));
        if active == base {
            self.comment("boundary conditions: as in the base state");
            Ok(())
        } else if active.is_superset(base) {
            self.line("*BOUNDARY");
            self.boundary_lines(model, r, active.difference(base))
        } else {
            self.line("*BOUNDARY, OP=NEW");
            self.boundary_lines(model, r, active)
        }
    }

    /// The keyword of one load family with `OP=NEW`, when the step has entries for it or the
    /// previous one had some to remove.
    fn op_new(&mut self, keyword: &str, now: bool, before: bool) -> bool {
        if now || before {
            self.line(&format!("{keyword}, OP=NEW"));
        }
        now
    }

    fn step_loads(&mut self, loads: &StepLoads, before: Written) -> Result<()> {
        if self.op_new("*CLOAD", !loads.concentrated.is_empty(), before.concentrated) {
            for ((node, dof), v) in &loads.concentrated {
                self.numbers_line(&format!("{node}, {dof}"), &[*v])?;
            }
        }
        if self.op_new("*DSLOAD", !loads.pressure.is_empty(), before.pressure) {
            for (surface, p) in &loads.pressure {
                self.numbers_line(&format!("{surface}, P"), &[*p])?;
            }
        }
        if self.op_new("*DLOAD", !loads.element.is_empty() || !loads.centrifugal.is_empty(), before.element) {
            for (target, kind, v) in &loads.element {
                let target = match target {
                    ElementTarget::Set(s) => s.to_string(),
                    ElementTarget::Element(e) => e.to_string(),
                };
                match kind {
                    BodyKind::Gravity => {
                        let g = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
                        self.numbers_line(&format!("{target}, GRAV"), &[g, v[0] / g, v[1] / g, v[2] / g])?;
                    }
                    BodyKind::Force => {
                        let labels: &[&str] = if loads.axisymmetric { &["BR", "BZ"] } else { &["BX", "BY", "BZ"] };
                        for (label, c) in labels.iter().zip(v) {
                            if *c != 0.0 {
                                self.numbers_line(&format!("{target}, {label}"), &[*c])?;
                            }
                        }
                    }
                }
            }
            for (set, w2, point, axis) in &loads.centrifugal {
                // CENTRIF: ω², a point on the axis and the axis direction.
                self.numbers_line(&format!("{set}, CENTRIF"), &[*w2, point[0], point[1], point[2], axis[0], axis[1], axis[2]])?;
            }
        }
        Ok(())
    }

    fn output(&mut self, step: &Step) {
        let pick = |names: [(OutputVariable, &'static str); 2]| -> Vec<&'static str> {
            names.into_iter().filter(|(v, _)| step.output.contains(v)).map(|(_, n)| n).collect()
        };
        let node = pick([(OutputVariable::Displacement, "U"), (OutputVariable::ReactionForce, "RF")]);
        let element = pick([(OutputVariable::Stress, "S"), (OutputVariable::Strain, "E")]);
        if node.is_empty() && element.is_empty() {
            self.line("*OUTPUT, FIELD, OP=NEW, FREQUENCY=0");
            return;
        }
        self.line("*OUTPUT, FIELD, OP=NEW");
        if !node.is_empty() {
            self.line("*NODE OUTPUT");
            self.line(&node.join(", "));
        }
        if !element.is_empty() {
            self.line("*ELEMENT OUTPUT");
            self.line(&element.join(", "));
        }
    }

    fn steps(&mut self, model: &Model, r: &Resolved) -> Result<()> {
        let steps = model.steps();
        // The initial state: zero-valued boundary conditions active in every step.
        let initial: BTreeSet<String> = if steps.is_empty() {
            BTreeSet::new()
        } else {
            let in_every_step = |n: &String| steps.iter().all(|s| s.boundary_conditions.contains(n));
            let all = model.boundary_conditions().iter().filter(|(n, bc)| holds_at_zero(bc) && in_every_step(n)).map(|(n, _)| n.clone()).collect();
            held_conditions(model, &all)
        };
        if !initial.is_empty() {
            self.comment("Initial state: boundary conditions active in every step");
            self.line("*BOUNDARY");
            self.boundary_lines(model, r, &initial)?;
        }
        // Gravity on the whole model acts on every element, through the section sets.
        let whole_model_gravity = model.loads().values().any(|l| matches!(l, Load::Gravity { set: None, .. } | Load::Centrifugal { set: None, .. }));
        let everything: (Vec<&str>, Vec<u32>) = if whole_model_gravity {
            (model.sections().values().map(|s| s.set.as_str()).collect(), model.elements().keys().copied().collect())
        } else {
            Default::default()
        };
        let mut base_bcs = initial;
        let mut base_loads: BTreeSet<String> = BTreeSet::new();
        let mut written = Written::default();
        for step in steps {
            self.comment(&format!("Step {}", step.name));
            match step.procedure {
                Procedure::Fluid { .. } => return Err(crate::fluid_step(FORMAT, &step.name)),
                Procedure::Static {} | Procedure::NonlinearStatic(_) => {
                    if let Procedure::NonlinearStatic(n) = &step.procedure {
                        if n.imperfection.is_some() {
                            return Err(Error::Unsupported {
                                format: FORMAT,
                                what: format!(
                                    "step `{}`: an imperfection from a buckling mode (FieldCraft computes the mode when it solves)",
                                    step.name
                                ),
                            });
                        }
                        let nlgeom = if n.nonlinear_geometry { "YES" } else { "NO" };
                        self.line(&format!("*STEP, NAME={}, NLGEOM={nlgeom}, INC={}", step.name, n.max_increments));
                        if let Some(arc) = &n.arc_length {
                            // Riks: initial arc, period, minimum and maximum arcs, the load
                            // proportionality factor and the displacement that end the step.
                            if arc.constraint != fieldcraft_model::ArcConstraint::Cylindrical {
                                return Err(Error::Unsupported {
                                    format: FORMAT,
                                    what: format!("step `{}`: a spherical arc-length constraint", step.name),
                                });
                            }
                            self.line("*STATIC, RIKS");
                            self.numbers_line("", &[n.initial_increment, 1.0, n.min_increment, n.max_increment, arc.max_load_factor])?;
                            if let Some(l) = &arc.displacement_limit {
                                self.out.pop();
                                let _ = write!(self.out, ", {}, {}, ", l.node, l.direction);
                                push_number(&mut self.out, l.displacement)?;
                                self.out.push('\n');
                            }
                        } else {
                            self.line("*STATIC");
                            self.numbers_line("", &[n.initial_increment, 1.0, n.min_increment, n.max_increment])?;
                        }
                        if n.line_search {
                            self.line("*CONTROLS, PARAMETERS=LINE SEARCH");
                            self.line("5,");
                        }
                        self.comment("iterations and convergence: Abaqus's own controls (its residual norms differ from FieldCraft's)");
                    } else {
                        self.line(&format!("*STEP, NAME={}, NLGEOM=NO", step.name));
                        self.line("*STATIC");
                        self.line("1.0, 1.0");
                    }
                    self.step_boundaries(model, r, &base_bcs, &step.boundary_conditions)?;
                    if step.loads == base_loads {
                        self.comment("loads: as in the previous general step");
                    } else {
                        let loads = StepLoads::new(model, r, step, &everything);
                        self.step_loads(&loads, written)?;
                        written = loads.written();
                    }
                    base_bcs.clone_from(&step.boundary_conditions);
                    base_loads.clone_from(&step.loads);
                }
                Procedure::Buckle { modes, shift } => {
                    // Abaqus's buckling base state is the end of the previous general step;
                    // FieldCraft's steps all start from the undeformed, unloaded model.
                    if !base_loads.is_empty()
                        || steps
                            .iter()
                            .take_while(|s| !std::ptr::eq(*s, step))
                            .any(|s| !matches!(s.procedure, Procedure::Frequency { .. } | Procedure::Buckle { .. }))
                    {
                        return Err(Error::Unsupported {
                            format: FORMAT,
                            what: format!(
                                "step `{}`: a buckling step after a general step (Abaqus would buckle from that step's loaded state)",
                                step.name
                            ),
                        });
                    }
                    if shift.is_some() {
                        return Err(Error::Unsupported { format: FORMAT, what: format!("step `{}`: a buckling shift", step.name) });
                    }
                    self.line(&format!("*STEP, NAME={}", step.name));
                    self.line("*BUCKLE, EIGENSOLVER=LANCZOS");
                    self.line(&format!("{modes}"));
                    self.step_boundaries(model, r, &base_bcs, &step.boundary_conditions)?;
                    let loads = StepLoads::new(model, r, step, &everything);
                    self.step_loads(&loads, written)?;
                    written = loads.written();
                }
                Procedure::Frequency { modes, min_hz, max_hz } => {
                    self.line(&format!("*STEP, NAME={}", step.name));
                    self.line("*FREQUENCY, EIGENSOLVER=LANCZOS");
                    // Modes, then the frequency range; empty trailing fields are left off.
                    let _ = write!(self.out, "{modes}");
                    let range: Vec<Option<f64>> = [min_hz, max_hz].into_iter().rev().skip_while(Option::is_none).collect();
                    for f in range.into_iter().rev() {
                        self.out.push_str(", ");
                        if let Some(f) = f {
                            push_number(&mut self.out, f)?;
                        }
                    }
                    self.out.push('\n');
                    self.step_boundaries(model, r, &base_bcs, &step.boundary_conditions)?;
                    if !step.loads.is_empty() {
                        self.comment("loads have no effect on natural frequencies and are not written");
                    }
                }
                Procedure::ExplicitDynamic(_) | Procedure::HeatTransfer { .. } | Procedure::CoupledTemperatureDisplacement {} => {
                    return Err(Error::Unsupported { format: FORMAT, what: format!("step `{}`: {} steps", step.name, step.procedure.name()) });
                }
            }
            self.output(step);
            self.line("*END STEP");
        }
        Ok(())
    }
}
