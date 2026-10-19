//! Nastran bulk data: the model as an input file for a linear static (SOL 101) or normal modes
//! (SOL 103) run, written from the public Quick Reference Guide in our own words.
//!
//! A deck has three parts: executive control (the solution sequence), case control (one subcase
//! per step, selecting its constraint set, load set or eigenvalue method and its output) and bulk
//! data (the entries). What the model holds becomes:
//!
//! - nodes and reference points: `GRID`;
//! - solids: `CTETRA` (tet4, tet10), `CPENTA` (wedge6), `CPYRAM` (pyramid5), `CHEXA` (hex8, hex20) with `PSOLID`;
//!   shells: `CTRIA3`, `CQUAD4`, `CTRIA6`, `CQUAD8` with `PSHELL`; beams: `CBEAM` with `PBEAML`;
//!   rods: `CROD` with `PROD`. Nastran numbers a 20-node hexahedron's mid-side nodes bottom
//!   edges, vertical edges, top edges, where the model has bottom, top, vertical
//!   ([`grids`] reorders them); every other kind keeps the model's order;
//! - the element formulation: `PSOLID`'s integration network and scheme (see [`solid_rule`]);
//!   formulations Nastran's linear elements have no counterpart for are reported, not
//!   approximated;
//! - isotropic materials: `MAT1`, with `MATS1` for von Mises plasticity (bilinear by its
//!   hardening slope, a hardening table as a `TABLES1` stress–strain curve; Johnson–Cook, damage
//!   and erosion are reported, not written); a linear solution reads the elastic entry alone;
//! - orthotropic engineering constants: `MAT9` (the stiffness matrix in the basic axes) for
//!   solids and `MAT8` for shells lying in planes of constant Z, with the material X axis tied to
//!   the basic X axis through a `CORD2R`;
//! - boundary conditions: per step, one constraint set of `SPC1` (held degrees of freedom) and
//!   `SPC` (prescribed values); in a modes step every constrained degree of freedom is held;
//!   symmetry about a plane not normal to X, Y or Z: one `MPC` set per step, an equation per node
//!   and held direction (n · u = 0, and t · θ = 0 about the plane's two axes at shell and beam
//!   grids), its dependent component one no `SPC` or rigid element constrains;
//! - forces, moments and tractions: `FORCE` and `MOMENT` at the nodal values the model resolves,
//!   so every writer and the solver share one distribution; pressure: `PLOAD4` per element face;
//!   gravity: `GRAV`; a step with several loads combines them with `LOAD`;
//! - kinematic coupling: `RBE2` (the reference point independent); distributing coupling: `RBE3`
//!   (equal weights on the set's translations; the reference point dependent);
//! - frequency steps: `EIGRL` (Lanczos) with the step's mode count and range;
//! - 2D models: plane stress as membranes, `CTRIA3`, `CQUAD4`, `CTRIA6`, `CQUAD8` with a
//!   `PSHELL` that has only a membrane material, every grid holding its out-of-plane degrees of
//!   freedom 3–6 (`GRID` field PS); edge pressures as the consistent `FORCE`s the model's edge
//!   integration gives. Plane strain and axisymmetric elements are not written: their cards
//!   (`PPLANE`/`PLPLANE` with `CPLSTN*`/`CQUADX`, `CTRIAX6`, `PAXI`) differ between Nastran
//!   dialects, so a deck could not be read the same everywhere;
//! - centrifugal loads: `RFORCE` about an axis through the origin of the basic system (speed in
//!   revolutions per time unit);
//! - point masses: `CONM1`, the node's 6 × 6 mass matrix (the masses per direction on the
//!   translations, the inertia tensor's components on the rotations, so no product-of-inertia sign
//!   convention is involved); springs and dashpots along a global degree of freedom: `CELAS2` and
//!   `CDAMP2` (the second grid left blank for one to the ground). Their ids follow the rigid
//!   elements'. A spring or dashpot along a skew direction and the three-node truss (rod3) have
//!   no counterpart here and are reported.
//!
//! Field widths, continuations and number forms are in [`card`]. Model content Nastran cannot
//! express is reported as [`Error::Unsupported`], every item at once, and nothing is written.
//! The same model always gives the same bytes.

mod card;

use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;

use fieldcraft_model::{
    BeamProfile, Coupling, CouplingKind, Dofs, Elastic, ElementKind, Formulation, Hardening, Integration, Load, MAX_LABEL, Model, OutputVariable,
    Plastic, Procedure, Resolved, SectionKind, Step,
};

use fieldcraft_model::units::Preset;

use crate::{Error, Result};
use card::Card;

/// The format's name in errors.
pub const FORMAT: &str = "Nastran";

/// Most `GRAV` entries one run may select.
const MAX_GRAV: usize = 9;

/// Grids per `SPC1` entry; longer lists are split over entries of the same set.
const SPC1_GRIDS: usize = 1000;

/// The model as a Nastran deck, with every step (a deck of some steps:
/// [`crate::Format::write_steps`]).
pub fn write(model: &Model, units: Preset) -> Result<String> {
    let expressed = crate::express(model, units)?;
    write_resolved(&expressed, &expressed.resolve()?, units)
}

/// [`write`] for a model already expressed in `units` ([`crate::express`]) and resolved.
pub fn write_resolved(model: &Model, r: &Resolved, units: Preset) -> Result<String> {
    crate::refuse_fluid(FORMAT, model)?;
    let mut deck = Deck::new(model, r, model.steps().iter().collect());
    deck.plan()?;
    deck.emit(units)
}

/// The Nastran grid order of an element's nodes (model order in, Nastran order out).
pub fn grids(kind: ElementKind, nodes: &[u32]) -> impl Iterator<Item = u32> + '_ {
    order(kind).iter().filter_map(|&i| nodes.get(i).copied())
}

/// For each Nastran grid position, the position in the model's node list.
pub fn order(kind: ElementKind) -> &'static [usize] {
    const IDENTITY: [usize; 20] = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19];
    // Model: 9–12 bottom edges, 13–16 top edges, 17–20 vertical edges. Nastran: bottom, vertical, top.
    const HEX20: [usize; 20] = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 16, 17, 18, 19, 12, 13, 14, 15];
    match kind {
        ElementKind::Hex20 => &HEX20,
        k => IDENTITY.get(..k.node_count()).unwrap_or(&IDENTITY),
    }
}

/// The connection entry's name for an element kind.
pub fn element_card(kind: ElementKind) -> &'static str {
    match kind {
        ElementKind::Rod2 | ElementKind::Rod3 => "CROD",
        ElementKind::Beam2 => "CBEAM",
        ElementKind::Tri3 => "CTRIA3",
        ElementKind::Tri6 => "CTRIA6",
        ElementKind::Quad4 => "CQUAD4",
        ElementKind::Quad8 => "CQUAD8",
        ElementKind::Tet4 | ElementKind::Tet10 => "CTETRA",
        ElementKind::Wedge6 => "CPENTA",
        ElementKind::Pyramid5 => "CPYRAM",
        ElementKind::Hex8 | ElementKind::Hex20 => "CHEXA",
    }
}

/// `PSOLID`'s integration network (`IN`) and scheme (`ISOP`) for an element kind and
/// formulation, or why Nastran has none.
///
/// - Full integration is the standard isoparametric rule: two points per direction for linear
///   kinds, three for quadratic ones (`TWO`/`THREE`, `FULL`).
/// - Reduced integration of a 20-node hexahedron is the 2 × 2 × 2 standard rule.
/// - Incompatible modes of an 8-node hexahedron are Nastran's bubble functions (`BUBBLE`, with
///   its default reduced-shear scheme).
/// - One-point reduced integration of linear solids (with hourglass control) and the hybrid
///   (displacement–pressure) formulation have no linear Nastran counterpart.
pub fn solid_rule(kind: ElementKind, f: Formulation) -> std::result::Result<(Option<&'static str>, Option<&'static str>), String> {
    if f.hybrid {
        return Err("the hybrid (displacement–pressure) formulation".into());
    }
    let quadratic = kind.is_quadratic();
    match (kind, f.incompatible_modes, f.integration) {
        (ElementKind::Hex8, true, _) => Ok((Some("BUBBLE"), None)),
        (_, true, _) => Err(format!("incompatible modes on {} elements", kind.name())),
        (ElementKind::Hex20, false, Integration::Reduced) => Ok((Some("TWO"), Some("FULL"))),
        (_, false, Integration::Reduced) => Err(format!("reduced integration of {} elements", kind.name())),
        (_, false, Integration::Full) if quadratic => Ok((Some("THREE"), Some("FULL"))),
        (_, false, Integration::Full) => Ok((Some("TWO"), Some("FULL"))),
    }
}

/// The stiffness of orthotropic engineering constants for Nastran's `MAT9`, as
/// `[[C11, C12, C13], [C22, C23, C33], [G12, G23, G13]]` (the normal block's upper triangle, then
/// the shear moduli for xy, yz and zx in `MAT9`'s order), or `None` when it is not finite.
pub fn orthotropic_stiffness(young: [f64; 3], poisson: [f64; 3], shear: [f64; 3]) -> Option<[[f64; 3]; 3]> {
    let [e1, e2, e3] = young;
    let [n12, n13, n23] = poisson;
    // Compliance of the normal block (symmetric: ν_ij / E_i = ν_ji / E_j).
    let (s11, s22, s33) = (1.0 / e1, 1.0 / e2, 1.0 / e3);
    let (s12, s13, s23) = (-n12 / e1, -n13 / e1, -n23 / e2);
    let det = s11 * (s22 * s33 - s23 * s23) - s12 * (s12 * s33 - s23 * s13) + s13 * (s12 * s23 - s22 * s13);
    let c = [
        [(s22 * s33 - s23 * s23) / det, (s13 * s23 - s12 * s33) / det, (s12 * s23 - s13 * s22) / det],
        [0.0, (s11 * s33 - s13 * s13) / det, (s12 * s13 - s11 * s23) / det],
        [0.0, 0.0, (s11 * s22 - s12 * s12) / det],
    ];
    let [g12, g13, g23] = shear;
    let out = [[c[0][0], c[0][1], c[0][2]], [c[1][1], c[1][2], c[2][2]], [g12, g23, g13]];
    out.iter().flatten().all(|v| v.is_finite()).then_some(out)
}

/// Nastran's component code of a set of degrees of freedom, such as 123 or 156.
fn component(dofs: Dofs) -> i64 {
    dofs.list().into_iter().fold(0, |c, d| c * 10 + i64::from(d))
}

/// The next identification number of a counter.
fn next_id(counter: &mut u32) -> Result<u32> {
    let id = *counter;
    if id > MAX_LABEL {
        return Err(Error::Limit(format!("more than {MAX_LABEL} identification numbers of one kind")));
    }
    *counter += 1;
    Ok(id)
}

/// The id of `value` in `list`, adding it under a new id from `counter` when it is not there.
fn intern<T: PartialEq>(list: &mut Vec<(u32, T)>, counter: &mut u32, value: T) -> Result<u32> {
    if let Some((id, _)) = list.iter().find(|(_, v)| *v == value) {
        return Ok(*id);
    }
    let id = next_id(counter)?;
    list.push((id, value));
    Ok(id)
}

/// A material's entry identification numbers.
#[derive(Debug, Clone, Copy, Default)]
struct MaterialIds {
    /// `MAT1`, for isotropic materials (every element kind).
    isotropic: Option<u32>,
    /// `MAT9`, for orthotropic materials of solids.
    solid: Option<u32>,
    /// `MAT8`, for orthotropic materials of shells.
    shell: Option<u32>,
    /// `TABLES1` of a tabulated plasticity's stress–strain curve.
    table: Option<u32>,
}

/// Which element families use a material.
#[derive(Debug, Clone, Copy, Default)]
struct MaterialUse {
    solid: bool,
    shell: bool,
    line: bool,
}

/// One property entry, and what its elements' connection entries need besides their nodes.
enum Property {
    Solid {
        mid: u32,
        rule: (Option<&'static str>, Option<&'static str>),
    },
    Shell {
        mid: u32,
        thickness: f64,
        offset: Option<f64>,
        mcid: Option<u32>,
    },
    Beam {
        mid: u32,
        profile: BeamProfile,
        orientation: [f64; 3],
    },
    Rod {
        mid: u32,
        area: f64,
    },
    /// Plane stress: a shell with only a membrane material.
    Membrane {
        mid: u32,
        thickness: f64,
        mcid: Option<u32>,
    },
}

/// A constraint set: per node, the prescribed value of each degree of freedom.
type Constraints = BTreeMap<u32, [Option<f64>; 6]>;

/// One step's subcase: the sets it selects.
struct Subcase<'m> {
    step: &'m Step,
    spc: Option<u32>,
    /// The equations of its inclined symmetry planes.
    mpc: Option<u32>,
    load: Option<u32>,
    method: Option<u32>,
}

struct Deck<'m> {
    model: &'m Model,
    r: &'m Resolved,
    sol: u32,
    subcases: Vec<Subcase<'m>>,
    /// Everything the format cannot express.
    problems: BTreeSet<String>,
    next_sid: u32,
    /// Material name → ids.
    materials: BTreeMap<&'m str, MaterialIds>,
    /// The coordinate system tying shell material axes to the basic axes, when one is needed.
    material_csys: Option<u32>,
    /// Properties in id order (id = position + 1), with their section's name.
    properties: Vec<(&'m str, Property)>,
    /// Element → property id.
    element_pid: BTreeMap<u32, u32>,
    /// Rigid elements: (id, coupling name, coupling).
    rigid: Vec<(u32, &'m str, &'m Coupling)>,
    /// Point masses, springs and dashpots: (comment, card).
    discrete: Vec<(String, Card)>,
    /// Node → (dependent degrees of freedom, coupling that makes them dependent).
    dependent: BTreeMap<u32, (Dofs, &'m str)>,
    /// Constraint sets by id.
    spc_sets: Vec<(u32, Constraints)>,
    /// Multipoint constraint sets (inclined symmetry planes) by id.
    mpc_sets: Vec<(u32, Vec<crate::Equation>)>,
    /// Load name → load set id.
    load_sids: BTreeMap<&'m str, u32>,
    /// Combinations: (id, load set ids).
    combinations: Vec<(u32, Vec<u32>)>,
}

/// The solution sequence of a list of steps, or why there is none.
fn solution(steps: &[&Step]) -> std::result::Result<u32, String> {
    if let Some(s) = steps.iter().find(|s| matches!(s.procedure, Procedure::ExplicitDynamic(_))) {
        return Err(format!("step `{}`: explicit dynamic steps (Nastran's linear solutions have none)", s.name));
    }
    if let Some(s) = steps.iter().find(|s| matches!(s.procedure, Procedure::Buckle { .. })) {
        return Err(format!("step `{}`: buckling steps (Nastran's SOL 105; this writer writes SOL 101 and 103 so far)", s.name));
    }
    if let Some(s) = steps.iter().find(|s| matches!(s.procedure, Procedure::NonlinearStatic(_))) {
        return Err(format!("step `{}`: nonlinear static steps (Nastran's SOL 400; this writer writes SOL 101 and 103 so far)", s.name));
    }
    let statics = steps.iter().find(|s| matches!(s.procedure, Procedure::Static {}));
    let modes = steps.iter().find(|s| matches!(s.procedure, Procedure::Frequency { .. }));
    match (statics, modes) {
        (Some(_), None) => Ok(101),
        (None, Some(_)) => Ok(103),
        (Some(s), Some(m)) => Err(format!(
            "one run solves either static steps (SOL 101) or frequency steps (SOL 103), and step `{}` is static while step `{}` extracts frequencies; write them as separate decks (choose the steps to write)",
            s.name, m.name
        )),
        (None, None) => Err("there is no step to write".into()),
    }
}

impl<'m> Deck<'m> {
    fn new(model: &'m Model, r: &'m Resolved, steps: Vec<&'m Step>) -> Deck<'m> {
        let mut problems = BTreeSet::new();
        let sol = solution(&steps).unwrap_or_else(|p| {
            problems.insert(p);
            0
        });
        Deck {
            model,
            r,
            sol,
            subcases: steps.into_iter().map(|step| Subcase { step, spc: None, mpc: None, load: None, method: None }).collect(),
            problems,
            next_sid: 1,
            materials: BTreeMap::new(),
            material_csys: None,
            properties: Vec::new(),
            element_pid: BTreeMap::new(),
            rigid: Vec::new(),
            discrete: Vec::new(),
            dependent: BTreeMap::new(),
            spc_sets: Vec::new(),
            mpc_sets: Vec::new(),
            load_sids: BTreeMap::new(),
            combinations: Vec::new(),
        }
    }

    fn plan(&mut self) -> Result<()> {
        self.plan_materials()?;
        self.plan_sections()?;
        self.plan_couplings()?;
        self.plan_discrete()?;
        self.plan_constraints()?;
        self.plan_loads()?;
        self.problems.extend(crate::thermal_items(self.model));
        for name in self.model.contacts().keys() {
            self.problems.insert(format!("contact `{name}` (Nastran's linear solutions SOL 101 and 103 have no contact)"));
        }
        for name in self.model.ties().keys() {
            self.problems.insert(format!("tie `{name}` (this writer does not write glued contact yet)"));
        }
        for sub in &mut self.subcases {
            if matches!(sub.step.procedure, Procedure::Frequency { .. }) {
                sub.method = Some(next_id(&mut self.next_sid)?);
            }
        }
        crate::unsupported(FORMAT, std::mem::take(&mut self.problems))
    }

    /// Material ids, for the forms the sections need.
    fn plan_materials(&mut self) -> Result<()> {
        let model = self.model;
        let mut uses: BTreeMap<&str, MaterialUse> = BTreeMap::new();
        for s in model.sections().values() {
            let u = uses.entry(s.material.as_str()).or_default();
            match s.kind {
                SectionKind::Solid {} => u.solid = true,
                SectionKind::Shell { .. } | SectionKind::PlaneStress { .. } | SectionKind::PlaneStrain { .. } | SectionKind::Axisymmetric {} => {
                    u.shell = true
                }
                SectionKind::Beam { .. } | SectionKind::Rod { .. } => u.line = true,
            }
        }
        let mut next_mid = 1u32;
        let mut take = |wanted: bool| {
            wanted.then(|| {
                next_mid += 1;
                next_mid - 1
            })
        };
        for (name, u) in uses {
            let ids = match model.material(name).and_then(|m| m.elastic) {
                Some(Elastic::Isotropic { .. }) => MaterialIds { isotropic: take(true), ..MaterialIds::default() },
                Some(Elastic::EngineeringConstants { .. }) => {
                    if u.line {
                        self.problems
                            .insert(format!("material `{name}`: orthotropic constants on beams or rods (their properties take isotropic materials)"));
                    }
                    MaterialIds { solid: take(u.solid), shell: take(u.shell), ..MaterialIds::default() }
                }
                None => MaterialIds::default(),
            };
            let mut ids = ids;
            if let Some(m) = model.material(name) {
                match &m.plastic {
                    Some(Plastic::JohnsonCook { .. }) => {
                        self.problems
                            .insert(format!("material `{name}`: Johnson–Cook plasticity (MATS1 has no rate- and temperature-dependent flow stress)"));
                    }
                    Some(Plastic::Table { .. }) => ids.table = Some(next_id(&mut self.next_sid)?),
                    _ => {}
                }
                if m.damage.is_some() || m.erosion.is_some() {
                    self.problems
                        .insert(format!("material `{name}`: damage and erosion (Nastran's linear and implicit solutions have no element failure)"));
                }
            }
            self.materials.insert(name, ids);
        }
        Ok(())
    }

    /// A new property for `elements`.
    fn add_property(&mut self, section: &'m str, property: Property, elements: impl IntoIterator<Item = u32>) -> Result<()> {
        self.properties.push((section, property));
        let pid = u32::try_from(self.properties.len()).ok().filter(|p| *p <= MAX_LABEL).ok_or_else(|| Error::Limit("too many properties".into()))?;
        for e in elements {
            self.element_pid.insert(e, pid);
        }
        Ok(())
    }

    /// Properties and each element's property.
    fn plan_sections(&mut self) -> Result<()> {
        let (model, r) = (self.model, self.r);
        for (sname, sec) in model.sections() {
            let ids = self.materials.get(sec.material.as_str()).copied().unwrap_or_default();
            let elements = r.sets.get(&sec.set).map(|s| s.elements.as_slice()).unwrap_or_default();
            let kind_of = |e: u32| model.element(e).map(|el| el.kind);
            match sec.kind {
                SectionKind::Solid {} => {
                    let Some(mid) = ids.isotropic.or(ids.solid) else { continue };
                    // One property per integration rule the section's element kinds need.
                    let kinds: BTreeSet<ElementKind> = elements.iter().filter_map(|e| kind_of(*e)).collect();
                    let mut rules = BTreeMap::new();
                    for kind in kinds {
                        match solid_rule(kind, sec.formulation) {
                            Ok(rule) => {
                                rules.entry(rule).or_insert_with(Vec::new).push(kind);
                            }
                            Err(why) => {
                                self.problems.insert(format!("section `{sname}`: {why} (Nastran's linear solids have no such formulation)"));
                            }
                        }
                    }
                    for (rule, kinds) in rules {
                        let members = elements.iter().copied().filter(|e| kind_of(*e).is_some_and(|k| kinds.contains(&k)));
                        self.add_property(sname, Property::Solid { mid, rule }, members)?;
                    }
                }
                SectionKind::Shell { thickness, offset, points } => {
                    if sec.formulation.integration == Integration::Reduced {
                        self.problems.insert(format!("section `{sname}`: reduced integration of shells (Nastran's shells have one formulation)"));
                    }
                    if points < 2 {
                        self.problems.insert(format!(
                            "section `{sname}`: one integration point through the thickness (no bending stiffness; Nastran integrates a linear section exactly)"
                        ));
                    }
                    let (mid, mcid) = match ids {
                        MaterialIds { isotropic: Some(m), .. } => (m, None),
                        MaterialIds { shell: Some(m), .. } => {
                            if !elements.iter().all(|e| model.element(*e).is_some_and(|el| self.in_constant_z(&el.nodes, el.kind))) {
                                self.problems.insert(format!(
                                    "section `{sname}`: orthotropic shells that do not lie in planes of constant Z (the material axes are the global ones)"
                                ));
                            }
                            (m, Some(*self.material_csys.get_or_insert(1)))
                        }
                        _ => continue,
                    };
                    // The nodes lie on the reference surface; Nastran's offset runs from them to
                    // the mid-surface along the element normal.
                    let offset = (offset != 0.0).then_some(-offset * thickness);
                    self.add_property(sname, Property::Shell { mid, thickness, offset, mcid }, elements.iter().copied())?;
                }
                SectionKind::Beam { profile, orientation } => {
                    let Some(mid) = ids.isotropic else { continue };
                    self.add_property(sname, Property::Beam { mid, profile, orientation }, elements.iter().copied())?;
                }
                SectionKind::Rod { area } => {
                    if elements.iter().any(|e| kind_of(*e) == Some(ElementKind::Rod3)) {
                        self.problems.insert(format!("section `{sname}`: three-node trusses (rod3; Nastran's rods have two grids)"));
                    }
                    let Some(mid) = ids.isotropic else { continue };
                    self.add_property(sname, Property::Rod { mid, area }, elements.iter().copied())?;
                }
                SectionKind::PlaneStress { thickness } => {
                    if sec.formulation.integration == Integration::Reduced || sec.formulation.incompatible_modes {
                        self.problems.insert(format!(
                            "section `{sname}`: reduced integration or incompatible modes of plane stress elements (Nastran's membranes have one formulation)"
                        ));
                    }
                    let (mid, mcid) = match ids {
                        MaterialIds { isotropic: Some(m), .. } => (m, None),
                        MaterialIds { shell: Some(m), .. } => (m, Some(*self.material_csys.get_or_insert(1))),
                        _ => continue,
                    };
                    self.add_property(sname, Property::Membrane { mid, thickness, mcid }, elements.iter().copied())?;
                }
                SectionKind::PlaneStrain { .. } | SectionKind::Axisymmetric {} => {
                    self.problems.insert(format!(
                        "section `{sname}`: {} elements (Nastran's plane strain and axisymmetric cards differ between dialects; not written)",
                        sec.kind.name()
                    ));
                }
            }
        }
        Ok(())
    }

    /// Whether a shell element's corners lie in a plane of constant Z.
    fn in_constant_z(&self, nodes: &[u32], kind: ElementKind) -> bool {
        let corners = nodes.get(..kind.corner_count()).unwrap_or_default();
        let mut z = corners.iter().map(|n| self.model.node(*n).map(|n| n.xyz[2]));
        let first = z.next().flatten();
        first.is_some() && z.all(|v| v == first)
    }

    /// Rigid element ids and the degrees of freedom couplings make dependent.
    fn plan_couplings(&mut self) -> Result<()> {
        let (model, r) = (self.model, self.r);
        let mut next = model.elements().keys().next_back().copied().unwrap_or(0);
        let mut overlaps: BTreeMap<(&str, &str), BTreeSet<u32>> = BTreeMap::new();
        for (name, c) in model.couplings() {
            next = next
                .checked_add(1)
                .filter(|n| *n <= MAX_LABEL)
                .ok_or_else(|| Error::Limit(format!("coupling `{name}`: no element identification number left below {}", MAX_LABEL as u64 + 1)))?;
            self.rigid.push((next, name.as_str(), c));
            let dependents: &[u32] = match c.kind {
                CouplingKind::Kinematic => r.sets.get(&c.set).map(|s| s.nodes.as_slice()).unwrap_or_default(),
                CouplingKind::Distributing => r.reference_point_nodes.get(&c.reference_point).map(std::slice::from_ref).unwrap_or_default(),
            };
            for &node in dependents {
                match self.dependent.get(&node).copied() {
                    Some((had, other)) if had.bits() & c.dofs.bits() != 0 => {
                        overlaps.entry((other, name.as_str())).or_default().insert(node);
                    }
                    Some((had, other)) => {
                        self.dependent.insert(node, (Dofs::from_bits(had.bits() | c.dofs.bits()), other));
                    }
                    None => {
                        self.dependent.insert(node, (c.dofs, name));
                    }
                }
            }
        }
        for ((a, b), nodes) in overlaps {
            self.problems.insert(format!(
                "couplings `{a}` and `{b}` both make nodes {} follow another node (Nastran makes a degree of freedom dependent once)",
                labels(&nodes)
            ));
        }
        Ok(())
    }

    /// Point masses, springs and dashpots, numbered after the rigid elements.
    fn plan_discrete(&mut self) -> Result<()> {
        let r = self.r;
        let mut next = self.rigid.last().map(|e| e.0).or_else(|| self.model.elements().keys().next_back().copied()).unwrap_or(0);
        let mut id = || -> Result<u32> {
            next = next
                .checked_add(1)
                .filter(|n| *n <= MAX_LABEL)
                .ok_or_else(|| Error::Limit(format!("no element id left below {} for point masses, springs and dashpots", MAX_LABEL as u64 + 1)))?;
            Ok(next)
        };
        for m in &r.masses {
            let [mx, my, mz] = m.mass;
            let [i11, i22, i33, i12, i13, i23] = m.inertia;
            // M11, M21, M22, M31, M32, M33, then the rotations' rows: M41…M44, M51…M55, M61…M66.
            let card = Card::new("CONM1")
                .int(id()?)
                .int(m.node)
                .blank()
                .reals([mx, 0.0, my, 0.0, 0.0, mz, 0.0, 0.0, 0.0, i11, 0.0, 0.0, 0.0, i12, i22, 0.0, 0.0, 0.0, i13, i23, i33]);
            self.discrete.push((format!("point mass {}", m.name), card));
        }
        for s in &r.springs {
            let what = s.kind();
            let Some(dof) = crate::axis_dof(s.direction) else {
                self.problems.insert(format!("{what} `{}`: a direction along no global axis (CELAS and CDAMP act along grid components)", s.name));
                continue;
            };
            let card = Card::new(if s.dashpot { "CDAMP2" } else { "CELAS2" }).int(id()?).real(s.value).int(s.first).int(dof);
            let card = match s.second {
                Some(b) => card.int(b).int(dof),
                None => card,
            };
            self.discrete.push((format!("{what} {}", s.name), card));
        }
        Ok(())
    }

    /// One constraint set per distinct step constraint, checked against coupled nodes.
    fn plan_constraints(&mut self) -> Result<()> {
        let (model, r) = (self.model, self.r);
        let rotating = model.rotating_nodes();
        // Conflicts, one problem per pair of definitions: the nodes they meet on.
        let mut clashes: BTreeMap<(&str, &str, &str), BTreeSet<u32>> = BTreeMap::new();
        let mut on_dependent: BTreeMap<(&str, &str), BTreeSet<u32>> = BTreeMap::new();
        for sub in &mut self.subcases {
            let step = sub.step;
            let modes = matches!(step.procedure, Procedure::Frequency { .. });
            let mut set: Constraints = BTreeMap::new();
            let mut by: BTreeMap<(u32, u8), &str> = BTreeMap::new();
            for name in &step.boundary_conditions {
                let Some(bc) = model.boundary_condition(name) else { continue };
                let dofs = bc.dofs();
                // Inclined symmetry planes hold no global degree of freedom: they are equations.
                if dofs.is_empty() {
                    continue;
                }
                // In a modes step every constrained degree of freedom is held.
                let values: Vec<(u8, f64)> = dofs.list().into_iter().filter_map(|d| Some((d, if modes { 0.0 } else { bc.value(d)? }))).collect();
                for &node in r.target_nodes(bc.target()) {
                    if let Some((dep, coupling)) = self.dependent.get(&node)
                        && dep.bits() & dofs.bits() != 0
                    {
                        on_dependent.entry((name.as_str(), coupling)).or_default().insert(node);
                    }
                    let slots = set.entry(node).or_default();
                    for &(dof, value) in &values {
                        let Some(slot) = usize::from(dof).checked_sub(1).and_then(|i| slots.get_mut(i)) else { continue };
                        match *slot {
                            Some(v) if v != value => {
                                let other = by.get(&(node, dof)).copied().unwrap_or_default();
                                clashes.entry((step.name.as_str(), other, name.as_str())).or_default().insert(node);
                            }
                            _ => {
                                *slot = Some(value);
                                by.insert((node, dof), name);
                            }
                        }
                    }
                }
            }
            let dependent = &self.dependent;
            let held = |node: u32, dof: u8| match (
                dependent.get(&node),
                set.get(&node).and_then(|v| v.get(usize::from(dof).checked_sub(1)?).copied().flatten()),
            ) {
                (Some((d, _)), _) if d.contains(dof) => crate::Held::Sometimes,
                (_, Some(v)) => crate::Held::Always(v),
                _ => crate::Held::Free,
            };
            let equations = crate::symmetry_equations(model, r, &step.boundary_conditions, &rotating, held, &mut self.problems);
            if !equations.is_empty() {
                sub.mpc = Some(intern(&mut self.mpc_sets, &mut self.next_sid, equations)?);
            }
            if !set.is_empty() {
                sub.spc = Some(intern(&mut self.spc_sets, &mut self.next_sid, set)?);
            }
        }
        for ((step, a, b), nodes) in clashes {
            self.problems.insert(format!(
                "step `{step}`: boundary conditions `{a}` and `{b}` prescribe different values at nodes {} (Nastran takes one value per degree of freedom)",
                labels(&nodes)
            ));
        }
        for ((bc, coupling), nodes) in on_dependent {
            self.problems.insert(format!(
                "boundary condition `{bc}` constrains nodes {} whose motion coupling `{coupling}` prescribes (Nastran cannot constrain a dependent degree of freedom)",
                labels(&nodes)
            ));
        }
        Ok(())
    }

    /// Load set ids for the static steps' loads, and their combinations.
    fn plan_loads(&mut self) -> Result<()> {
        let (model, r) = (self.model, self.r);
        let mut gravity = 0usize;
        for sub in &mut self.subcases {
            if !matches!(sub.step.procedure, Procedure::Static {}) {
                continue;
            }
            let mut sids = Vec::new();
            for name in &sub.step.loads {
                if let Some(sid) = self.load_sids.get(name.as_str()) {
                    sids.push(*sid);
                    continue;
                }
                let Some(load) = model.load(name) else { continue };
                match load {
                    Load::BodyForce { .. } => {
                        self.problems.insert(format!("load `{name}`: a body force (Nastran has no force-per-volume load for these elements)"));
                    }
                    Load::Gravity { set: Some(s), .. } if r.sets.get(s).map(|s| s.elements.len()) != Some(model.elements().len()) => {
                        self.problems.insert(format!("load `{name}`: gravity on part of the model (set `{s}`); Nastran's gravity acts on all of it"));
                    }
                    Load::Gravity { .. } => gravity += 1,
                    Load::Centrifugal { set: Some(s), .. } if r.sets.get(s).map(|s| s.elements.len()) != Some(model.elements().len()) => {
                        self.problems
                            .insert(format!("load `{name}`: a rotation of part of the model (set `{s}`); Nastran's RFORCE acts on all of it"));
                    }
                    Load::Centrifugal { point, axis, .. } => {
                        let c = [
                            point[1] * axis[2] - point[2] * axis[1],
                            point[2] * axis[0] - point[0] * axis[2],
                            point[0] * axis[1] - point[1] * axis[0],
                        ];
                        if c.iter().any(|v| *v != 0.0) {
                            self.problems.insert(format!(
                                "load `{name}`: a rotation about an axis that misses the origin (written as RFORCE about the origin)"
                            ));
                        }
                    }
                    _ => {}
                }
                let sid = next_id(&mut self.next_sid)?;
                self.load_sids.insert(name.as_str(), sid);
                sids.push(sid);
            }
            sub.load = match sids.as_slice() {
                [] => None,
                [one] => Some(*one),
                _ => Some(intern(&mut self.combinations, &mut self.next_sid, sids)?),
            };
        }
        if gravity > MAX_GRAV {
            self.problems.insert(format!("{gravity} gravity loads (a run selects at most {MAX_GRAV})"));
        }
        Ok(())
    }

    fn emit(&self, units: Preset) -> Result<String> {
        let mut out = String::new();
        out.push_str("$ Nastran bulk data written by FieldCraft\n");
        let _ = writeln!(out, "$ units: {}", units.label());
        let _ = writeln!(out, "$ in one consistent unit system: {}; masses are masses (WTMASS 1.0)", crate::units_in_words(units));
        let _ = writeln!(out, "SOL {}\nCEND\nTITLE = FieldCraft model\nECHO = NONE", self.sol);
        for (i, sub) in self.subcases.iter().enumerate() {
            let _ = writeln!(out, "SUBCASE {}", i + 1);
            let title: String = format!("  SUBTITLE = {}", sub.step.name).chars().take(72).collect();
            let _ = writeln!(out, "{title}");
            for (key, sid) in [("SPC", sub.spc), ("MPC", sub.mpc), ("LOAD", sub.load), ("METHOD", sub.method)] {
                if let Some(sid) = sid {
                    let _ = writeln!(out, "  {key} = {sid}");
                }
            }
            for v in &sub.step.output {
                let request = match v {
                    OutputVariable::Displacement => "DISPLACEMENT",
                    OutputVariable::ReactionForce if sub.mpc.is_some() => "SPCFORCES = ALL\n  MPCFORCES",
                    OutputVariable::ReactionForce => "SPCFORCES",
                    OutputVariable::Stress => "STRESS",
                    OutputVariable::Strain => "STRAIN",
                    // Heat transfer steps are refused (thermal_items).
                    OutputVariable::Temperature | OutputVariable::HeatFlux => continue,
                };
                let _ = writeln!(out, "  {request} = ALL");
            }
        }
        out.push_str("BEGIN BULK\n");
        Card::new("PARAM").text("POST").int(-1).write(&mut out)?;
        Card::new("PARAM").text("AUTOSPC").text("YES").write(&mut out)?;
        // The deck's mass unit is consistent with its force, length and time units.
        Card::new("PARAM").text("WTMASS").real(1.0).write(&mut out)?;
        if self.sol == 103 {
            // Consistent (coupled) mass, as the other formats use by default.
            Card::new("PARAM").text("COUPMASS").int(1).write(&mut out)?;
        }
        self.emit_mesh(&mut out)?;
        self.emit_constraints(&mut out)?;
        self.emit_loads(&mut out)?;
        if self.subcases.iter().any(|s| s.method.is_some()) {
            comment(&mut out, "eigenvalue extraction");
        }
        for sub in &self.subcases {
            if let (Some(sid), Procedure::Frequency { modes, min_hz, max_hz }) = (sub.method, &sub.step.procedure) {
                Card::new("EIGRL").int(sid).opt(*min_hz).opt(*max_hz).int(*modes).write(&mut out)?;
            }
        }
        out.push_str("ENDDATA\n");
        Ok(out)
    }

    /// Grids, materials, properties, elements and couplings.
    fn emit_mesh(&self, out: &mut String) -> Result<()> {
        comment(out, "nodes");
        for (label, node) in self.model.nodes() {
            let grid = Card::new("GRID").int(*label).blank().reals(node.xyz);
            // A 2D model's grids hold their out-of-plane translation and rotations.
            let grid = if self.r.is_planar() { grid.blank().int(3456) } else { grid };
            grid.write(out)?;
        }
        if !self.model.reference_points().is_empty() {
            comment(out, "reference points");
        }
        for (name, rp) in self.model.reference_points() {
            comment(out, &format!("reference point {name}: GRID {}", rp.label));
            Card::new("GRID").int(rp.label).blank().reals(rp.xyz).write(out)?;
        }
        if let Some(cid) = self.material_csys {
            comment(out, "material axes of orthotropic shells: the basic axes");
            Card::new("CORD2R").int(cid).blank().reals([0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0]).write(out)?;
        }

        comment(out, "materials");
        for (name, ids) in &self.materials {
            let Some(m) = self.model.material(name) else { continue };
            let rho = m.density;
            let card = match (m.elastic, ids) {
                (Some(Elastic::Isotropic { young, poisson }), MaterialIds { isotropic: Some(mid), .. }) => {
                    Some(Card::new("MAT1").int(*mid).real(young).blank().real(poisson).opt(rho))
                }
                _ => None,
            };
            if let (Some(Elastic::EngineeringConstants { young, poisson, shear }), Some(mid)) = (m.elastic, ids.solid) {
                let [[c11, c12, c13], [c22, c23, c33], [g12, g23, g13]] = orthotropic_stiffness(young, poisson, shear)
                    .ok_or_else(|| Error::Limit(format!("material `{name}`: its stiffness is not a finite number")))?;
                let card = Card::new("MAT9")
                    .int(mid)
                    .reals([c11, c12, c13, 0.0, 0.0, 0.0, c22, c23, 0.0, 0.0, 0.0, c33, 0.0, 0.0, 0.0, g12, 0.0, 0.0, g23, 0.0, g13])
                    .opt(rho);
                write_named(out, &format!("material {name}"), mid, card)?;
            }
            if let (Some(Elastic::EngineeringConstants { young: [e1, e2, _], poisson: [nu12, ..], shear: [g12, g13, g23] }), Some(mid)) =
                (m.elastic, ids.shell)
            {
                write_named(out, &format!("material {name}"), mid, Card::new("MAT8").int(mid).reals([e1, e2, nu12, g12, g13, g23]).opt(rho))?;
            }
            if let (Some(card), Some(mid)) = (card, ids.isotropic) {
                write_named(out, &format!("material {name}"), mid, card)?;
                if let (Some(plastic), Some(Elastic::Isotropic { young, .. })) = (&m.plastic, m.elastic) {
                    plasticity(out, name, mid, young, plastic, ids.table)?;
                }
            }
        }

        comment(out, "properties");
        for (pid, (section, p)) in (1u32..).zip(&self.properties) {
            let card = match p {
                Property::Solid { mid, rule: (network, scheme) } => Card::new("PSOLID").int(pid).int(*mid).blank().opt(*network).blank().opt(*scheme),
                Property::Shell { mid, thickness, .. } => Card::new("PSHELL").int(pid).int(*mid).real(*thickness).int(*mid).blank().int(*mid),
                Property::Beam { mid, profile, .. } => {
                    let (shape, dims) = beam_dimensions(*profile);
                    Card::new("PBEAML").int(pid).int(*mid).blank().text(shape).end_line().reals(dims)
                }
                Property::Rod { mid, area } => Card::new("PROD").int(pid).int(*mid).real(*area),
                Property::Membrane { mid, thickness, .. } => Card::new("PSHELL").int(pid).int(*mid).real(*thickness),
            };
            write_named(out, &format!("section {section}"), pid, card)?;
        }

        comment(out, "elements");
        for (label, el) in self.model.elements() {
            let Some(pid) = self.element_pid.get(label).copied() else { continue };
            let Some((_, property)) = usize::try_from(pid).ok().and_then(|p| self.properties.get(p - 1)) else { continue };
            let card = Card::new(element_card(el.kind)).int(*label).int(pid).ints(grids(el.kind, &el.nodes));
            let card = match property {
                Property::Beam { orientation, .. } => card.reals(self.beam_vector(&el.nodes, *orientation)),
                // CQUAD8 has four corner thicknesses before the material axes and the offset.
                Property::Shell { offset, mcid, .. } if el.kind == ElementKind::Quad8 => card.blank().blank().blank().blank().opt(*mcid).opt(*offset),
                Property::Shell { offset, mcid, .. } => card.opt(*mcid).opt(*offset),
                Property::Membrane { mcid, .. } if el.kind == ElementKind::Quad8 => card.blank().blank().blank().blank().opt(*mcid),
                Property::Membrane { mcid, .. } => card.opt(*mcid),
                Property::Solid { .. } | Property::Rod { .. } => card,
            };
            card.write(out)?;
        }

        if !self.rigid.is_empty() {
            comment(out, "couplings");
        }
        for (eid, name, c) in &self.rigid {
            let rp = self.r.reference_point_nodes.get(&c.reference_point).copied().unwrap_or_default();
            let nodes = self.r.sets.get(&c.set).map(|s| s.nodes.as_slice()).unwrap_or_default();
            let card = match c.kind {
                CouplingKind::Kinematic => Card::new("RBE2").int(*eid).int(rp).int(component(c.dofs)).ints(nodes.iter().copied()),
                // Equal weights on the set's translations.
                CouplingKind::Distributing => {
                    Card::new("RBE3").int(*eid).blank().int(rp).int(component(c.dofs)).real(1.0).int(123).ints(nodes.iter().copied())
                }
            };
            write_named(out, &format!("coupling {name}"), *eid, card)?;
        }
        if !self.discrete.is_empty() {
            comment(out, "point masses, springs and dashpots");
        }
        for (what, card) in &self.discrete {
            comment(out, what);
            card.write(out)?;
        }
        Ok(())
    }

    /// Constraint sets: held degrees of freedom by `SPC1`, prescribed values by `SPC`.
    fn emit_constraints(&self, out: &mut String) -> Result<()> {
        if !self.spc_sets.is_empty() {
            comment(out, "boundary conditions");
        }
        for (sid, set) in &self.spc_sets {
            for sub in self.subcases.iter().filter(|s| s.spc == Some(*sid)) {
                comment(out, &format!("step {}: SPC {sid}", sub.step.name));
            }
            let mut held: BTreeMap<u8, Vec<u32>> = BTreeMap::new();
            let mut prescribed = Vec::new();
            for (node, values) in set {
                let mut zero = 0u8;
                // Nonzero values: one entry per value, with every degree of freedom it applies to.
                let mut by_value: Vec<(f64, u8)> = Vec::new();
                for (bit, v) in values.iter().enumerate().filter_map(|(i, v)| Some((1u8 << i, (*v)?))) {
                    match by_value.iter_mut().find(|(w, _)| *w == v) {
                        _ if v == 0.0 => zero |= bit,
                        Some((_, m)) => *m |= bit,
                        None => by_value.push((v, bit)),
                    }
                }
                if zero != 0 {
                    held.entry(zero).or_default().push(*node);
                }
                prescribed.extend(by_value.into_iter().map(|(v, m)| (*node, m, v)));
            }
            for (mask, nodes) in &held {
                for chunk in nodes.chunks(SPC1_GRIDS) {
                    Card::new("SPC1").int(*sid).int(component(Dofs::from_bits(*mask))).ints(chunk.iter().copied()).write(out)?;
                }
            }
            for (node, mask, value) in prescribed {
                Card::new("SPC").int(*sid).int(node).int(component(Dofs::from_bits(mask))).real(value).write(out)?;
            }
        }
        for (sid, equations) in &self.mpc_sets {
            for sub in self.subcases.iter().filter(|s| s.mpc == Some(*sid)) {
                comment(out, &format!("step {}: inclined symmetry planes, MPC {sid} (the first grid component is dependent)", sub.step.name));
            }
            for e in equations {
                // SID, then (G, C, A) triples: two on the first line (its last field blank), two
                // per continuation after a blank first field.
                let mut card = Card::new("MPC").int(*sid);
                for (i, (dof, c)) in e.terms.iter().enumerate() {
                    if i % 2 == 0 && i > 0 {
                        card = card.blank().blank();
                    }
                    card = card.int(e.node).int(i64::from(*dof)).real(*c);
                }
                card.write(out)?;
            }
        }
        Ok(())
    }

    /// Load sets and their combinations.
    fn emit_loads(&self, out: &mut String) -> Result<()> {
        if !self.load_sids.is_empty() {
            comment(out, "loads");
        }
        for (name, sid) in &self.load_sids {
            let Some(load) = self.model.load(name) else { continue };
            let entry = |card: Card, out: &mut String, first: &mut bool| -> Result<()> {
                if std::mem::take(first) {
                    comment(out, &format!("load {name}: {} {sid}", card.name()));
                }
                card.write(out)
            };
            let mut first = true;
            match load {
                Load::Force { .. } | Load::Traction { .. } | Load::Moment { .. } => {
                    let card = if matches!(load, Load::Moment { .. }) { "MOMENT" } else { "FORCE" };
                    let values = self.r.nodal_loads.get(*name).map(Vec::as_slice).unwrap_or_default();
                    for (node, v) in values.iter().filter(|(_, v)| v.iter().any(|c| *c != 0.0)) {
                        entry(Card::new(card).int(*sid).int(*node).blank().real(1.0).reals(*v), out, &mut first)?;
                    }
                }
                Load::Pressure { .. } if self.r.is_planar() => {
                    let forces = crate::edge_pressure_loads(self.model, self.r, [*name])?;
                    let mut by_node: BTreeMap<u32, [f64; 3]> = BTreeMap::new();
                    for ((node, dof), v) in forces {
                        if let Some(c) = by_node.entry(node).or_default().get_mut(usize::from(dof) - 1) {
                            *c = v;
                        }
                    }
                    for (node, v) in by_node {
                        entry(Card::new("FORCE").int(*sid).int(node).blank().real(1.0).reals(v), out, &mut first)?;
                    }
                }
                Load::Pressure { set, magnitude } => {
                    for f in self.r.sets.get(set).map(|s| s.faces.as_slice()).unwrap_or_default() {
                        if let Some(card) = self.pressure(*sid, *f, *magnitude) {
                            entry(card, out, &mut first)?;
                        }
                    }
                }
                Load::Gravity { acceleration, .. } => entry(Card::new("GRAV").int(*sid).blank().real(1.0).reals(*acceleration), out, &mut first)?,
                Load::Centrifugal { angular_velocity, axis, .. } => {
                    // RFORCE SID G CID A R1 R2 R3: about the basic origin (G blank), A in
                    // revolutions per time unit along the direction R.
                    let a = angular_velocity / (2.0 * std::f64::consts::PI);
                    entry(Card::new("RFORCE").int(*sid).blank().blank().real(a).reals(*axis), out, &mut first)?
                }
                Load::BodyForce { .. } | Load::SurfaceFlux { .. } | Load::Convection { .. } | Load::HeatSource { .. } | Load::Temperature { .. } => {}
            }
        }
        for (sid, loads) in &self.combinations {
            let mut card = Card::new("LOAD").int(*sid).real(1.0);
            for l in loads {
                card = card.real(1.0).int(*l);
            }
            card.write(out)?;
        }
        Ok(())
    }

    /// A `PLOAD4` for one element face, positive pushing into the element.
    fn pressure(&self, sid: u32, f: fieldcraft_model::ElementFace, magnitude: f64) -> Option<Card> {
        let el = self.model.element(f.element)?;
        let card = Card::new("PLOAD4").int(sid).int(f.element);
        if el.kind.dimension() == 2 {
            // Nastran pushes along the element normal; face 1 is the side the normal points to.
            let p = if f.face == 1 { -magnitude } else { magnitude };
            return Some(card.real(p));
        }
        let corners = el.kind.corner_count();
        let local = el.kind.faces().get(usize::from(f.face).checked_sub(1)?)?;
        let face_corners = local.get(..local.iter().take_while(|&&i| i < corners).count())?;
        let node = |i: usize| el.nodes.get(i).copied();
        let card = card.real(magnitude).blank().blank().blank().int(node(*face_corners.first()?)?);
        Some(match (el.kind, face_corners.len()) {
            // Tetrahedra: the corner off the face.
            (ElementKind::Tet4 | ElementKind::Tet10, _) => card.int(node((0..corners).find(|i| !face_corners.contains(i))?)?),
            // Quadrilateral faces: the corner diagonally opposite G1.
            (_, 4) => card.int(node(*face_corners.get(2)?)?),
            // Triangular faces of wedges: G3 stays blank.
            _ => card,
        })
    }

    /// The `CBEAM` orientation vector: the profile's local 1 direction (the section's
    /// orientation) is Nastran's element z axis, so v = axis × orientation.
    fn beam_vector(&self, nodes: &[u32], o: [f64; 3]) -> [f64; 3] {
        let p = |i: usize| nodes.get(i).and_then(|n| self.model.node(*n)).map(|n| n.xyz).unwrap_or_default();
        let (a, b) = (p(0), p(1));
        let t = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
        [t[1] * o[2] - t[2] * o[1], t[2] * o[0] - t[0] * o[2], t[0] * o[1] - t[1] * o[0]]
    }
}

/// `PBEAML`'s shape and dimensions. Nastran's sections have their height along the element y
/// axis and their width (flanges) along z, which [`Deck::beam_vector`] makes the profile's
/// local 1 direction.
fn beam_dimensions(p: BeamProfile) -> (&'static str, Vec<f64>) {
    match p {
        BeamProfile::Rectangle { width, height } => ("BAR", vec![width, height]),
        BeamProfile::Circle { radius } => ("ROD", vec![radius]),
        BeamProfile::Pipe { radius, thickness } if radius - thickness > 0.0 => ("TUBE", vec![radius, radius - thickness]),
        BeamProfile::Pipe { radius, .. } => ("ROD", vec![radius]),
        BeamProfile::I { height, width, web, flange } => ("I", vec![height, width, width, web, flange, flange]),
        BeamProfile::Box { width, height, thickness } => ("BOX", vec![width, height, thickness, thickness]),
    }
}

/// `MATS1` (von Mises yield, `YF` 1) of an isotropic material's plasticity: a bilinear law by its
/// work-hardening slope `H` (stress against plastic strain) and initial yield stress `LIMIT1`,
/// isotropic (`HR` 1) or kinematic (`HR` 2); a hardening table as a `TABLES1` stress–strain
/// curve through the origin, with the total strain σ/E + ε_p of each point.
fn plasticity(out: &mut String, name: &str, mid: u32, young: f64, plastic: &Plastic, table: Option<u32>) -> Result<()> {
    let what = format!("plasticity of {name}");
    match (plastic, table) {
        (Plastic::Bilinear { yield_stress, hardening_modulus, hardening }, _) => {
            let hr = if *hardening == Hardening::Kinematic { 2 } else { 1 };
            let card = Card::new("MATS1").int(mid).blank().text("PLASTIC").real(*hardening_modulus).int(1).int(hr).real(*yield_stress);
            write_named(out, &what, mid, card)
        }
        (Plastic::Table { points }, Some(tid)) => {
            write_named(out, &what, mid, Card::new("MATS1").int(mid).int(tid).text("PLASTIC").blank().int(1).int(1))?;
            let curve = std::iter::once([0.0, 0.0]).chain(points.iter().map(|p| [p.stress / young + p.strain, p.stress])).flatten();
            write_named(out, &format!("stress-strain curve of {name}"), tid, Card::new("TABLES1").int(tid).end_line().reals(curve).text("ENDT"))
        }
        (Plastic::Table { .. }, None) | (Plastic::JohnsonCook { .. }, _) => {
            Err(Error::Model(fieldcraft_model::Error::Invalid(format!("material `{name}`: its plasticity was not planned for MATS1"))))
        }
    }
}

/// An entry preceded by a comment naming what it holds (`material Steel: MAT1 1`).
fn write_named(out: &mut String, what: &str, id: u32, card: Card) -> Result<()> {
    comment(out, &format!("{what}: {} {id}", card.name()));
    card.write(out)
}

/// The first few labels of a list, and how many more there are.
fn labels(nodes: &BTreeSet<u32>) -> String {
    const SHOWN: usize = 5;
    let shown: Vec<String> = nodes.iter().take(SHOWN).map(u32::to_string).collect();
    match nodes.len().saturating_sub(SHOWN) {
        0 => shown.join(", "),
        more => format!("{}, and {more} more", shown.join(", ")),
    }
}

/// A comment line, at most 80 columns.
fn comment(out: &mut String, text: &str) {
    out.extend(format!("$ {text}").chars().take(80));
    out.push('\n');
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn orthotropic_stiffness_of_isotropic_constants_is_lames() {
        let (e, nu) = (210_000.0, 0.3);
        let g = e / (2.0 * (1.0 + nu));
        let c = orthotropic_stiffness([e; 3], [nu; 3], [g; 3]).unwrap();
        let k = e / ((1.0 + nu) * (1.0 - 2.0 * nu));
        let close = |a: f64, b: f64| (a - b).abs() < 1e-9 * b.abs();
        assert!(close(c[0][0], k * (1.0 - nu)) && close(c[1][0], k * (1.0 - nu)) && close(c[1][2], k * (1.0 - nu)));
        assert!(close(c[0][1], k * nu) && close(c[0][2], k * nu) && close(c[1][1], k * nu));
        assert_eq!(c[2], [g; 3]);
    }

    #[test]
    fn hex20_mid_side_nodes_swap_blocks() {
        let nodes: Vec<u32> = (1..=20).collect();
        let g: Vec<u32> = grids(ElementKind::Hex20, &nodes).collect();
        assert_eq!(&g[..12], &nodes[..12]);
        assert_eq!(&g[12..16], &[17, 18, 19, 20]);
        assert_eq!(&g[16..], &[13, 14, 15, 16]);
        for k in ElementKind::ALL {
            let mut sorted = order(k).to_vec();
            sorted.sort_unstable();
            assert_eq!(sorted, (0..k.node_count()).collect::<Vec<_>>(), "{k:?} is a permutation");
        }
    }

    #[test]
    fn formulations_map_or_are_named() {
        let f = |integration, incompatible_modes, hybrid| Formulation { integration, incompatible_modes, hybrid };
        assert_eq!(solid_rule(ElementKind::Hex8, f(Integration::Full, false, false)), Ok((Some("TWO"), Some("FULL"))));
        assert_eq!(solid_rule(ElementKind::Hex8, f(Integration::Full, true, false)), Ok((Some("BUBBLE"), None)));
        assert_eq!(solid_rule(ElementKind::Hex20, f(Integration::Reduced, false, false)), Ok((Some("TWO"), Some("FULL"))));
        assert_eq!(solid_rule(ElementKind::Tet10, f(Integration::Full, false, false)), Ok((Some("THREE"), Some("FULL"))));
        assert!(solid_rule(ElementKind::Hex8, f(Integration::Reduced, false, false)).is_err());
        assert!(solid_rule(ElementKind::Tet4, f(Integration::Full, false, true)).is_err());
    }
}
