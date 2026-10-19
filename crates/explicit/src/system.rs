//! From the resolved model to the explicit system: elements (solids and trusses) with their
//! materials and stable lengths, lumped masses with the point masses, springs and dashpots
//! ([`crate::discrete`]), prescribed velocities, rigid walls, initial velocities and the load
//! vector of an explicit dynamic step.

use fieldcraft_fem::constitutive::MaterialModel;
use fieldcraft_fem::material::Elasticity;
use fieldcraft_fem::plasticity::J2;
use fieldcraft_fem::uniform::{Hex8, Tet4, modulus_bound};
use fieldcraft_fem::{DofKind, DofLayout, DofSet, face};
use fieldcraft_model::{
    BoundaryCondition, Elastic, ElementKind, ExplicitDynamic, Friction, HourglassKind, Integration, Load, Model, Resolved, SectionKind, Step,
    StressRate, Target,
};

use crate::contact::Interface;
use crate::discrete::{Discrete, Rod};
use crate::elements::{HexScales, add_at, frequency_bound, hex_scales};
use crate::{Error, Limits, Result};

/// Degrees of freedom per node: the three translations, node-major (see [`System`]).
pub const NODE_DOFS: usize = 3;

/// Most nodes whose histories a step records.
pub const MAX_HISTORY_NODES: usize = 1000;

/// A material as the elements use it.
#[derive(Debug, Clone)]
pub struct Material {
    pub name: String,
    pub elasticity: Elasticity,
    pub density: f64,
    /// Dilatational wave speed √(M/ρ), with M the largest normal stiffness (λ + 2μ when
    /// isotropic).
    pub wave_speed: f64,
    /// The bound of [`fieldcraft_fem::uniform::modulus_bound`].
    pub modulus: f64,
    /// The largest normal stiffness M (λ + 2μ when isotropic): c = √(M/ρ).
    pub stiffest: f64,
    /// The hypoelastic(–plastic) update of an isotropic material; `None` for orthotropic
    /// elasticity, updated with [`Material::elasticity`].
    pub j2: Option<J2>,
    /// The finite-strain material point of a hyperelastic material (its stress from the
    /// deformation gradient; `None` for the rate-form materials).
    pub hyperelastic: Option<MaterialModel>,
    /// The equivalent plastic strain at which elements erode.
    pub erosion: Option<f64>,
    /// The (isotropic) expansion coefficient and reference temperature, when the material
    /// expands.
    pub expansion: Option<(f64, f64)>,
}

/// A one-point hexahedron.
#[derive(Debug, Clone)]
pub struct Hex {
    pub label: u32,
    /// Position in [`System::elements`].
    pub index: usize,
    pub nodes: [usize; 8],
    pub material: usize,
    pub kinematics: Hex8,
    pub gamma: [[f64; 8]; 4],
    /// Row-sum lumped masses of the element's nodes.
    pub mass: [f64; 8],
    /// Size, hourglass coefficient and frequency bound of the undeformed element.
    pub rest: HexScales,
    /// The element's stable length: its critical increment without bulk viscosity is
    /// `length / wave_speed`.
    pub length: f64,
    /// The rate of its thermal strain, α ΔT / end time (every normal direction).
    pub thermal_rate: f64,
    /// Its density's factor from mass scaling (1 without).
    pub mass_factor: f64,
}

/// A nodal-pressure tetrahedron.
#[derive(Debug, Clone)]
pub struct Tet {
    pub label: u32,
    pub index: usize,
    pub nodes: [usize; 4],
    pub material: usize,
    pub kinematics: Tet4,
    pub mass: [f64; 4],
    /// The frequency bound ω² of the undeformed element.
    pub omega2: f64,
    pub length: f64,
    /// The rate of its thermal strain, α ΔT / end time (every normal direction).
    pub thermal_rate: f64,
    /// Its density's factor from mass scaling (1 without).
    pub mass_factor: f64,
}

/// A fixed rigid plane.
#[derive(Debug, Clone)]
pub struct Wall {
    pub point: [f64; 3],
    /// Unit normal, towards the side the nodes stay on.
    pub normal: [f64; 3],
    pub nodes: Vec<usize>,
    pub friction: Friction,
}

/// One velocity constraint of the symmetry planes not normal to an axis at a node: c · v = 0
/// over the node's translations, `q` its part on the components no prescribed velocity holds,
/// of unit length and perpendicular to the node's other rows' (Gram–Schmidt over the planes that
/// meet at the node). [`SymmetryRow::project`] removes q (c · v) from the velocity: with the same
/// lumped mass on a node's three translations this is the mass-weighted closest velocity that
/// satisfies the planes and the prescribed values (the constrained update of a lumped-mass
/// central difference scheme), exactly, in one pass over the rows.
#[derive(Debug, Clone)]
pub struct SymmetryRow {
    pub node: usize,
    pub q: [f64; 3],
    pub c: [f64; 3],
}

impl SymmetryRow {
    pub fn project(&self, v: &mut [f64]) {
        let at = |k: usize| NODE_DOFS * self.node + k;
        let along: f64 = (0..3).map(|k| self.c[k] * v.get(at(k)).copied().unwrap_or(0.0)).sum();
        for k in 0..3 {
            if let Some(x) = v.get_mut(at(k)) {
                *x -= along * self.q[k];
            }
        }
    }
}

/// A free part below this fraction of its row's size is round-off: the row repeats the node's
/// earlier ones, or the prescribed velocities hold it.
const REPEATED: f64 = 1e-10;

/// Prescribed velocities satisfy a constraint row, and a node lies on a rigid wall's side, within
/// this fraction of the size of the terms compared: coordinates and velocities are exact to a few
/// units in the last place.
const ROUND_OFF: f64 = 1e-9; // unit-factors: relative to the terms compared

/// The discrete system of an explicit dynamic step.
///
/// Its elements are solids, so every node carries the three translations and the layout is
/// node-major: degree of freedom 3 n + axis, which the integrator's per-node loops rely on
/// (checked when the system is built).
#[derive(Debug, Clone)]
pub struct System {
    /// Node labels with degrees of freedom (those of the elements), sorted.
    pub labels: Vec<u32>,
    pub layout: DofLayout,
    pub xyz: Vec<[f64; 3]>,
    /// Element labels, sorted.
    pub elements: Vec<u32>,
    pub materials: Vec<Material>,
    pub hexes: Vec<Hex>,
    pub tets: Vec<Tet>,
    pub rods: Vec<Rod>,
    /// Row-sum lumped mass per node.
    pub mass: Vec<f64>,
    /// Point masses, springs and dashpots.
    pub discrete: Discrete,
    /// The springs' and dashpots' stable increment on the lumped mass, and a link that sets it,
    /// when there are any.
    pub discrete_increment: Option<(f64, String)>,
    /// What mass scaling added: per node, and per section (in name order).
    pub added_mass: Vec<f64>,
    pub mass_scaling: Option<crate::result::MassScalingSummary>,
    /// Prescribed velocity per degree of freedom (`None` = free).
    pub prescribed: Vec<Option<f64>>,
    pub walls: Vec<Wall>,
    /// The symmetry planes not normal to an axis, per node.
    pub planes: Vec<SymmetryRow>,
    /// The model's contacts ([`crate::contact`]).
    pub contacts: Vec<Interface>,
    /// Velocity per degree of freedom at the start (prescribed velocities included).
    pub initial_velocity: Vec<f64>,
    /// Constant external force per degree of freedom.
    pub external: Vec<f64>,
    /// Nodes whose histories are recorded.
    pub history_nodes: Vec<usize>,
    pub hourglass: HourglassKind,
    pub hourglass_coefficient: f64,
    /// Large deformation (updated Lagrangian) or the undeformed mesh's gradients.
    pub nonlinear: bool,
    pub stress_rate: StressRate,
    pub warnings: Vec<String>,
}

fn invalid(m: impl Into<String>) -> Error {
    Error::Invalid(m.into())
}

/// A material's constants for the explicit elements.
fn material(model: &Model, name: &str) -> Result<Material> {
    let m = model.material(name).ok_or_else(|| invalid(format!("material `{name}` is missing")))?;
    let elastic = m.elastic.as_ref().ok_or_else(|| invalid(format!("material `{name}` has no elastic constants")))?;
    let density = m
        .density
        .filter(|d| *d > 0.0 && d.is_finite())
        .ok_or_else(|| invalid(format!("material `{name}` has no density, and explicit dynamics needs one")))?;
    let d = fieldcraft_fem::elasticity(elastic).map_err(|e| invalid(format!("material `{name}`: {e}")))?;
    let stiffest = d[0][0].max(d[1][1]).max(d[2][2]);
    // The model keeps plasticity, damage and erosion to isotropic materials.
    let j2 = match *elastic {
        Elastic::Isotropic { .. } => {
            Some(J2::from_model(elastic, m.plastic.as_ref(), m.damage.as_ref()).map_err(|e| invalid(format!("material `{name}`: {e}")))?)
        }
        Elastic::EngineeringConstants { .. } => None,
    };
    let expansion = match m.thermal.and_then(|t| t.expansion.map(|a| (a.axes(), t.reference_temperature))) {
        None => None,
        Some(([a, b, c], t_ref)) if a == b && b == c => Some((a, t_ref)),
        Some(_) => return Err(Error::NotSupported(format!("orthotropic thermal expansion (material `{name}`) in explicit dynamics"))),
    };
    Ok(Material {
        name: name.to_string(),
        expansion,
        elasticity: d,
        density,
        wave_speed: (stiffest / density).sqrt(),
        modulus: modulus_bound(elastic, &d),
        stiffest,
        j2,
        erosion: m.erosion.map(|e| e.plastic_strain),
        hyperelastic: m
            .hyperelastic
            .as_ref()
            .map(|h| MaterialModel::from_model(elastic, Some(h)))
            .transpose()
            .map_err(|e| invalid(format!("material `{name}`: {e}")))?,
    })
}

/// The section's elements of a kind the explicit solver has, or why not.
fn check_formulation(kind: ElementKind, section_name: &str, f: &fieldcraft_model::Formulation) -> Result<()> {
    let refuse = |why: &str| Err(Error::NotSupported(format!("section `{section_name}`: {why}")));
    if f.hybrid {
        return refuse("hybrid elements in explicit dynamics");
    }
    match kind {
        ElementKind::Hex8 if f.incompatible_modes || f.integration == Integration::Full => refuse(
            "fully integrated or incompatible-mode hex8 in explicit dynamics (explicit hex8 elements are one-point with hourglass control: set the section's integration to reduced)",
        ),
        ElementKind::Hex8 => Ok(()),
        ElementKind::Tet4 if f.integration == Integration::Reduced => refuse("reduced integration of tet4 elements"),
        ElementKind::Tet4 => Ok(()),
        ElementKind::Rod2 | ElementKind::Rod3 => Ok(()),
        other => refuse(&format!("{} elements in explicit dynamics (it has one-point hex8 and nodal-pressure tet4 so far)", other.name())),
    }
}

impl System {
    pub fn new(model: &Model, resolved: &Resolved, step: &Step, settings: &ExplicitDynamic, limits: &Limits) -> Result<System> {
        let max_dofs = limits.max_dofs;
        if let Some(name) = model.couplings().keys().next() {
            return Err(Error::NotSupported(format!("couplings (coupling `{name}`) in explicit dynamics")));
        }
        if let Some(name) = model.ties().keys().next() {
            return Err(Error::NotSupported(format!("ties (tie `{name}`) in explicit dynamics")));
        }

        // Elements by label, with their section's material.
        let mut materials: Vec<Material> = Vec::new();
        let mut chosen: Vec<(u32, ElementKind, Vec<u32>, usize, Option<f64>)> = Vec::new();
        for (name, section) in model.sections() {
            let area = match section.kind {
                SectionKind::Solid {} => None,
                SectionKind::Rod { area } => Some(area),
                _ => {
                    return Err(Error::NotSupported(format!(
                        "{} sections (section `{name}`); explicit dynamics has solid and truss elements so far",
                        section.kind.name()
                    )));
                }
            };
            let mi = match materials.iter().position(|m| m.name == section.material) {
                Some(i) => i,
                None => {
                    materials.push(material(model, &section.material)?);
                    materials.len() - 1
                }
            };
            let set = resolved.sets.get(&section.set).ok_or_else(|| invalid(format!("set `{}` is not resolved", section.set)))?;
            for &label in &set.elements {
                let el = model.element(label).ok_or_else(|| invalid(format!("element {label} is missing")))?;
                if area.is_some() != matches!(el.kind, ElementKind::Rod2 | ElementKind::Rod3) {
                    return Err(invalid(format!(
                        "section `{name}`: a {} section cannot hold {} element {label}",
                        section.kind.name(),
                        el.kind.name()
                    )));
                }
                check_formulation(el.kind, name, &section.formulation)?;
                if el.kind == ElementKind::Tet4 && materials.get(mi).is_some_and(|m| m.j2.is_none()) {
                    return Err(Error::NotSupported(format!(
                        "nodal-pressure tet4 elements of the orthotropic material `{}` (section `{name}`)",
                        section.material
                    )));
                }
                chosen.push((label, el.kind, el.nodes.clone(), mi, area));
            }
        }
        chosen.sort_by_key(|e| e.0);
        let mut labels: Vec<u32> = chosen.iter().flat_map(|e| e.2.iter().copied()).collect();
        // Mesh nodes only discrete elements use (reference points are refused by Discrete::new).
        labels.extend(resolved.discrete_labels().filter(|l| model.node(*l).is_some()));
        labels.sort_unstable();
        labels.dedup();
        let layout = DofLayout::uniform(labels.iter().copied(), DofSet::TRANSLATIONS);
        if layout.dofs() > max_dofs {
            return Err(Error::Limit(format!("the model has {} degrees of freedom; the limit is {max_dofs}", layout.dofs())));
        }
        if layout.stride() != Some(NODE_DOFS) {
            return Err(Error::Internal("explicit dynamics needs three translations per node".into()));
        }
        if labels.is_empty() {
            return Err(invalid("the model has no elements or point masses to analyse"));
        }
        let xyz: Vec<[f64; 3]> =
            labels.iter().map(|l| model.node(*l).map(|n| n.xyz).ok_or_else(|| invalid(format!("node {l} is missing")))).collect::<Result<_>>()?;
        let (n, dofs) = (labels.len(), layout.dofs());
        let mut sys = System {
            layout,
            elements: chosen.iter().map(|e| e.0).collect(),
            materials,
            hexes: Vec::new(),
            tets: Vec::new(),
            rods: Vec::new(),
            mass: vec![0.0; n],
            added_mass: Vec::new(),
            mass_scaling: None,
            discrete: Discrete::default(),
            discrete_increment: None,
            prescribed: vec![None; dofs],
            walls: Vec::new(),
            planes: Vec::new(),
            contacts: Vec::new(),
            initial_velocity: vec![0.0; dofs],
            external: vec![0.0; dofs],
            history_nodes: Vec::new(),
            hourglass: settings.hourglass.kind,
            hourglass_coefficient: settings.hourglass.coefficient,
            nonlinear: settings.nonlinear_geometry,
            stress_rate: settings.stress_rate,
            warnings: resolved.warnings.iter().map(ToString::to_string).collect(),
            labels,
            xyz,
        };
        // Contact first: nodes that start inside a surface are moved before the elements are
        // built on the coordinates.
        sys.contacts = crate::contact::interfaces(model, resolved, &mut sys.xyz, &sys.labels, &sys.elements, limits.max_contact_faces)?;
        for c in &sys.contacts {
            if let Some(deepest) = c.adjusted.iter().max_by(|a, b| a.depth.total_cmp(&b.depth)) {
                sys.warnings.push(format!(
                    "contact `{}`: {} nodes started inside the other surface and were moved onto it (the deepest, node {}, by {:.3e})",
                    c.name,
                    c.adjusted.len(),
                    deepest.node,
                    deepest.depth
                ));
            }
        }
        for (index, (label, kind, nodes, mi, area)) in chosen.into_iter().enumerate() {
            let idx: Vec<usize> = nodes.iter().map(|l| sys.node(*l).ok_or_else(|| Error::Internal(format!("node {l}")))).collect::<Result<_>>()?;
            let mat = sys.materials.get(mi).ok_or_else(|| Error::Internal("material index".into()))?;
            let (rho, c) = (mat.density, mat.wave_speed);
            let at = |e| Error::element(label, e);
            match kind {
                ElementKind::Hex8 => {
                    let nodes: [usize; 8] = idx.try_into().map_err(|_| invalid(format!("element {label} does not have 8 nodes")))?;
                    let x = nodes.map(|i| sys.xyz.get(i).copied().unwrap_or_default());
                    let k = Hex8::hex8(&x).map_err(at)?;
                    let gamma = k.hourglass_vectors(&x);
                    let mass = k.nodal_volume.map(|v| rho * v);
                    let rest = hex_scales(&k, &gamma, &mass, mat, rho, settings.hourglass);
                    sys.lump(&nodes, &mass);
                    let length = 2.0 * c / rest.omega2.sqrt();
                    sys.hexes.push(Hex {
                        label,
                        index,
                        nodes,
                        material: mi,
                        kinematics: k,
                        gamma,
                        mass,
                        rest,
                        length,
                        thermal_rate: 0.0,
                        mass_factor: 1.0,
                    });
                }
                ElementKind::Tet4 => {
                    let nodes: [usize; 4] = idx.try_into().map_err(|_| invalid(format!("element {label} does not have 4 nodes")))?;
                    let x = nodes.map(|i| sys.xyz.get(i).copied().unwrap_or_default());
                    let k = Tet4::tet4(&x).map_err(at)?;
                    let mass = k.nodal_volume.map(|v| rho * v);
                    let omega2 = frequency_bound(&k, k.gradient_norm(), mat.modulus, &mass);
                    sys.lump(&nodes, &mass);
                    sys.tets.push(Tet {
                        label,
                        index,
                        nodes,
                        material: mi,
                        kinematics: k,
                        mass,
                        omega2,
                        length: 2.0 * c / omega2.sqrt(),
                        thermal_rate: 0.0,
                        mass_factor: 1.0,
                    });
                }
                ElementKind::Rod2 | ElementKind::Rod3 => {
                    let area = area.ok_or_else(|| Error::Internal("truss area".into()))?;
                    let source = model.material(&mat.name).ok_or_else(|| Error::Internal("truss material".into()))?;
                    if source.plastic.is_some() || source.damage.is_some() {
                        return Err(Error::NotSupported(format!(
                            "plastic or damaging materials on trusses (material `{}`, element {label}): explicit trusses are elastic so far",
                            mat.name
                        )));
                    }
                    if kind == ElementKind::Rod3 && settings.nonlinear_geometry {
                        return Err(Error::NotSupported(format!(
                            "three-node trusses with nonlinear geometry (element {label}); set the step's nonlinear_geometry to false"
                        )));
                    }
                    let young =
                        fieldcraft_fem::truss::young_modulus(source.elastic.as_ref().ok_or_else(|| Error::Internal("truss material".into()))?)
                            .map_err(|e| invalid(format!("material `{}` of element {label}: {e}", mat.name)))?;
                    let x: Vec<[f64; 3]> = idx.iter().map(|i| sys.xyz.get(*i).copied().unwrap_or_default()).collect();
                    let rod = Rod::new(label, index, kind, idx.clone(), &x, mi, crate::discrete::RodSection { young, density: rho, area })?;
                    let mass: Vec<f64> = rod.nodal_volume.iter().map(|v| rho * v).collect();
                    sys.lump(&idx, &mass);
                    sys.rods.push(rod);
                }
                other => return Err(Error::Internal(format!("{} element kind", other.name()))),
            }
        }
        if let Some(scaling) = &settings.mass_scaling {
            sys.scale_mass(model, resolved, scaling, settings)?;
        }
        sys.discrete = Discrete::new(resolved, sys.dofs(), &sys.xyz, |l| sys.node(l))?;
        sys.discrete_increment = sys.discrete.critical_increment(&sys.dof_mass(), sys.nonlinear);
        sys.prescribe(model, resolved, step, settings.end_time)?;
        sys.load(model, resolved, step)?;
        sys.temperatures(model, resolved, step, settings.end_time)?;
        sys.conditions(resolved, settings)?;
        sys.check_hyperelastic()?;
        Ok(sys)
    }

    /// Mass scaling to a target increment (see [`fieldcraft_model::MassScaling`]): each solid
    /// element (of the set) whose stable increment times the step's scale is below the target
    /// gets its density times (target / increment)², its lumped masses and frequency bound with
    /// it; the added mass is kept per node and per section.
    fn scale_mass(&mut self, model: &Model, resolved: &Resolved, scaling: &fieldcraft_model::MassScaling, settings: &ExplicitDynamic) -> Result<()> {
        let chosen: Option<&[u32]> = match &scaling.set {
            Some(s) => Some(resolved.sets.get(s).map(|r| r.elements.as_slice()).ok_or_else(|| invalid(format!("set `{s}` is not resolved")))?),
            None => None,
        };
        let wanted = |label: u32| chosen.is_none_or(|l| l.binary_search(&label).is_ok());
        let factor = |length: f64, material: usize, materials: &[Material]| -> f64 {
            let c = materials.get(material).map_or(f64::INFINITY, |m| m.wave_speed);
            let dt = settings.dt_scale * length / c;
            if dt > 0.0 && dt < scaling.target_increment { (scaling.target_increment / dt).powi(2) } else { 1.0 }
        };
        let hourglass = fieldcraft_model::Hourglass { kind: settings.hourglass.kind, coefficient: settings.hourglass.coefficient };
        let mut added = vec![0.0; self.labels.len()];
        let mut by_element: Vec<(u32, f64)> = Vec::new();
        for h in self.hexes.iter_mut().filter(|h| wanted(h.label)) {
            let f = factor(h.length, h.material, &self.materials);
            if f == 1.0 {
                continue;
            }
            let extra = h.mass.map(|m| m * (f - 1.0));
            for (n, x) in h.nodes.iter().zip(extra) {
                if let Some(a) = added.get_mut(*n) {
                    *a += x;
                }
            }
            h.mass = h.mass.map(|m| m * f);
            h.mass_factor = f;
            if let Some(mat) = self.materials.get(h.material) {
                h.rest = hex_scales(&h.kinematics, &h.gamma, &h.mass, mat, mat.density * f, hourglass);
            }
            by_element.push((h.label, extra.iter().sum()));
        }
        for t in self.tets.iter_mut().filter(|t| wanted(t.label)) {
            let f = factor(t.length, t.material, &self.materials);
            if f == 1.0 {
                continue;
            }
            let extra = t.mass.map(|m| m * (f - 1.0));
            for (n, x) in t.nodes.iter().zip(extra) {
                if let Some(a) = added.get_mut(*n) {
                    *a += x;
                }
            }
            t.mass = t.mass.map(|m| m * f);
            t.omega2 /= f;
            t.mass_factor = f;
            by_element.push((t.label, extra.iter().sum()));
        }
        for (m, a) in self.mass.iter_mut().zip(&added) {
            *m += a;
        }
        let total_added: f64 = added.iter().sum();
        let total: f64 = self.mass.iter().sum();
        let mut sections: Vec<(String, f64)> = Vec::new();
        for (name, section) in model.sections() {
            let elements = resolved.sets.get(&section.set).map(|r| r.elements.as_slice()).unwrap_or_default();
            let mass: f64 = by_element.iter().filter(|(l, _)| elements.binary_search(l).is_ok()).map(|(_, m)| m).sum();
            if mass > 0.0 {
                sections.push((name.clone(), mass));
            }
        }
        if total_added > 0.0 {
            self.warnings.push(format!(
                "mass scaling to an increment of {:.3e} added {total_added:.4e} of mass ({:.2}% of the model's) to {} elements",
                scaling.target_increment,
                100.0 * total_added / total.max(f64::MIN_POSITIVE),
                by_element.len()
            ));
        }
        self.mass_scaling = Some(crate::result::MassScalingSummary {
            target_increment: scaling.target_increment,
            elements: by_element.len(),
            added_mass: total_added,
            fraction: total_added / total.max(f64::MIN_POSITIVE),
            sections,
        });
        self.added_mass = added;
        Ok(())
    }

    /// Refuse what hyperelastic elements cannot take yet: small-displacement steps (the
    /// material is a finite-strain law) and thermal strain.
    fn check_hyperelastic(&self) -> Result<()> {
        let solids =
            self.hexes.iter().map(|h| (h.label, h.material, h.thermal_rate)).chain(self.tets.iter().map(|t| (t.label, t.material, t.thermal_rate)));
        for (label, material, thermal) in solids {
            let Some(m) = self.materials.get(material).filter(|m| m.hyperelastic.is_some()) else { continue };
            if !self.nonlinear {
                return Err(Error::NotSupported(format!(
                    "hyperelastic material `{}` (element {label}) without nonlinear geometry; hyperelastic laws are finite-strain, set the step's nonlinear_geometry to true",
                    m.name
                )));
            }
            if thermal != 0.0 {
                return Err(Error::NotSupported(format!(
                    "thermal strain in hyperelastic material `{}` (element {label}) in explicit dynamics",
                    m.name
                )));
            }
        }
        Ok(())
    }

    /// Add an element's row-sum lumped masses.
    fn lump(&mut self, nodes: &[usize], mass: &[f64]) {
        for (node, v) in nodes.iter().zip(mass) {
            if let Some(m) = self.mass.get_mut(*node) {
                *m += v;
            }
        }
    }

    pub fn node(&self, label: u32) -> Option<usize> {
        self.labels.binary_search(&label).ok()
    }

    pub fn dofs(&self) -> usize {
        self.layout.dofs()
    }

    /// The lumped mass per degree of freedom: the elements' and the point masses'.
    pub fn dof_mass(&self) -> Vec<f64> {
        (0..self.dofs())
            .map(|d| self.mass.get(crate::discrete::node_direction(d).0).copied().unwrap_or(0.0) + self.discrete.mass.get(d).copied().unwrap_or(0.0))
            .collect()
    }

    /// Prescribed velocities from the step's boundary conditions: zero for held degrees of
    /// freedom, value / end time for a displacement reached over the step.
    fn prescribe(&mut self, model: &Model, resolved: &Resolved, step: &Step, end_time: f64) -> Result<()> {
        let mut normals: std::collections::BTreeMap<usize, Vec<[f64; 3]>> = std::collections::BTreeMap::new();
        for name in &step.boundary_conditions {
            let bc = model.boundary_condition(name).ok_or_else(|| invalid(format!("boundary condition `{name}` is missing")))?;
            if let Target::ReferencePoint(rp) = bc.target() {
                return Err(Error::NotSupported(format!("boundary conditions on reference points (`{name}` on `{rp}`)")));
            }
            if bc.is_thermal() {
                return Err(Error::NotSupported(format!("temperature boundary conditions (`{name}`) in explicit dynamics")));
            }
            if let BoundaryCondition::Symmetry { plane, .. } = bc
                && let Some(normal) = plane.inclined()
            {
                for node in resolved.target_nodes(bc.target()).iter().filter_map(|l| self.node(*l)) {
                    normals.entry(node).or_default().push(normal);
                }
                continue;
            }
            for &label in resolved.target_nodes(bc.target()) {
                let Some(node) = self.node(label) else { continue };
                for kind in self.layout.kinds(node).iter() {
                    let dof = kind.number();
                    let Some(value) = bc.value(dof) else { continue };
                    let velocity = value / end_time;
                    if !velocity.is_finite() {
                        return Err(invalid(format!("boundary condition `{name}`: the velocity {value} / {end_time} is not finite")));
                    }
                    let index = self.layout.dof(node, kind).ok_or_else(|| Error::Internal("dof index".into()))?;
                    let slot = self.prescribed.get_mut(index).ok_or_else(|| Error::Internal("dof index".into()))?;
                    match *slot {
                        Some(old) if old != velocity => {
                            return Err(invalid(format!(
                                "node {label} direction {dof} moves at both {old} and {velocity} (boundary condition `{name}` and another)"
                            )));
                        }
                        _ => *slot = Some(velocity),
                    }
                }
            }
        }
        // Each node's planes: orthonormal free parts, or rows the prescribed velocities satisfy.
        for (node, rows) in normals {
            let p = [0, 1, 2].map(|k| self.prescribed.get(NODE_DOFS * node + k).copied().flatten());
            let mut done: Vec<SymmetryRow> = Vec::new();
            for n in rows {
                let mut c = n;
                for r in &done {
                    let f: f64 = (0..3).map(|k| r.q[k] * if p[k].is_none() { c[k] } else { 0.0 }).sum();
                    for k in 0..3 {
                        c[k] -= f * r.c[k];
                    }
                }
                let q = [0, 1, 2].map(|k| if p[k].is_none() { c[k] } else { 0.0 });
                let length = q.iter().map(|x| x * x).sum::<f64>().sqrt();
                if length > REPEATED {
                    done.push(SymmetryRow { node, q: q.map(|x| x / length), c: c.map(|x| x / length) });
                    continue;
                }
                // Nothing free is left: the prescribed velocities must satisfy the row.
                let along: f64 = (0..3).map(|k| c[k] * p[k].unwrap_or(0.0)).sum();
                let size: f64 = (0..3).map(|k| (n[k] * p[k].unwrap_or(0.0)).abs()).sum();
                if along.abs() > ROUND_OFF * size {
                    let label = self.labels.get(node).copied().unwrap_or_default();
                    return Err(invalid(format!("node {label} moves off a symmetry plane (normal {n:?}) by its prescribed velocities")));
                }
            }
            self.planes.extend(done);
        }
        Ok(())
    }

    /// The thermal strain rates of the step's temperature loads: each element's mean nodal
    /// temperature change from its material's reference (nodes no load reaches stay at it) is
    /// reached linearly over the step, like a prescribed displacement, so the thermal strain
    /// α ΔT grows at the rate α ΔT / end time (no thermal shock at the start).
    fn temperatures(&mut self, model: &Model, resolved: &Resolved, step: &Step, end_time: f64) -> Result<()> {
        let mut field: Vec<Option<f64>> = vec![None; self.labels.len()];
        let mut any = false;
        for name in &step.loads {
            let Some(Load::Temperature { target, value }) = model.load(name) else { continue };
            any = true;
            for &label in resolved.target_nodes(target) {
                let Some(slot) = self.node(label).and_then(|n| field.get_mut(n)) else { continue };
                match *slot {
                    Some(old) if old != *value => {
                        return Err(invalid(format!("node {label} is given both {old} and {value} (temperature load `{name}` and another)")));
                    }
                    _ => *slot = Some(*value),
                }
            }
        }
        if !any {
            return Ok(());
        }
        if !self.rods.is_empty() {
            self.warnings
                .push("trusses take no thermal strain in explicit dynamics yet; the temperature loads strain the solid elements only".into());
        }
        let materials = &self.materials;
        let rate = |nodes: &[usize], material: usize| -> f64 {
            let Some((alpha, reference)) = materials.get(material).and_then(|m| m.expansion) else { return 0.0 };
            let dt: f64 = nodes.iter().map(|n| field.get(*n).copied().flatten().map_or(0.0, |t| t - reference)).sum::<f64>() / nodes.len() as f64;
            alpha * dt / end_time
        };
        for h in &mut self.hexes {
            h.thermal_rate = rate(&h.nodes, h.material);
        }
        for t in &mut self.tets {
            t.thermal_rate = rate(&t.nodes, t.material);
        }
        if self.hexes.iter().map(|h| h.thermal_rate).chain(self.tets.iter().map(|t| t.thermal_rate)).any(|r| !r.is_finite()) {
            return Err(invalid("the thermal strain rate is not finite"));
        }
        Ok(())
    }

    /// The constant load vector: forces and tractions as the model distributes them, consistent
    /// pressure loads, gravity and body forces with the lumped (row-sum) weights.
    fn load(&mut self, model: &Model, resolved: &Resolved, step: &Step) -> Result<()> {
        let mut f = std::mem::take(&mut self.external);
        for name in &step.loads {
            let load = model.load(name).ok_or_else(|| invalid(format!("load `{name}` is missing")))?;
            let off_mesh = |l: u32| invalid(format!("load `{name}`: node {l} belongs to no element"));
            match load {
                Load::Force { target: Target::ReferencePoint(rp), .. } => {
                    return Err(Error::NotSupported(format!("loads on reference points (`{name}` on `{rp}`)")));
                }
                Load::Moment { .. } => {
                    return Err(Error::NotSupported(format!("moments (`{name}`): solid nodes have no rotations")));
                }
                Load::SurfaceFlux { .. } | Load::Convection { .. } | Load::HeatSource { .. } => {
                    return Err(Error::NotSupported(format!("{} loads (`{name}`) in explicit dynamics", load.kind())));
                }
                // Thermal strains, set by `temperatures`.
                Load::Temperature { .. } => {}
                Load::Force { .. } | Load::Traction { .. } => {
                    for (label, v) in resolved.nodal_loads.get(name).map(Vec::as_slice).unwrap_or_default() {
                        let node = self.node(*label).ok_or_else(|| off_mesh(*label))?;
                        add_at(&mut f, node, v);
                    }
                }
                Load::Pressure { set, magnitude } => {
                    for fc in resolved.sets.get(set).map(|s| s.faces.as_slice()).unwrap_or_default() {
                        let labels = model
                            .face_nodes(*fc)
                            .ok_or_else(|| invalid(format!("load `{name}`: element {} face {} is missing", fc.element, fc.face)))?;
                        let idx: Vec<usize> = labels.iter().map(|l| self.node(*l).ok_or_else(|| off_mesh(*l))).collect::<Result<_>>()?;
                        let xyz: Vec<[f64; 3]> = idx.iter().filter_map(|i| self.xyz.get(*i).copied()).collect();
                        let forces = face::pressure(&xyz, *magnitude).map_err(|e| Error::element(fc.element, e))?;
                        for (i, v) in idx.iter().zip(forces) {
                            add_at(&mut f, *i, &v);
                        }
                    }
                }
                Load::Gravity { acceleration, set } => {
                    let chosen: Option<&[u32]> = set.as_ref().map(|s| resolved.sets.get(s).map(|r| r.elements.as_slice()).unwrap_or_default());
                    self.each_element(chosen, |nodes, weights, rho| {
                        for (node, w) in nodes.iter().zip(weights) {
                            add_at(&mut f, *node, &acceleration.map(|a| rho * w * a));
                        }
                    });
                    // Point masses take gravity on the whole model.
                    if set.is_none() {
                        for (d, (x, m)) in f.iter_mut().zip(&self.discrete.mass).enumerate() {
                            *x += m * acceleration.get(crate::discrete::node_direction(d).1).copied().unwrap_or(0.0);
                        }
                    }
                }
                Load::BodyForce { set, vector } => {
                    let chosen = resolved.sets.get(set).map(|r| r.elements.as_slice()).unwrap_or_default();
                    self.each_element(Some(chosen), |nodes, weights, _| {
                        for (node, w) in nodes.iter().zip(weights) {
                            add_at(&mut f, *node, &vector.map(|b| w * b));
                        }
                    });
                }
                Load::Centrifugal { angular_velocity, point, axis, set } => {
                    let chosen: Option<&[u32]> = set.as_ref().map(|s| resolved.sets.get(s).map(|r| r.elements.as_slice()).unwrap_or_default());
                    let xyz = &self.xyz;
                    self.each_element(chosen, |nodes, weights, rho| {
                        for (node, w) in nodes.iter().zip(weights) {
                            let x = xyz.get(*node).copied().unwrap_or_default();
                            let a = fieldcraft_model::centrifugal_acceleration(*angular_velocity, *point, *axis, x);
                            add_at(&mut f, *node, &a.map(|c| rho * w * c));
                        }
                    });
                }
            }
        }
        self.external = f;
        Ok(())
    }

    /// Visit the elements (all, or those listed): their nodes, nodal volumes and density.
    fn each_element(&self, only: Option<&[u32]>, mut visit: impl FnMut(&[usize], &[f64], f64)) {
        let wanted = |label: u32| only.is_none_or(|l| l.binary_search(&label).is_ok());
        let density = |m: usize| self.materials.get(m).map_or(0.0, |m| m.density);
        for h in self.hexes.iter().filter(|h| wanted(h.label)) {
            visit(&h.nodes, &h.kinematics.nodal_volume, density(h.material));
        }
        for t in self.tets.iter().filter(|t| wanted(t.label)) {
            visit(&t.nodes, &t.kinematics.nodal_volume, density(t.material));
        }
        for r in self.rods.iter().filter(|r| wanted(r.label)) {
            visit(&r.nodes, &r.nodal_volume, density(r.material));
        }
    }

    /// The analysed nodes of a resolved set.
    fn nodes_of(&self, resolved: &Resolved, set: &str) -> Vec<usize> {
        resolved.sets.get(set).map(|s| s.nodes.as_slice()).unwrap_or_default().iter().filter_map(|l| self.node(*l)).collect()
    }

    /// Initial velocities, rigid walls and history nodes of the step.
    fn conditions(&mut self, resolved: &Resolved, settings: &ExplicitDynamic) -> Result<()> {
        for iv in &settings.initial_velocities {
            for node in self.nodes_of(resolved, &iv.set) {
                for (k, c) in iv.velocity.iter().enumerate() {
                    if let Some(v) =
                        DofKind::translation(k).and_then(|kind| self.layout.dof(node, kind)).and_then(|d| self.initial_velocity.get_mut(d))
                    {
                        *v = *c;
                    }
                }
            }
        }
        // A prescribed velocity holds from the start, and so do the symmetry planes.
        for (v, p) in self.initial_velocity.iter_mut().zip(&self.prescribed) {
            if let Some(p) = p {
                *v = *p;
            }
        }
        for row in &self.planes {
            row.project(&mut self.initial_velocity);
        }
        // The model's extent (relative, so the test reads the same in any unit of length).
        let size = self.xyz.iter().flatten().fold(0.0f64, |m, c| m.max(c.abs()));
        for w in &settings.rigid_walls {
            let len = w.normal.iter().map(|c| c * c).sum::<f64>().sqrt();
            let normal = w.normal.map(|c| c / len);
            let nodes = self.nodes_of(resolved, &w.set);
            for &i in &nodes {
                let x = self.xyz.get(i).copied().unwrap_or_default();
                let gap: f64 = (0..3).map(|k| (x[k] - w.point[k]) * normal[k]).sum();
                if gap < -ROUND_OFF * size {
                    let label = self.labels.get(i).copied().unwrap_or(0);
                    return Err(invalid(format!("node {label} of set `{}` starts behind its rigid wall (gap {gap:e})", w.set)));
                }
            }
            self.walls.push(Wall { point: w.point, normal, nodes, friction: w.friction });
        }
        if let Some(set) = &settings.history_nodes {
            let nodes = self.nodes_of(resolved, set);
            if nodes.len() > MAX_HISTORY_NODES {
                return Err(Error::Limit(format!("set `{set}` has {} nodes; histories are recorded for at most {MAX_HISTORY_NODES}", nodes.len())));
            }
            self.history_nodes = nodes;
        }
        Ok(())
    }
}
