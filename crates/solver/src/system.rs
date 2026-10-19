//! From the resolved model to the discrete system: the solid, 2D and truss elements with their
//! materials and formulations, point masses, springs and dashpots ([`crate::discrete`]), the
//! reference points couplings or discrete elements use, the degrees of freedom, prescribed values and the
//! load vector of a step.
//!
//! Degrees of freedom ([`DofLayout`]): the kinds each mesh node's elements declare (three
//! translations for solids, Ux and Uy for 2D elements), nodes by label, then six per coupled reference point (points by
//! name): translations along X, Y, Z and rotations about them.

use fieldcraft_fem::beam::{Beam, BeamSection};
use fieldcraft_fem::constitutive::MaterialModel;
use fieldcraft_fem::dense::Matrix;
use fieldcraft_fem::element::{Formulation, Solid, ThermalStrain};
use fieldcraft_fem::finite::{Kinematics, Temperature};
use fieldcraft_fem::plane::{Plane, Theory};
use fieldcraft_fem::plasticity::J2;
use fieldcraft_fem::shell::{Shell, ShellSection};
use fieldcraft_fem::truss::{Truss, young_modulus};
use fieldcraft_fem::{DofKind, DofLayout, DofSet, Elasticity, ElementDofs, PointStress, elasticity, face, plane};
use fieldcraft_model::nodal::edge_pressure;
use std::collections::BTreeSet;
use std::sync::Arc;

use fieldcraft_model::{Elastic, ElementKind, Hyperelastic, Integration, Load, Model, Resolved, SectionKind, Step, Target, centrifugal_acceleration};

use crate::analysis::MassKind;
use crate::{Error, Result};

/// An element ready for the kernels: a solid, or a 2D element (`theory`).
#[derive(Debug, Clone)]
pub struct Elem {
    pub label: u32,
    pub kind: ElementKind,
    /// Indices into [`Mesh::labels`].
    pub nodes: Vec<usize>,
    pub xyz: Vec<[f64; 3]>,
    pub formulation: Formulation,
    pub d: Elasticity,
    pub density: Option<f64>,
    /// The 2D theory of an element of a 2D section; `None` for solids.
    pub theory: Option<Theory>,
    /// The expansion coefficients along X, Y, Z and the reference temperature, when the
    /// material expands.
    pub expansion: Option<([f64; 3], f64)>,
    /// The large-strain option of the element's material (nonlinear steps).
    pub hyperelastic: Option<Hyperelastic>,
    /// The element's material's plasticity (nonlinear steps), shared by the elements of a
    /// section.
    pub plastic: Option<Arc<J2>>,
}

/// An element's kernel: a solid or a 2D element.
#[derive(Debug, Clone, Copy)]
pub enum Kernel<'a> {
    Solid(Solid<'a>),
    Plane(Plane<'a>),
}

impl Kernel<'_> {
    /// The element's unknowns, in the order of its matrices.
    pub fn dofs(&self) -> ElementDofs {
        match self {
            Kernel::Solid(s) => s.dofs(),
            Kernel::Plane(p) => p.dofs(),
        }
    }

    pub fn stiffness(&self, d: &Elasticity) -> fieldcraft_fem::Result<fieldcraft_fem::dense::Matrix> {
        match self {
            Kernel::Solid(s) => s.stiffness(d),
            Kernel::Plane(p) => p.stiffness(d),
        }
    }

    pub fn mass(&self, density: f64) -> fieldcraft_fem::Result<fieldcraft_fem::dense::Matrix> {
        match self {
            Kernel::Solid(s) => s.mass(density),
            Kernel::Plane(p) => p.mass(density),
        }
    }

    pub fn lumped_mass(&self, density: f64) -> fieldcraft_fem::Result<Vec<f64>> {
        match self {
            Kernel::Solid(s) => s.lumped_mass(density),
            Kernel::Plane(p) => p.lumped_mass(density),
        }
    }

    /// Consistent nodal forces of a uniform force per volume, as 3-vectors (Z zero for 2D).
    pub fn body_force(&self, b: [f64; 3]) -> fieldcraft_fem::Result<Vec<[f64; 3]>> {
        match self {
            Kernel::Solid(s) => s.body_force(b),
            Kernel::Plane(p) => Ok(p.body_force_field(|_| b)?.into_iter().map(|v| [v[0], v[1], 0.0]).collect()),
        }
    }

    /// Consistent nodal forces of a force per volume b(x), as 3-vectors (Z zero for 2D).
    pub fn body_force_field(&self, b: impl Fn([f64; 3]) -> [f64; 3]) -> fieldcraft_fem::Result<Vec<[f64; 3]>> {
        match self {
            Kernel::Solid(s) => s.body_force_field(b),
            Kernel::Plane(p) => Ok(p.body_force_field(b)?.into_iter().map(|v| [v[0], v[1], 0.0]).collect()),
        }
    }

    /// Integration-point stresses for nodal displacements given as 3-vectors, with a thermal
    /// strain when given.
    pub fn stresses_with(&self, d: &Elasticity, u: &[[f64; 3]], thermal: Option<&ThermalStrain>) -> fieldcraft_fem::Result<Vec<PointStress>> {
        match self {
            Kernel::Solid(s) => s.stresses_with(d, u, thermal),
            Kernel::Plane(p) => p.stresses_with(d, &u.iter().map(|v| [v[0], v[1]]).collect::<Vec<_>>(), thermal),
        }
    }

    /// The equivalent nodal forces of a thermal strain, as 3-vectors (Z zero for 2D).
    pub fn thermal_load(&self, d: &Elasticity, thermal: &ThermalStrain) -> fieldcraft_fem::Result<Vec<[f64; 3]>> {
        match self {
            Kernel::Solid(s) => s.thermal_load(d, thermal),
            Kernel::Plane(p) => Ok(p.thermal_load(d, thermal)?.into_iter().map(|v| [v[0], v[1], 0.0]).collect()),
        }
    }

    /// Nonlinear internal forces and tangent ([`fieldcraft_fem::finite`]) for the element's
    /// displacements `u` (its degrees of freedom in order: 3 per node, or 2 for 2D elements),
    /// total Lagrangian or with small strain.
    #[allow(clippy::too_many_arguments)]
    pub fn nonlinear(
        &self,
        m: &MaterialModel,
        u: &[f64],
        kin: Kinematics,
        thermal: Option<&Temperature>,
        committed: &[f64],
        updated: &mut [f64],
        tangent: bool,
    ) -> fieldcraft_fem::Result<fieldcraft_fem::finite::Finite> {
        match self {
            Kernel::Solid(s) => s.nonlinear(m, u.as_chunks::<3>().0, kin, thermal, committed, updated, tangent),
            Kernel::Plane(p) => p.nonlinear(m, u.as_chunks::<2>().0, kin, thermal, committed, updated, tangent),
        }
    }

    /// Cauchy stress and strain (logarithmic at finite strain) at the integration points, as
    /// [`Self::nonlinear`].
    pub fn nonlinear_stresses(
        &self,
        m: &MaterialModel,
        u: &[f64],
        kin: Kinematics,
        thermal: Option<&Temperature>,
        committed: &[f64],
    ) -> fieldcraft_fem::Result<Vec<PointStress>> {
        match self {
            Kernel::Solid(s) => s.nonlinear_stresses(m, u.as_chunks::<3>().0, kin, thermal, committed),
            Kernel::Plane(p) => p.nonlinear_stresses(m, u.as_chunks::<2>().0, kin, thermal, committed),
        }
    }

    /// The initial-stress stiffness of stresses at the stiffness points.
    pub fn geometric_stiffness(&self, stress: &[[f64; 6]]) -> fieldcraft_fem::Result<fieldcraft_fem::dense::Matrix> {
        match self {
            Kernel::Solid(s) => s.geometric_stiffness(stress),
            Kernel::Plane(p) => p.geometric_stiffness(stress),
        }
    }

    pub fn rule(&self) -> fieldcraft_fem::quadrature::Rule {
        match self {
            Kernel::Solid(s) => s.rule(),
            Kernel::Plane(p) => p.rule(),
        }
    }

    pub fn ty(&self) -> &'static fieldcraft_fem::ElementType {
        match self {
            Kernel::Solid(s) => s.ty,
            Kernel::Plane(p) => p.ty,
        }
    }
}

impl Elem {
    pub fn kernel(&self) -> Result<Kernel<'_>> {
        let k = match self.theory {
            None => Solid::new(self.kind, &self.xyz, self.formulation).map(Kernel::Solid),
            Some(t) => Plane::new(self.kind, &self.xyz, self.formulation, t).map(Kernel::Plane),
        };
        k.map_err(|e| Error::element(self.label, e))
    }

    /// The element's thermal strain in a temperature field (nodes the field leaves out are at
    /// the reference temperature); `None` when the material does not expand or nothing changes.
    pub fn thermal_strain(&self, field: &[Option<f64>]) -> Option<ThermalStrain> {
        thermal_strain(self.expansion, &self.nodes, field)
    }
}

/// A truss element (rod2, rod3) ready for the kernel.
#[derive(Debug, Clone)]
pub struct Bar {
    pub label: u32,
    pub kind: ElementKind,
    /// Indices into [`Mesh::labels`] (layout slots).
    pub nodes: Vec<usize>,
    pub xyz: Vec<[f64; 3]>,
    pub young: f64,
    pub area: f64,
    pub density: Option<f64>,
    /// The expansion coefficients along X, Y, Z and the reference temperature, when the
    /// material expands.
    pub expansion: Option<([f64; 3], f64)>,
}

impl Bar {
    pub fn truss(&self) -> Result<Truss<'_>> {
        Truss::new(self.kind, &self.xyz).map_err(|e| Error::element(self.label, e))
    }

    /// The member's thermal strain in a temperature field, as for [`Elem::thermal_strain`].
    pub fn thermal_strain(&self, field: &[Option<f64>]) -> Option<ThermalStrain> {
        thermal_strain(self.expansion, &self.nodes, field)
    }
}

/// The thermal strain of an element of nodes `nodes` in a temperature field (nodes the field
/// leaves out are at the reference temperature); `None` when the material does not expand or
/// nothing changes.
fn thermal_strain(expansion: Option<([f64; 3], f64)>, nodes: &[usize], field: &[Option<f64>]) -> Option<ThermalStrain> {
    let (alpha, reference) = expansion?;
    let delta: Vec<f64> = nodes.iter().map(|n| field.get(*n).copied().flatten().map_or(0.0, |t| t - reference)).collect();
    delta.iter().any(|d| *d != 0.0).then_some(ThermalStrain { alpha, delta })
}

/// A solid or 2D element of a section, by model labels, before the mesh numbers its nodes.
struct SolidEntry {
    label: u32,
    kind: ElementKind,
    nodes: Vec<u32>,
    formulation: Formulation,
    d: Elasticity,
    density: Option<f64>,
    theory: Option<Theory>,
    dofs: ElementDofs,
    expansion: Option<([f64; 3], f64)>,
    hyperelastic: Option<Hyperelastic>,
    plastic: Option<Arc<J2>>,
}

/// What a shell or beam element is made of.
#[derive(Debug, Clone, Copy)]
pub enum Section {
    Shell(ShellSection),
    Beam { section: BeamSection, orientation: [f64; 3] },
}

/// A shell or beam element ready for its kernel: six degrees of freedom per node.
#[derive(Debug, Clone)]
pub struct Structural {
    pub label: u32,
    pub kind: ElementKind,
    /// Indices into [`Mesh::labels`] (layout slots).
    pub nodes: Vec<usize>,
    pub xyz: Vec<[f64; 3]>,
    pub density: Option<f64>,
    pub section: Section,
}

/// The kernel of a shell or beam element.
#[derive(Debug, Clone)]
pub enum StructuralKernel {
    Shell(Shell),
    Beam(Beam),
}

impl Structural {
    /// The unknowns every node of a shell or beam carries.
    pub const KINDS: DofSet = DofSet::STRUCTURAL;

    pub fn kernel(&self) -> Result<StructuralKernel> {
        let k = match self.section {
            Section::Shell(s) => Shell::new(self.kind, &self.xyz, s).map(StructuralKernel::Shell),
            Section::Beam { section, orientation } => Beam::new(&self.xyz, orientation, section).map(StructuralKernel::Beam),
        };
        k.map_err(|e| Error::element(self.label, e))
    }
}

impl StructuralKernel {
    pub fn stiffness(&self) -> fieldcraft_fem::Result<Matrix> {
        match self {
            StructuralKernel::Shell(s) => s.stiffness(),
            StructuralKernel::Beam(b) => Ok(b.stiffness()),
        }
    }

    /// The consistent or lumped mass over every element degree of freedom.
    pub fn mass(&self, density: f64, kind: MassKind) -> fieldcraft_fem::Result<Matrix> {
        Ok(match (self, kind) {
            (StructuralKernel::Shell(s), MassKind::Consistent) => s.mass(density)?,
            (StructuralKernel::Shell(s), MassKind::Lumped) => {
                let diagonal = s.lumped_mass(density)?;
                let mut m = Matrix::zeros(diagonal.len());
                for (i, v) in diagonal.iter().enumerate() {
                    m.add(i, i, *v);
                }
                m
            }
            (StructuralKernel::Beam(b), MassKind::Consistent) => b.mass(density),
            (StructuralKernel::Beam(b), MassKind::Lumped) => b.lumped_mass(density),
        })
    }

    /// Consistent nodal forces and moments of a uniform force per volume, six per node.
    pub fn body_force(&self, b: [f64; 3]) -> fieldcraft_fem::Result<Vec<f64>> {
        match self {
            StructuralKernel::Shell(s) => s.body_force(b),
            StructuralKernel::Beam(beam) => Ok(beam.body_force(b)),
        }
    }
}

/// The translations of a model's nodes: X and Y in a 2D model, else X, Y and Z. Trusses and
/// translational discrete elements declare these.
pub fn translation_kinds(planar: bool) -> DofSet {
    if planar { DofSet::PLANE } else { DofSet::TRANSLATIONS }
}

/// Keep the rows and columns of a matrix over three translations per node (node by node) that
/// `kinds` has: the in-plane block of a truss in a 2D model, or the matrix itself.
pub fn translation_block(data: &[f64], nodes: usize, kinds: DofSet) -> Vec<f64> {
    if kinds == DofSet::TRANSLATIONS {
        return data.to_vec();
    }
    let keep: Vec<usize> = (0..3 * nodes).filter(|i| DofKind::translation(i % 3).is_some_and(|k| kinds.contains(k))).collect();
    let width = 3 * nodes;
    keep.iter().flat_map(|i| keep.iter().map(move |j| data.get(i * width + j).copied().unwrap_or(0.0))).collect()
}

/// A reference point some coupling ties to the mesh.
#[derive(Debug, Clone, PartialEq)]
pub struct Point {
    pub name: String,
    pub label: u32,
    pub xyz: [f64; 3],
}

/// The analysed mesh: the nodes that carry degrees of freedom (those of the solid elements),
/// sorted by label, the elements sorted by label, the coupled reference points by name, and the
/// numbering of their degrees of freedom (mesh node `n` is layout slot `n`, reference point `k`
/// is [`DofLayout::point_slot`]`(k)`).
#[derive(Debug, Clone)]
pub struct Mesh {
    pub labels: Vec<u32>,
    pub xyz: Vec<[f64; 3]>,
    pub elements: Vec<Elem>,
    pub points: Vec<Point>,
    pub layout: DofLayout,
    /// The 2D theory when this is a 2D model (a model is 2D as a whole).
    pub theory: Option<Theory>,
    pub trusses: Vec<Bar>,
    /// Shells and beams, by label.
    pub structural: Vec<Structural>,
    /// What the solver changed or ignored in the model's sections.
    pub warnings: Vec<String>,
    /// Point masses, springs and dashpots ([`crate::discrete`]).
    pub discrete: crate::discrete::Discrete,
}

impl Mesh {
    pub fn node(&self, label: u32) -> Option<usize> {
        self.labels.binary_search(&label).ok()
    }

    pub fn truss(&self, label: u32) -> Option<&Bar> {
        self.trusses.binary_search_by_key(&label, |e| e.label).ok().and_then(|i| self.trusses.get(i))
    }

    /// The index of a shell or beam element in [`Self::structural`].
    pub fn structural_index(&self, label: u32) -> Option<usize> {
        self.structural.binary_search_by_key(&label, |e| e.label).ok()
    }

    /// The kinds trusses and translational discrete elements carry at their nodes.
    pub fn translation_kinds(&self) -> DofSet {
        translation_kinds(self.theory.is_some())
    }

    pub fn element(&self, label: u32) -> Option<&Elem> {
        self.elements.binary_search_by_key(&label, |e| e.label).ok().and_then(|i| self.elements.get(i))
    }

    pub fn point(&self, name: &str) -> Option<usize> {
        self.points.binary_search_by(|p| p.name.as_str().cmp(name)).ok()
    }

    /// The degrees of freedom of the mesh nodes; the reference points' follow them.
    pub fn mesh_dofs(&self) -> usize {
        self.layout.node_dofs()
    }

    /// Every degree of freedom.
    pub fn dofs(&self) -> usize {
        self.layout.dofs()
    }

    /// The degree of freedom of `kind` at layout slot `slot` (a mesh node or reference point).
    pub fn dof(&self, slot: usize, kind: DofKind) -> Result<usize> {
        self.layout.dof(slot, kind).ok_or_else(|| Error::Internal(format!("slot {slot} has no {kind:?} degree of freedom")))
    }

    /// The kind of a degree of freedom.
    pub fn kind(&self, dof: usize) -> Option<DofKind> {
        self.layout.locate(dof).map(|(_, kind)| kind)
    }

    /// What a degree of freedom is, in the user's terms.
    pub fn describe_dof(&self, dof: usize) -> String {
        let Some((slot, kind)) = self.layout.locate(dof) else { return format!("degree of freedom {dof}") };
        match self.labels.get(slot) {
            Some(label) => format!("node {label} in {}", kind.describe()),
            None => {
                let name = self.points.get(slot - self.labels.len()).map_or("?", |p| p.name.as_str());
                format!("reference point `{name}` in {}", kind.describe())
            }
        }
    }

    /// Build it from the model and its resolved view, refusing what the solver cannot represent.
    pub fn new(model: &Model, resolved: &Resolved, max_dofs: usize) -> Result<Mesh> {
        let mut elements = Vec::new();
        let mut trusses = Vec::new();
        let mut structural = Vec::new();
        let mut warnings = Vec::new();
        // A model is 2D as a whole (the model check refuses a mix).
        let planar = model.sections().values().find_map(|s| Theory::from_section(&s.kind));
        if matches!(planar, Some(Theory::Axisymmetric))
            && let Some(what) = model
                .sections()
                .values()
                .any(|s| matches!(s.kind, SectionKind::Rod { .. }))
                .then_some("trusses")
                .or_else(|| (!model.masses().is_empty() || !model.springs().is_empty()).then_some("point masses, springs and dashpots"))
        {
            return Err(Error::NotSupported(format!("{what} in axisymmetric models")));
        }
        for (name, section) in model.sections() {
            let theory = Theory::from_section(&section.kind);
            let rod = match section.kind {
                SectionKind::Rod { area } => Some(area),
                _ => None,
            };
            let shell_or_beam = matches!(section.kind, SectionKind::Shell { .. } | SectionKind::Beam { .. });
            if shell_or_beam && planar.is_some() {
                return Err(Error::NotSupported(format!("{} sections in 2D models (section `{name}`)", section.kind.name())));
            }
            if !matches!(section.kind, SectionKind::Solid {}) && theory.is_none() && rod.is_none() && !shell_or_beam {
                return Err(Error::NotSupported(format!(
                    "{} sections (section `{name}`); the solver has solid, 2D (plane stress, plane strain, axisymmetric) and truss elements so far",
                    section.kind.name()
                )));
            }
            if theory.is_some() != planar.is_some() && rod.is_none() {
                return Err(Error::NotSupported("2D and solid sections in one model".into()));
            }
            let material = model
                .material(&section.material)
                .ok_or_else(|| Error::Invalid(format!("section `{name}`: material `{}` is missing", section.material)))?;
            let elastic =
                material.elastic.as_ref().ok_or_else(|| Error::Invalid(format!("material `{}` has no elastic constants", section.material)))?;
            if section.formulation.hybrid && !matches!(elastic, Elastic::Isotropic { .. }) {
                return Err(Error::NotSupported(format!(
                    "hybrid elements of the orthotropic material `{}` (section `{name}`): the pressure split needs isotropic elasticity",
                    section.material
                )));
            }
            let set = resolved.sets.get(&section.set).ok_or_else(|| Error::Invalid(format!("set `{}` is not resolved", section.set)))?;
            let bad = |e: fieldcraft_fem::Error| match e {
                fieldcraft_fem::Error::NotSupported(m) => Error::NotSupported(format!("section `{name}`: {m}")),
                e => Error::Invalid(format!("section `{name}`, material `{}`: {e}", section.material)),
            };
            let shared = match &section.kind {
                SectionKind::Shell { thickness, offset, .. } => {
                    if section.formulation.integration == Integration::Reduced {
                        warnings.push(format!(
                            "section `{name}`: shells are MITC elements with full integration; reduced integration is solved the same way"
                        ));
                    }
                    Some(Section::Shell(ShellSection::new(*thickness, *offset, elastic).map_err(bad)?))
                }
                SectionKind::Beam { profile, orientation } => {
                    Some(Section::Beam { section: BeamSection::new(profile, elastic).map_err(bad)?, orientation: *orientation })
                }
                _ => None,
            };
            if let Some(body) = shared {
                for &label in &set.elements {
                    let el = model.element(label).ok_or_else(|| Error::Invalid(format!("element {label} is missing")))?;
                    structural.push((label, el.kind, el.nodes.clone(), material.density, body));
                }
                continue;
            }
            if let Some(area) = rod {
                let young =
                    young_modulus(elastic).map_err(|e| Error::Invalid(format!("material `{}` of section `{name}`: {e}", section.material)))?;
                for &label in &set.elements {
                    let el = model.element(label).ok_or_else(|| Error::Invalid(format!("element {label} is missing")))?;
                    let expansion = material.thermal.and_then(|t| t.expansion.map(|a| (a.axes(), t.reference_temperature)));
                    trusses.push((label, el.kind, el.nodes.clone(), young, area, material.density, expansion));
                }
                continue;
            }
            let d = elasticity(elastic).map_err(|e| Error::Invalid(format!("material `{}`: {e}", section.material)))?;
            let plastic = match &material.plastic {
                Some(p) => Some(Arc::new(
                    J2::from_model(elastic, Some(p), material.damage.as_ref())
                        .map_err(|e| Error::Invalid(format!("material `{}`: {e}", section.material)))?,
                )),
                None => None,
            };
            for &label in &set.elements {
                let el = model.element(label).ok_or_else(|| Error::Invalid(format!("element {label} is missing")))?;
                let (formulation, dofs) = match theory {
                    None => Formulation::from_section(el.kind, &section.formulation).map(|f| (f, f.dofs())),
                    Some(_) => plane::formulation(el.kind, &section.formulation).map(|f| (f, ElementDofs::new(DofSet::PLANE))),
                }
                .map_err(|e| Error::element(label, e))?;
                let expansion = material.thermal.and_then(|t| t.expansion.map(|a| (a.axes(), t.reference_temperature)));
                elements.push(SolidEntry {
                    label,
                    kind: el.kind,
                    nodes: el.nodes.clone(),
                    formulation,
                    d,
                    density: material.density,
                    theory,
                    dofs,
                    expansion,
                    hyperelastic: material.hyperelastic,
                    plastic: plastic.clone(),
                });
            }
        }
        if let (Some(_), Some(name)) = (planar, model.couplings().keys().next()) {
            return Err(Error::NotSupported(format!("couplings in 2D models (coupling `{name}`)")));
        }
        elements.sort_by_key(|e| e.label);
        trusses.sort_by_key(|e| e.0);
        structural.sort_by_key(|e| e.0);
        let translations = translation_kinds(planar.is_some());
        // What discrete elements need: at mesh nodes (some only they use), and at reference
        // points, which then get their six degrees of freedom like coupled ones.
        let needs: Vec<(u32, DofSet)> = crate::discrete::needs(resolved, translations).collect();
        let mut labels: Vec<u32> = elements.iter().flat_map(|e| e.nodes.iter().copied()).collect();
        labels.extend(trusses.iter().flat_map(|e| e.2.iter().copied()));
        labels.extend(structural.iter().flat_map(|e| e.2.iter().copied()));
        labels.extend(needs.iter().map(|n| n.0).filter(|l| model.node(*l).is_some()));
        labels.sort_unstable();
        labels.dedup();
        let mut coupled: BTreeSet<&str> = model.couplings().values().map(|c| c.reference_point.as_str()).collect();
        let discrete_labels: BTreeSet<u32> = needs.iter().map(|n| n.0).collect();
        coupled.extend(model.reference_points().iter().filter(|(_, p)| discrete_labels.contains(&p.label)).map(|(n, _)| n.as_str()));
        let points: Vec<Point> = coupled
            .into_iter()
            .map(|name| {
                let p = model.reference_point(name).ok_or_else(|| Error::Invalid(format!("reference point `{name}` is missing")))?;
                Ok(Point { name: name.to_string(), label: p.label, xyz: p.xyz })
            })
            .collect::<Result<_>>()?;
        let mut layout = DofLayout::builder();
        for e in &elements {
            layout.element(&e.nodes, e.dofs);
        }
        for e in &trusses {
            layout.element(&e.2, fieldcraft_fem::ElementDofs::new(translations));
        }
        for e in &structural {
            layout.element(&e.2, ElementDofs::new(Structural::KINDS));
        }
        for (label, kinds) in needs.iter().filter(|n| model.node(n.0).is_some()) {
            layout.node(*label, *kinds);
        }
        for p in &points {
            layout.point(p.label, DofSet::STRUCTURAL);
        }
        let layout = layout.build();
        let dofs = layout.dofs();
        if dofs > max_dofs {
            return Err(Error::Limit(format!("the model has {dofs} degrees of freedom; the limit is {max_dofs}")));
        }
        let xyz: Vec<[f64; 3]> = labels
            .iter()
            .map(|l| model.node(*l).map(|n| n.xyz).ok_or_else(|| Error::Invalid(format!("node {l} is missing"))))
            .collect::<Result<_>>()?;
        let mut mesh = Mesh {
            labels,
            xyz,
            elements: Vec::with_capacity(elements.len()),
            points,
            layout,
            theory: planar,
            trusses: Vec::with_capacity(trusses.len()),
            structural: Vec::with_capacity(structural.len()),
            warnings,
            discrete: crate::discrete::Discrete::default(),
        };
        for e in elements {
            let idx: Vec<usize> = e.nodes.iter().filter_map(|n| mesh.node(*n)).collect();
            let exyz: Vec<[f64; 3]> = idx.iter().filter_map(|i| mesh.xyz.get(*i).copied()).collect();
            mesh.elements.push(Elem {
                label: e.label,
                kind: e.kind,
                nodes: idx,
                xyz: exyz,
                formulation: e.formulation,
                d: e.d,
                density: e.density,
                theory: e.theory,
                expansion: e.expansion,
                hyperelastic: e.hyperelastic,
                plastic: e.plastic,
            });
        }
        for (label, kind, nodes, young, area, density, expansion) in trusses {
            let idx: Vec<usize> = nodes.iter().filter_map(|n| mesh.node(*n)).collect();
            let exyz: Vec<[f64; 3]> = idx.iter().filter_map(|i| mesh.xyz.get(*i).copied()).collect();
            mesh.trusses.push(Bar { label, kind, nodes: idx, xyz: exyz, young, area, density, expansion });
        }
        for (label, kind, nodes, density, section) in structural {
            let idx: Vec<usize> = nodes.iter().filter_map(|n| mesh.node(*n)).collect();
            let exyz: Vec<[f64; 3]> = idx.iter().filter_map(|i| mesh.xyz.get(*i).copied()).collect();
            mesh.structural.push(Structural { label, kind, nodes: idx, xyz: exyz, density, section });
        }
        mesh.discrete = crate::discrete::Discrete::new(&mesh, resolved)?;
        Ok(mesh)
    }

    /// The layout slots (mesh nodes, reference points) every element and discrete element
    /// couples, for the matrix pattern.
    pub fn blocks(&self) -> impl Iterator<Item = std::borrow::Cow<'_, [usize]>> {
        let elements = self
            .elements
            .iter()
            .map(|e| e.nodes.as_slice())
            .chain(self.trusses.iter().map(|t| t.nodes.as_slice()))
            .chain(self.structural.iter().map(|s| s.nodes.as_slice()));
        elements.map(std::borrow::Cow::Borrowed).chain(self.discrete.couplings().map(|dofs| {
            let mut slots: Vec<usize> = dofs.iter().filter_map(|d| self.layout.locate(*d).map(|(slot, _)| slot)).collect();
            slots.dedup();
            std::borrow::Cow::Owned(slots)
        }))
    }

    /// The slots of the matrix pattern: the mesh nodes, and the reference points too when a
    /// discrete element acts on them (couplings tie points through constraints, not the matrix).
    pub fn block_count(&self) -> usize {
        if self.discrete.couplings().flatten().any(|d| d >= self.mesh_dofs()) { self.layout.slots() } else { self.labels.len() }
    }
}

/// Loads and reactions balance when they differ by less than this fraction of
/// [`balance_scale`]; heat supplied and removed, of the gross heat flows.
pub const BALANCE_TOLERANCE: f64 = 1e-6; // unit-factors: relative to the balance's scale

/// The scale of the loads-and-reactions balance check: the largest load or reaction, nodal ones
/// too, since reactions that cancel out (a part squeezed between held faces) still set it.
pub(crate) fn balance_scale(applied: &[f64; 3], total: &[f64; 3], reaction: &[[f64; 3]]) -> f64 {
    let largest = reaction.iter().flatten().fold(0.0f64, |m, v| m.max(v.abs()));
    applied.iter().chain(total).fold(largest, |m, v| m.max(v.abs()))
}

/// Refuse contact in an analysis that cannot take it (linear static, frequency, buckling): a
/// model's contacts apply to all its steps.
pub fn refuse_contact(model: &Model, analysis: &str) -> Result<()> {
    match model.contacts().keys().next() {
        Some(name) => Err(Error::NotSupported(format!(
            "contact (`{name}`) in {analysis} steps: a model's contacts apply to all its steps, so solve its nonlinear static or explicit dynamic steps, or remove the contacts to solve this one"
        ))),
        None => Ok(()),
    }
}

/// The prescribed value of every degree of freedom (`None` = free) under a step's boundary
/// conditions; `values: false` holds them at zero (frequency steps). Rotations of solid nodes do
/// not exist and are ignored, like the directions of a reference point no coupling ties
/// (`active`); a boundary condition on a reference point that no coupling ties to the mesh holds
/// nothing and is reported in `warnings`.
pub fn prescribed(
    model: &Model,
    resolved: &Resolved,
    mesh: &Mesh,
    step: &Step,
    values: bool,
    active: &[bool],
    warnings: &mut Vec<String>,
) -> Result<Vec<Option<f64>>> {
    let mut out: Vec<Option<f64>> = vec![None; mesh.dofs()];
    for name in &step.boundary_conditions {
        let bc = model.boundary_condition(name).ok_or_else(|| Error::Invalid(format!("boundary condition `{name}` is missing")))?;
        if let Target::ReferencePoint(rp) = bc.target()
            && mesh.point(rp).is_none()
        {
            warnings.push(format!(
                "boundary condition `{name}` acts on reference point `{rp}`, which no coupling ties to the mesh and no spring or point mass uses; it holds nothing"
            ));
            continue;
        }
        for slot in resolved.target_nodes(bc.target()).iter().filter_map(|l| mesh.layout.slot(*l)) {
            for kind in mesh.layout.kinds(slot).iter() {
                let Some(v) = bc.value(kind.number()) else { continue };
                let index = mesh.dof(slot, kind)?;
                if !active.get(index).copied().unwrap_or(false) {
                    continue;
                }
                let v = if values { v } else { 0.0 };
                let slot = out.get_mut(index).ok_or_else(|| Error::Internal("dof index".into()))?;
                match *slot {
                    Some(old) if old != v => {
                        return Err(Error::Invalid(format!(
                            "{} is prescribed both {old} and {v} (boundary condition `{name}` and another)",
                            mesh.describe_dof(index)
                        )));
                    }
                    _ => *slot = Some(v),
                }
            }
        }
    }
    // Axisymmetric: a node on the axis cannot move radially.
    if mesh.theory == Some(Theory::Axisymmetric) {
        for (n, p) in mesh.xyz.iter().enumerate() {
            let index = mesh.dof(n, DofKind::Ux)?;
            if p[0] != 0.0 || !active.get(index).copied().unwrap_or(false) {
                continue;
            }
            let slot = out.get_mut(index).ok_or_else(|| Error::Internal("dof index".into()))?;
            match *slot {
                Some(v) if v != 0.0 && values => {
                    return Err(Error::Invalid(format!(
                        "{} is on the axis, where the radial displacement is zero, not {v}",
                        mesh.describe_dof(index)
                    )));
                }
                _ => *slot = Some(0.0),
            }
        }
    }
    Ok(out)
}

/// Add a vector's components to the degrees of freedom `dofs` (the first `dofs.len()`
/// components).
fn add(f: &mut [f64], dofs: &[usize], v: [f64; 3]) {
    for (d, c) in dofs.iter().zip(v) {
        if let Some(x) = f.get_mut(*d) {
            *x += c;
        }
    }
}

/// The translations (or rotations) of layout slot `slot`, if it has them: X, Y, Z, or X and Y
/// for a node of a 2D model.
fn vector_dofs(mesh: &Mesh, slot: usize, rotation: bool) -> Option<Vec<usize>> {
    let kind = |c| if rotation { DofKind::rotation(c) } else { DofKind::translation(c) };
    let dofs: Vec<Option<usize>> = [0, 1, 2].map(|c| kind(c).and_then(|k| mesh.layout.dof(slot, k))).into();
    match dofs.as_slice() {
        [Some(a), Some(b), Some(c)] => Some(vec![*a, *b, *c]),
        [Some(a), Some(b), None] if !rotation && mesh.theory.is_some() => Some(vec![*a, *b]),
        _ => None,
    }
}

/// The translations of mesh node `node`.
pub(crate) fn translations(mesh: &Mesh, node: usize) -> Result<Vec<usize>> {
    vector_dofs(mesh, node, false).ok_or_else(|| Error::Internal(format!("mesh node {node} has no translations")))
}

/// The load vector of a static step: nodal forces and tractions as the model distributes them,
/// forces and moments on coupled reference points, consistent pressure loads, gravity and body
/// forces.
pub fn loads(model: &Model, resolved: &Resolved, mesh: &Mesh, step: &Step, active: &[bool]) -> Result<Vec<f64>> {
    let mut f = vec![0.0; mesh.dofs()];
    for name in &step.loads {
        let load = model.load(name).ok_or_else(|| Error::Invalid(format!("load `{name}` is missing")))?;
        match load {
            Load::Force { .. } | Load::Moment { .. } | Load::Traction { .. } => {
                if let Some(rp) = load.reference_point()
                    && mesh.point(rp).is_none()
                {
                    return Err(Error::Invalid(format!(
                        "load `{name}` acts on reference point `{rp}`, which no coupling ties to the mesh and no spring or point mass uses"
                    )));
                }
                // Moments act on the rotations, which reference points and the nodes of shells
                // and beams have.
                let rotation = matches!(load, Load::Moment { .. });
                for (label, v) in resolved.nodal_loads.get(name).map(Vec::as_slice).unwrap_or_default() {
                    let slot =
                        mesh.layout.slot(*label).ok_or_else(|| Error::Invalid(format!("load `{name}`: node {label} belongs to no solid element")))?;
                    let dofs = vector_dofs(mesh, slot, rotation).ok_or_else(|| {
                        Error::Invalid(format!(
                            "moment `{name}`: node {label} has no rotations (only shell and beam nodes do); apply it to a coupled reference point"
                        ))
                    })?;
                    for (&dof, x) in dofs.iter().zip(v) {
                        if *x != 0.0 && !active.get(dof).copied().unwrap_or(false) {
                            return Err(Error::Invalid(format!(
                                "load `{name}`: {} is not coupled to the mesh (its coupling leaves that direction free)",
                                mesh.describe_dof(dof)
                            )));
                        }
                    }
                    add(&mut f, &dofs, *v);
                }
            }
            Load::Pressure { set, magnitude } if mesh.theory.is_some() => {
                if resolved.sets.get(set).is_some_and(|s| !s.faces.is_empty()) {
                    return Err(Error::Invalid(format!(
                        "load `{name}`: a 2D model's pressures act on element edges, and set `{set}` has element faces"
                    )));
                }
                for g in resolved.planar_edges(set) {
                    let e = mesh.element(g.element).ok_or_else(|| Error::Invalid(format!("load `{name}`: element {} is missing", g.element)))?;
                    let labels = model
                        .edge_nodes(g)
                        .ok_or_else(|| Error::Invalid(format!("load `{name}`: element {} edge {} is missing", g.element, g.edge)))?;
                    let idx: Vec<usize> = labels.iter().filter_map(|l| mesh.node(*l)).collect();
                    let xyz: Vec<[f64; 3]> = idx.iter().filter_map(|i| mesh.xyz.get(*i).copied()).collect();
                    let theory = e.theory.ok_or_else(|| Error::Internal("2D element without a theory".into()))?;
                    let forces = edge_pressure(&xyz, theory.extent(), *magnitude)
                        .map_err(|err| Error::element(g.element, fieldcraft_fem::Error::Invalid(err.to_string())))?;
                    for (i, v) in idx.iter().zip(forces) {
                        add(&mut f, &translations(mesh, *i)?, v);
                    }
                }
            }
            Load::Pressure { set, magnitude } => {
                let faces = resolved.sets.get(set).map(|s| s.faces.as_slice()).unwrap_or_default();
                for fc in faces {
                    let labels = model
                        .face_nodes(*fc)
                        .ok_or_else(|| Error::Invalid(format!("load `{name}`: element {} face {} is missing", fc.element, fc.face)))?;
                    let idx: Vec<usize> = labels
                        .iter()
                        .map(|l| mesh.node(*l).ok_or_else(|| Error::Invalid(format!("load `{name}`: node {l} belongs to no solid element"))))
                        .collect::<Result<_>>()?;
                    let xyz: Vec<[f64; 3]> = idx.iter().filter_map(|i| mesh.xyz.get(*i).copied()).collect();
                    let forces = face::pressure(&xyz, *magnitude).map_err(|e| Error::element(fc.element, e))?;
                    for (i, v) in idx.iter().zip(forces) {
                        add(&mut f, &translations(mesh, *i)?, v);
                    }
                }
            }
            Load::Gravity { acceleration, set } => {
                let chosen: Option<&[u32]> = set.as_ref().map(|s| resolved.sets.get(s).map(|r| r.elements.as_slice()).unwrap_or_default());
                for e in &mesh.elements {
                    if chosen.is_some_and(|c| c.binary_search(&e.label).is_err()) {
                        continue;
                    }
                    let rho = e.density.ok_or_else(|| Error::Invalid(format!("load `{name}`: element {} has no density", e.label)))?;
                    body(&mut f, mesh, e, acceleration.map(|a| rho * a))?;
                }
                for t in &mesh.trusses {
                    if chosen.is_some_and(|c| c.binary_search(&t.label).is_err()) {
                        continue;
                    }
                    let rho = t.density.ok_or_else(|| Error::Invalid(format!("load `{name}`: element {} has no density", t.label)))?;
                    bar_body(&mut f, mesh, t, acceleration.map(|a| rho * a))?;
                }
                for e in &mesh.structural {
                    if chosen.is_some_and(|c| c.binary_search(&e.label).is_err()) {
                        continue;
                    }
                    let rho = e.density.ok_or_else(|| Error::Invalid(format!("load `{name}`: element {} has no density", e.label)))?;
                    structural_body(&mut f, mesh, e, acceleration.map(|a| rho * a))?;
                }
                // Point masses take gravity on the whole model.
                if set.is_none() {
                    mesh.discrete.add_gravity(&mut f, *acceleration, mesh);
                }
            }
            // Heat loads act in heat transfer; temperature fields through their thermal strain
            // (add_thermal_loads).
            Load::SurfaceFlux { .. } | Load::Convection { .. } | Load::HeatSource { .. } | Load::Temperature { .. } => {}
            Load::BodyForce { set, vector } => {
                for label in resolved.sets.get(set).map(|r| r.elements.as_slice()).unwrap_or_default() {
                    if let Some(t) = mesh.truss(*label) {
                        bar_body(&mut f, mesh, t, *vector)?;
                        continue;
                    }
                    if let Some(e) = mesh.structural_index(*label).and_then(|i| mesh.structural.get(i)) {
                        structural_body(&mut f, mesh, e, *vector)?;
                        continue;
                    }
                    let e = mesh.element(*label).ok_or_else(|| Error::Invalid(format!("load `{name}`: element {label} is not a solid element")))?;
                    body(&mut f, mesh, e, *vector)?;
                }
            }
            Load::Centrifugal { angular_velocity, point, axis, set } => {
                let chosen: Option<&[u32]> = set.as_ref().map(|s| resolved.sets.get(s).map(|r| r.elements.as_slice()).unwrap_or_default());
                if let Some(t) = mesh.trusses.iter().find(|t| chosen.is_none_or(|c| c.binary_search(&t.label).is_ok())) {
                    return Err(Error::NotSupported(format!("centrifugal loads on trusses (load `{name}`, element {})", t.label)));
                }
                if let Some(e) = mesh.structural.iter().find(|e| chosen.is_none_or(|c| c.binary_search(&e.label).is_ok())) {
                    return Err(Error::NotSupported(format!("centrifugal loads on shells and beams (load `{name}`, element {})", e.label)));
                }
                if set.is_none() && !mesh.discrete.masses.is_empty() {
                    return Err(Error::NotSupported(format!("centrifugal loads on point masses (load `{name}`)")));
                }
                for e in &mesh.elements {
                    if chosen.is_some_and(|c| c.binary_search(&e.label).is_err()) {
                        continue;
                    }
                    let rho = e.density.ok_or_else(|| Error::Invalid(format!("load `{name}`: element {} has no density", e.label)))?;
                    let forces = e
                        .kernel()?
                        .body_force_field(|x| centrifugal_acceleration(*angular_velocity, *point, *axis, x).map(|a| rho * a))
                        .map_err(|err| Error::element(e.label, err))?;
                    for (node, v) in e.nodes.iter().zip(forces) {
                        add(&mut f, &translations(mesh, *node)?, v);
                    }
                }
            }
        }
    }
    Ok(f)
}

/// A force `b` per volume on a truss.
fn bar_body(f: &mut [f64], mesh: &Mesh, t: &Bar, b: [f64; 3]) -> Result<()> {
    let forces = t.truss()?.body_force(b, t.area).map_err(|err| Error::element(t.label, err))?;
    for (node, v) in t.nodes.iter().zip(forces) {
        add(f, &translations(mesh, *node)?, v);
    }
    Ok(())
}

/// Consistent forces and moments of a uniform force per volume on a shell or beam.
fn structural_body(f: &mut [f64], mesh: &Mesh, e: &Structural, b: [f64; 3]) -> Result<()> {
    let forces = e.kernel()?.body_force(b).map_err(|err| Error::element(e.label, err))?;
    let dofs =
        mesh.layout.element_dofs(&e.nodes, Structural::KINDS).ok_or_else(|| Error::Internal(format!("element {} degrees of freedom", e.label)))?;
    for (d, v) in dofs.iter().zip(&forces) {
        if let Some(x) = f.get_mut(*d) {
            *x += v;
        }
    }
    Ok(())
}

/// The temperature field of a static step's temperature loads, per mesh node; `None` when the
/// step has none.
pub fn temperature_loads(model: &Model, resolved: &Resolved, mesh: &Mesh, step: &Step) -> Result<Option<Vec<Option<f64>>>> {
    let mut field: Option<Vec<Option<f64>>> = None;
    for name in &step.loads {
        let Some(Load::Temperature { target, value }) = model.load(name) else { continue };
        let field = field.get_or_insert_with(|| vec![None; mesh.labels.len()]);
        for &label in resolved.target_nodes(target) {
            let Some(n) = mesh.node(label) else { continue };
            hold(field.get_mut(n), *value, |old| format!("node {label} is given both {old} and {value} (temperature load `{name}` and another)"))?;
        }
    }
    Ok(field)
}

/// Set a node's value, refusing a different one already there (`conflict` names both).
pub(crate) fn hold(slot: Option<&mut Option<f64>>, value: f64, conflict: impl FnOnce(f64) -> String) -> Result<()> {
    let slot = slot.ok_or_else(|| Error::Internal("node index".into()))?;
    match *slot {
        Some(old) if old != value => Err(Error::Invalid(conflict(old))),
        _ => {
            *slot = Some(value);
            Ok(())
        }
    }
}

/// Add the equivalent nodal forces of the thermal strains of a temperature field.
pub fn add_thermal_loads(f: &mut [f64], mesh: &Mesh, field: &[Option<f64>]) -> Result<()> {
    add_element_thermal_loads(f, mesh, field, |_| true)?;
    for t in &mesh.trusses {
        let Some(th) = t.thermal_strain(field) else { continue };
        let forces = t.truss()?.thermal_load(t.young, t.area, &th).map_err(|err| Error::element(t.label, err))?;
        for (node, v) in t.nodes.iter().zip(forces) {
            add(f, &translations(mesh, *node)?, v);
        }
    }
    Ok(())
}

/// Add the equivalent nodal forces of the thermal strains of the solid and 2D elements whose
/// index `keep` accepts.
pub fn add_element_thermal_loads(f: &mut [f64], mesh: &Mesh, field: &[Option<f64>], keep: impl Fn(usize) -> bool) -> Result<()> {
    for (_, e) in mesh.elements.iter().enumerate().filter(|(i, _)| keep(*i)) {
        let Some(th) = e.thermal_strain(field) else { continue };
        let forces = e.kernel()?.thermal_load(&e.d, &th).map_err(|err| Error::element(e.label, err))?;
        for (node, v) in e.nodes.iter().zip(forces) {
            add(f, &translations(mesh, *node)?, v);
        }
    }
    Ok(())
}

/// The uniform force per volume a step's gravity and body forces put on each shell and beam
/// element ([`Mesh::structural`] order): beam end forces are recovered with it.
pub fn structural_body_forces(model: &Model, resolved: &Resolved, mesh: &Mesh, step: &Step) -> Result<Vec<[f64; 3]>> {
    let mut out = vec![[0.0; 3]; mesh.structural.len()];
    for name in &step.loads {
        let load = model.load(name).ok_or_else(|| Error::Invalid(format!("load `{name}` is missing")))?;
        match load {
            Load::Gravity { acceleration, set } => {
                let chosen: Option<&[u32]> = set.as_ref().map(|s| resolved.sets.get(s).map(|r| r.elements.as_slice()).unwrap_or_default());
                for (e, b) in mesh.structural.iter().zip(out.iter_mut()) {
                    if chosen.is_some_and(|c| c.binary_search(&e.label).is_err()) {
                        continue;
                    }
                    let rho = e.density.ok_or_else(|| Error::Invalid(format!("load `{name}`: element {} has no density", e.label)))?;
                    *b = [0, 1, 2].map(|k| b[k] + rho * acceleration[k]);
                }
            }
            Load::BodyForce { set, vector } => {
                for label in resolved.sets.get(set).map(|r| r.elements.as_slice()).unwrap_or_default() {
                    if let Some(b) = mesh.structural_index(*label).and_then(|i| out.get_mut(i)) {
                        *b = [0, 1, 2].map(|k| b[k] + vector[k]);
                    }
                }
            }
            _ => {}
        }
    }
    Ok(out)
}

fn body(f: &mut [f64], mesh: &Mesh, e: &Elem, b: [f64; 3]) -> Result<()> {
    let forces = e.kernel()?.body_force(b).map_err(|err| Error::element(e.label, err))?;
    for (node, v) in e.nodes.iter().zip(forces) {
        add(f, &translations(mesh, *node)?, v);
    }
    Ok(())
}

/// The linear solvers use the elastic constants alone: a warning for each material in use whose
/// plasticity, damage or erosion they leave out.
pub fn inelastic_warnings(model: &Model) -> Vec<String> {
    let used: BTreeSet<&str> = model.sections().values().map(|s| s.material.as_str()).collect();
    used.into_iter()
        .filter(|name| model.material(name).is_some_and(|m| m.plastic.is_some() || m.damage.is_some() || m.erosion.is_some()))
        .map(|name| {
            format!(
                "material `{name}`: the linear solver uses its elastic constants only (its plasticity acts in nonlinear static steps and explicit dynamics, its damage and erosion in explicit dynamics)"
            )
        })
        .collect()
}
