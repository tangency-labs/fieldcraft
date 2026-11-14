//! Ordered structural contributions shared by sparse assembly and streamed operators.
//!
//! Support discovery does not construct numerical kernels. Operators retain only the mesh;
//! application holds at most one element matrix, never a global matrix or element cache.

use crate::analysis::{GRAIN, MassKind, ROUND};
use crate::discrete::{MassAt, SpringAt};
use crate::operator::LinearOperator;
use crate::sparse::{Csr, add_block, add_diagonal, add_entry, add_scalar};
use crate::system::{Bar, Elem, Mesh, Structural, translation_block};
use crate::{Error, Result};
use fieldcraft_fem::DofSet;

/// The requested numerical field; evaluating one never evaluates the other.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Field {
    Stiffness,
    Mass(MassKind),
}

/// Identity within an immutable mesh snapshot, in assembly order (families, then stored order).
/// Discrete indices are not persistent identities across mesh changes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Identity {
    Element(u32),
    Truss(u32),
    Structural(u32),
    Spring(usize),
    PointMass(usize),
}

#[derive(Clone, Copy)]
enum Source<'a> {
    Element(&'a Elem),
    Truss(&'a Bar),
    Structural(&'a Structural),
    Spring(&'a SpringAt),
    Mass(&'a MassAt),
}

/// A borrowed contribution. Its support can be inspected without evaluating a kernel.
#[derive(Clone, Copy)]
pub struct Contribution<'a> {
    mesh: &'a Mesh,
    source: Source<'a>,
    identity: Identity,
}

/// Symmetric local storage in [`Contribution::support`] order. `Block` is full row-major;
/// `Scalar` is a per-node square matrix applied independently to each component and `Diagonal`
/// is one value per node repeated for each component. [`Contribution::nodes`] supplies the
/// node/component layout for these compact masses; discrete contributions use `Block`.
pub enum LocalMatrix {
    Block(Vec<f64>),
    Scalar(Vec<f64>),
    Diagonal(Vec<f64>),
}

/// Contributions in exactly the order used by linear sparse assembly. Dashpots are excluded.
pub fn contributions(mesh: &Mesh) -> impl Iterator<Item = Contribution<'_>> {
    mesh.elements
        .iter()
        .map(|e| (Source::Element(e), Identity::Element(e.label)))
        .chain(mesh.trusses.iter().map(|e| (Source::Truss(e), Identity::Truss(e.label))))
        .chain(mesh.structural.iter().map(|e| (Source::Structural(e), Identity::Structural(e.label))))
        .chain(mesh.discrete.springs.iter().enumerate().filter(|(_, s)| !s.dashpot).map(|(i, s)| (Source::Spring(s), Identity::Spring(i))))
        .chain(mesh.discrete.masses.iter().enumerate().map(|(i, m)| (Source::Mass(m), Identity::PointMass(i))))
        .map(move |(source, identity)| Contribution { mesh, source, identity })
}

impl Contribution<'_> {
    pub fn identity(&self) -> Identity {
        self.identity
    }

    /// Node slots and per-node components of compact local matrices, in local DOF order.
    /// Discrete contributions return `None` and always use full `Block` matrices.
    pub fn nodes(&self) -> Option<(&[usize], DofSet)> {
        match self.source {
            Source::Element(e) => Some((&e.nodes, if e.theory.is_some() { DofSet::PLANE } else { DofSet::TRANSLATIONS })),
            Source::Truss(e) => Some((&e.nodes, self.mesh.translation_kinds())),
            Source::Structural(e) => Some((&e.nodes, Structural::KINDS)),
            _ => None,
        }
    }

    /// Global degrees of freedom in local matrix order, including repeated spring terms.
    pub fn support(&self) -> Result<Vec<usize>> {
        let mut dofs = Vec::new();
        self.support_into(&mut dofs)?;
        Ok(dofs)
    }

    /// Fill caller-owned support storage, retaining its capacity across contributions.
    pub fn support_into(&self, dofs: &mut Vec<usize>) -> Result<()> {
        dofs.clear();
        if let Some((nodes, kinds)) = self.nodes() {
            for n in nodes {
                for k in kinds.iter() {
                    dofs.push(self.mesh.layout.dof(*n, k).ok_or_else(|| Error::Internal("contribution degree of freedom".into()))?);
                }
            }
        } else {
            match self.source {
                Source::Spring(s) => dofs.extend(s.terms.iter().map(|t| t.0)),
                Source::Mass(m) => dofs.extend(m.translations.iter().map(|t| t.0).chain(m.rotations.iter().copied())),
                _ => {}
            }
        }
        Ok(())
    }

    fn active(&self, field: Field) -> bool {
        !matches!((self.source, field), (Source::Spring(_), Field::Mass(_)) | (Source::Mass(_), Field::Stiffness))
    }

    fn local_size(&self) -> Result<usize> {
        match self.source {
            Source::Spring(s) => Some(s.terms.len()),
            Source::Mass(m) => m.translations.len().checked_add(m.rotations.len()),
            _ => self.nodes().and_then(|(nodes, kinds)| nodes.len().checked_mul(kinds.len())),
        }
        .ok_or_else(|| Error::Internal("contribution size overflow".into()))
    }

    fn square(n: usize) -> Result<usize> {
        n.checked_mul(n)
            .filter(|size| *size <= isize::MAX as usize / std::mem::size_of::<f64>())
            .ok_or_else(|| Error::Internal("contribution matrix size overflow".into()))
    }

    // One entry convention for discrete local matrices and allocation-free CSR scatter.
    fn discrete_entries(&self, field: Field, mut add: impl FnMut(usize, usize, usize, usize, f64) -> Result<()>) -> Result<()> {
        if !self.active(field) {
            return Ok(());
        }
        match self.source {
            Source::Spring(s) => {
                for (i, (a, ci)) in s.terms.iter().enumerate() {
                    for (j, (b, cj)) in s.terms.iter().enumerate() {
                        add(i, j, *a, *b, s.value * ci * cj)?;
                    }
                }
            }
            Source::Mass(m) => {
                for (i, (d, v)) in m.translations.iter().enumerate() {
                    add(i, i, *d, *d, *v)?;
                }
                let offset = m.translations.len();
                for (i, a) in m.rotations.iter().enumerate() {
                    for (j, b) in m.rotations.iter().enumerate() {
                        let v = m.inertia.get(i).and_then(|row| row.get(j)).copied().ok_or_else(|| Error::Internal("point inertia size".into()))?;
                        if v != 0.0 {
                            add(offset + i, offset + j, *a, *b, v)?;
                        }
                    }
                }
            }
            _ => return Err(Error::Internal("not a discrete contribution".into())),
        }
        Ok(())
    }

    fn scatter_field(&self, field: Field, global: &mut Csr) -> Result<()> {
        if matches!(self.source, Source::Spring(_) | Source::Mass(_)) {
            self.discrete_entries(field, |_, _, a, b, v| add_entry(global, a, b, v))
        } else {
            if let Some(matrix) = self.local_matrix(field)? {
                self.scatter(&matrix, global)?;
            }
            Ok(())
        }
    }

    /// Evaluate only the requested matrix. `None` denotes an identically zero field.
    pub fn local_matrix(&self, field: Field) -> Result<Option<LocalMatrix>> {
        if !self.active(field) {
            return Ok(None);
        }
        let n = self.local_size()?;
        Self::square(n)?;
        let density = |rho: Option<f64>, label| rho.ok_or_else(|| Error::Invalid(format!("element {label} has no density")));
        Ok(Some(match self.source {
            Source::Element(e) => {
                let kernel = e.kernel()?;
                let at = |err| Error::element(e.label, err);
                match field {
                    Field::Stiffness => LocalMatrix::Block(kernel.stiffness(&e.d).map_err(at)?.data),
                    Field::Mass(MassKind::Consistent) => LocalMatrix::Scalar(kernel.mass(density(e.density, e.label)?).map_err(at)?.data),
                    Field::Mass(MassKind::Lumped) => LocalMatrix::Diagonal(kernel.lumped_mass(density(e.density, e.label)?).map_err(at)?),
                }
            }
            Source::Truss(e) => {
                let kernel = e.truss()?;
                let at = |err| Error::element(e.label, err);
                match field {
                    Field::Stiffness => {
                        let data = kernel.stiffness(e.young, e.area).map_err(at)?.data;
                        LocalMatrix::Block(if self.mesh.translation_kinds() == DofSet::TRANSLATIONS {
                            data
                        } else {
                            translation_block(&data, e.nodes.len(), self.mesh.translation_kinds())
                        })
                    }
                    Field::Mass(MassKind::Consistent) => LocalMatrix::Scalar(kernel.mass(density(e.density, e.label)?, e.area).map_err(at)?.data),
                    Field::Mass(MassKind::Lumped) => LocalMatrix::Diagonal(kernel.lumped_mass(density(e.density, e.label)?, e.area).map_err(at)?),
                }
            }
            Source::Structural(e) => {
                let kernel = e.kernel()?;
                let at = |err| Error::element(e.label, err);
                LocalMatrix::Block(match field {
                    Field::Stiffness => kernel.stiffness().map_err(at)?.data,
                    Field::Mass(kind) => kernel.mass(density(e.density, e.label)?, kind).map_err(at)?.data,
                })
            }
            Source::Spring(_) | Source::Mass(_) => {
                let mut data = vec![0.0; Self::square(n)?];
                self.discrete_entries(field, |i, j, _, _, v| {
                    *data.get_mut(i * n + j).ok_or_else(|| Error::Internal("local matrix index".into()))? = v;
                    Ok(())
                })?;
                LocalMatrix::Block(data)
            }
        }))
    }

    fn scatter(&self, matrix: &LocalMatrix, global: &mut Csr) -> Result<()> {
        if let Some((nodes, kinds)) = self.nodes() {
            return match matrix {
                LocalMatrix::Block(data) => add_block(global, &self.mesh.layout, nodes, kinds, data),
                LocalMatrix::Scalar(data) => add_scalar(global, &self.mesh.layout, nodes, kinds, data),
                LocalMatrix::Diagonal(data) => add_diagonal(global, &self.mesh.layout, nodes, kinds, data),
            };
        }
        Err(Error::Internal("discrete contribution requires direct scatter".into()))
    }

    /// Independent local action, overwriting the local output. Springs and point masses act
    /// directly; existing continuum/beam/shell kernels use their compact local matrices.
    pub fn apply_local(&self, field: Field, x: &[f64], y: &mut [f64]) -> Result<()> {
        let n = self.local_size()?;
        if x.len() != n || y.len() != n {
            return Err(Error::Invalid("local operator vector length".into()));
        }
        y.fill(0.0);
        match self.source {
            Source::Spring(s) => {
                if field == Field::Stiffness {
                    let stretch: f64 = s.terms.iter().zip(x).map(|((_, c), v)| c * v).sum();
                    for (out, (_, c)) in y.iter_mut().zip(&s.terms) {
                        *out = s.value * c * stretch;
                    }
                }
                return Ok(());
            }
            Source::Mass(m) => {
                if field != Field::Stiffness {
                    for ((out, input), (_, mass)) in y.iter_mut().zip(x).zip(&m.translations) {
                        *out = mass * input;
                    }
                    let offset = m.translations.len();
                    for (a, _) in m.rotations.iter().enumerate() {
                        for (b, _) in m.rotations.iter().enumerate() {
                            let v =
                                m.inertia.get(a).and_then(|row| row.get(b)).copied().ok_or_else(|| Error::Internal("point inertia size".into()))?;
                            let input = x.get(offset + b).copied().ok_or_else(|| Error::Internal("point mass local index".into()))?;
                            *y.get_mut(offset + a).ok_or_else(|| Error::Internal("point mass local index".into()))? += v * input;
                        }
                    }
                }
                return Ok(());
            }
            _ => {}
        }
        if let Some(matrix) = self.local_matrix(field)? {
            let per = self.nodes().map_or(1, |(_, k)| k.len());
            match matrix {
                LocalMatrix::Block(data) => {
                    if n == 0 || data.len() != Self::square(n)? {
                        return Err(Error::Internal("local matrix size".into()));
                    }
                    for (out, row) in y.iter_mut().zip(data.chunks_exact(n)) {
                        *out = row.iter().zip(x).map(|(a, b)| a * b).sum();
                    }
                }
                LocalMatrix::Scalar(data) => {
                    let nodes = n / per;
                    if data.len() != Self::square(nodes)? {
                        return Err(Error::Internal("local scalar matrix size".into()));
                    }
                    for i in 0..nodes {
                        for j in 0..nodes {
                            for c in 0..per {
                                y[i * per + c] += data[i * nodes + j] * x[j * per + c];
                            }
                        }
                    }
                }
                LocalMatrix::Diagonal(data) => {
                    if data.len() != n / per {
                        return Err(Error::Internal("local diagonal size".into()));
                    }
                    for (i, out) in y.iter_mut().enumerate() {
                        *out = data[i / per] * x[i];
                    }
                }
            }
        }
        Ok(())
    }
}

/// Caller-owned, element-bounded scratch reusable across operator applications.
#[derive(Default)]
pub struct OperatorScratch {
    dofs: Vec<usize>,
    local: Vec<f64>,
    action: Vec<f64>,
}

/// Full-space streamed stiffness or mass operator. Current continuum/beam/shell kernels
/// allocate their compact element matrices; no global matrix or element cache is retained.
pub struct StructuralOperator<'a> {
    mesh: &'a Mesh,
    field: Field,
}
impl<'a> StructuralOperator<'a> {
    pub fn new(mesh: &'a Mesh, field: Field) -> Self {
        Self { mesh, field }
    }

    /// Apply with reusable gather/scatter storage. Output is overwritten.
    pub fn apply_with_scratch(&self, x: &[f64], y: &mut [f64], scratch: &mut OperatorScratch) -> Result<()> {
        if x.len() != self.dimension() || y.len() != self.dimension() {
            return Err(Error::Invalid("structural operator vector length".into()));
        }
        y.fill(0.0);
        let OperatorScratch { dofs, local, action } = scratch;
        for contribution in contributions(self.mesh) {
            if !contribution.active(self.field) {
                continue;
            }
            contribution.support_into(dofs)?;
            local.clear();
            for d in dofs.iter() {
                local.push(x.get(*d).copied().ok_or_else(|| Error::Internal("contribution degree of freedom".into()))?);
            }
            action.resize(dofs.len(), 0.0);
            contribution.apply_local(self.field, local, action)?;
            for (d, v) in dofs.iter().zip(action.iter()) {
                *y.get_mut(*d).ok_or_else(|| Error::Internal("contribution degree of freedom".into()))? += v;
            }
        }
        Ok(())
    }
}
impl LinearOperator for StructuralOperator<'_> {
    fn dimension(&self) -> usize {
        self.mesh.dofs()
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) -> Result<()> {
        self.apply_with_scratch(x, y, &mut OperatorScratch::default())
    }
}

/// Bounded, deterministic parallel evaluation followed by ordered block/scalar scatter.
pub fn assemble_into(mesh: &Mesh, k: &mut Csr, mut mass: Option<(&mut Csr, MassKind)>) -> Result<()> {
    // Keep trusses and discrete elements sequential as in the original assembly path.
    let mut iter = contributions(mesh).peekable();
    let mut elements = Vec::new();
    while let Some(first) = iter.peek().copied() {
        let parallel = matches!(first.source, Source::Element(_) | Source::Structural(_));
        if !parallel {
            if let Some(c) = iter.next() {
                c.scatter_field(Field::Stiffness, k)?;
                if let Some((global, kind)) = mass.as_mut() {
                    c.scatter_field(Field::Mass(*kind), global)?;
                }
            }
            continue;
        }
        let family = std::mem::discriminant(&first.source);
        elements.clear();
        let limit = ROUND;
        while elements.len() < limit && iter.peek().is_some_and(|c| std::mem::discriminant(&c.source) == family) {
            if let Some(c) = iter.next() {
                elements.push(c);
            }
        }
        let matrices = fieldcraft_parallel::map_grained(elements.len(), GRAIN, |i| -> Result<_> {
            let c = elements.get(i).ok_or_else(|| Error::Internal("contribution index".into()))?;
            let stiffness = c.local_matrix(Field::Stiffness)?;
            let mm = mass.as_ref().map(|(_, kind)| c.local_matrix(Field::Mass(*kind))).transpose()?.flatten();
            Ok((stiffness, mm))
        });
        for (c, computed) in elements.iter().zip(matrices) {
            let (stiffness, mm) = computed?;
            if let Some(data) = stiffness {
                c.scatter(&data, k)?;
            }
            if let (Some((global, _)), Some(data)) = (mass.as_mut(), mm) {
                c.scatter(&data, global)?;
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::discrete::Discrete;
    use crate::sparse::NodeGraph;
    use crate::system::Section;
    use fieldcraft_fem::beam::BeamSection;
    use fieldcraft_fem::element::Formulation;
    use fieldcraft_fem::plane::Theory;
    use fieldcraft_fem::shell::ShellSection;
    use fieldcraft_fem::{DofKind, DofLayout, elasticity};
    use fieldcraft_model::{BeamProfile, Elastic, ElementKind};

    fn fixture(planar: bool) -> Mesh {
        let xyz = vec![[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]];
        let theory = planar.then_some(Theory::PlaneStress { thickness: 1.0 });
        let nodes = if planar { vec![0, 1, 2] } else { vec![0, 1, 2, 3] };
        let mut builder = DofLayout::builder();
        for label in 1..=4 {
            builder.node(label, if label == 4 { DofSet::TRANSLATIONS } else { DofSet::STRUCTURAL });
        }
        let layout = builder.build();
        let ux = |slot| layout.dof(slot, DofKind::Ux).unwrap();
        let uy = |slot| layout.dof(slot, DofKind::Uy).unwrap();
        let rotations = [DofKind::Rx, DofKind::Ry, DofKind::Rz].iter().map(|k| layout.dof(0, *k).unwrap()).collect();
        let discrete = Discrete {
            masses: vec![MassAt {
                name: "inertia".into(),
                translations: vec![(ux(0), 2.0), (uy(0), 3.0)],
                rotations,
                inertia: [[4., 0.3, -0.2], [0.3, 5., 0.4], [-0.2, 0.4, 6.]],
            }],
            springs: vec![SpringAt {
                name: "skew".into(),
                dashpot: false,
                grounded: false,
                terms: vec![(ux(0), -0.6), (uy(0), -0.8), (ux(1), 0.6), (uy(1), 0.8)],
                value: 7.0,
            }],
        };
        let elastic = Elastic::Isotropic { young: 1000.0, poisson: 0.25 };
        let structural = if planar {
            vec![]
        } else {
            vec![
                Structural {
                    label: 3,
                    kind: ElementKind::Tri3,
                    nodes: vec![0, 1, 2],
                    xyz: xyz[..3].to_vec(),
                    section: Section::Shell(ShellSection::new(0.1, 0.0, &elastic).unwrap()),
                    density: Some(2.0),
                },
                Structural {
                    label: 4,
                    kind: ElementKind::Beam2,
                    nodes: vec![0, 1],
                    xyz: xyz[..2].to_vec(),
                    section: Section::Beam {
                        section: BeamSection::new(&BeamProfile::Rectangle { width: 0.1, height: 0.2 }, &elastic).unwrap(),
                        orientation: [0., 0., 1.],
                    },
                    density: Some(2.0),
                },
            ]
        };
        Mesh {
            labels: vec![1, 2, 3, 4],
            xyz: xyz.clone(),
            layout,
            theory,
            elements: vec![Elem {
                label: 1,
                kind: if planar { ElementKind::Tri3 } else { ElementKind::Tet4 },
                xyz: nodes.iter().map(|i| xyz[*i]).collect(),
                nodes,
                formulation: Formulation::Full,
                d: elasticity(&elastic).unwrap(),
                density: Some(2.0),
                theory,
                expansion: None,
                hyperelastic: None,
                plastic: None,
            }],
            trusses: vec![Bar {
                label: 2,
                kind: ElementKind::Rod2,
                nodes: vec![1, 2],
                xyz: vec![xyz[1], xyz[2]],
                young: 1000.,
                area: 0.2,
                density: Some(2.0),
                expansion: None,
            }],
            structural,
            discrete,
            points: vec![],
            warnings: vec![],
        }
    }

    fn close(a: &[f64], b: &[f64]) {
        assert_eq!(a.len(), b.len());
        for (a, b) in a.iter().zip(b) {
            assert!((a - b).abs() <= 1e-10 * (1.0 + a.abs() + b.abs()), "{a} != {b}");
        }
    }

    #[test]
    fn streamed_actions_match_sparse_for_every_family_and_mass_kind() {
        for planar in [false, true] {
            let mesh = fixture(planar);
            let graph = NodeGraph::new(mesh.block_count(), mesh.blocks());
            let x: Vec<_> = (0..mesh.dofs()).map(|i| (i as f64 + 0.3).sin()).collect();
            for kind in [MassKind::Consistent, MassKind::Lumped] {
                let mut k = graph.zero_matrix(&mesh.layout);
                let mut m = graph.zero_matrix(&mesh.layout);
                assemble_into(&mesh, &mut k, Some((&mut m, kind))).unwrap();
                for (field, matrix) in [(Field::Stiffness, &k), (Field::Mass(kind), &m)] {
                    let mut expected = vec![0.; mesh.dofs()];
                    matrix.mul(&x, &mut expected);
                    let mut actual = vec![99.; mesh.dofs()];
                    StructuralOperator::new(&mesh, field).apply(&x, &mut actual).unwrap();
                    close(&actual, &expected);
                }
            }
            // Uniform translation is strain-free, including the ungrounded skew spring.
            let mut rigid = vec![0.; mesh.dofs()];
            for slot in 0..mesh.layout.nodes() {
                if let Some(d) = mesh.layout.dof(slot, DofKind::Ux) {
                    rigid[d] = 1.;
                }
            }
            let mut action = vec![0.; mesh.dofs()];
            StructuralOperator::new(&mesh, Field::Stiffness).apply(&rigid, &mut action).unwrap();
            close(&action, &vec![0.; mesh.dofs()]);
            // Total mass along X: continuum + truss + shell + beam + point mass.
            let expected_mass = if planar { 1.0 } else { 1.0 / 3.0 + 0.1 + 0.04 } + 0.4 * 2.0_f64.sqrt() + 2.0;
            for kind in [MassKind::Consistent, MassKind::Lumped] {
                StructuralOperator::new(&mesh, Field::Mass(kind)).apply(&rigid, &mut action).unwrap();
                let total: f64 = rigid.iter().zip(&action).map(|(a, b)| a * b).sum();
                close(&[total], &[expected_mass]);
            }
        }
    }

    #[test]
    fn support_survives_invalid_geometry_and_stiffness_needs_no_density() {
        let mut mesh = fixture(true);
        mesh.elements[0].density = None;
        mesh.trusses[0].density = None;
        let support = contributions(&mesh).next().unwrap().support().unwrap();
        let mut invalid = mesh.clone();
        invalid.elements[0].xyz.fill([0.0; 3]);
        let element = contributions(&invalid).next().unwrap();
        assert_eq!(element.support().unwrap(), support);
        assert!(element.local_matrix(Field::Stiffness).is_err());
        let mut out = vec![0.; mesh.dofs()];
        StructuralOperator::new(&mesh, Field::Stiffness).apply(&vec![0.; mesh.dofs()], &mut out).unwrap();
        assert!(StructuralOperator::new(&mesh, Field::Mass(MassKind::Consistent)).apply(&vec![0.; mesh.dofs()], &mut out).is_err());
        assert!(StructuralOperator::new(&mesh, Field::Stiffness).apply(&[], &mut out).is_err());
        let c = contributions(&mesh).find(|c| matches!(c.identity(), Identity::PointMass(_))).unwrap();
        let mut y = vec![0.; c.support().unwrap().len()];
        c.apply_local(Field::Mass(MassKind::Consistent), &[0., 0., 1., 2., 3.], &mut y).unwrap();
        close(&y, &[0., 0., 4., 11.5, 18.6]);
    }

    #[test]
    fn reusable_scratch_and_hostile_discrete_indices() {
        let mut mesh = fixture(false);
        let x = vec![0.25; mesh.dofs()];
        let mut expected = vec![0.; mesh.dofs()];
        let mut actual = vec![9.; mesh.dofs()];
        let mut scratch = OperatorScratch::default();
        let op = StructuralOperator::new(&mesh, Field::Mass(MassKind::Consistent));
        op.apply(&x, &mut expected).unwrap();
        for _ in 0..2 {
            op.apply_with_scratch(&x, &mut actual, &mut scratch).unwrap();
            close(&actual, &expected);
        }
        mesh.discrete.masses[0].rotations.push(usize::MAX);
        let c = contributions(&mesh).find(|c| matches!(c.identity(), Identity::PointMass(_))).unwrap();
        assert!(c.local_matrix(Field::Mass(MassKind::Consistent)).is_err());
        assert!(StructuralOperator::new(&mesh, Field::Mass(MassKind::Consistent)).apply_with_scratch(&x, &mut actual, &mut scratch).is_err());
        assert!(Contribution::square(usize::MAX).is_err());
    }
}
