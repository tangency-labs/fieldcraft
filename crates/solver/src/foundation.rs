//! Surface foundation assembly using the shared consistent FEM face kernel.
use crate::sparse::{Csr, add_block};
use crate::system::Mesh;
use crate::{Error, Result};
use fieldcraft_fem::DofSet;
use fieldcraft_model::{Model, Resolved};

#[derive(Debug, Clone)]
pub struct Face {
    pub name: String,
    pub nodes: Vec<usize>,
    pub xyz: Vec<[f64; 3]>,
    pub normal: f64,
    pub tangential: f64,
    pub normal_damping: f64,
    pub tangential_damping: f64,
    pub loss_factor: f64,
}
pub fn resolve(model: &Model, resolved: &Resolved, mesh: &Mesh) -> Result<Vec<Face>> {
    let mut out = Vec::new();
    for (name, f) in model.foundations() {
        if mesh.theory.is_some() {
            return Err(Error::NotSupported(format!("foundation `{name}` on a planar/axisymmetric model")));
        }
        let set = resolved.sets.get(&f.set).ok_or_else(|| Error::Invalid(format!("foundation `{name}` has no face set")))?;
        for face in &set.faces {
            let labels = model.face_nodes(*face).ok_or_else(|| Error::Invalid(format!("foundation `{name}` has an invalid face")))?;
            let nodes: Vec<_> = labels
                .iter()
                .map(|label| {
                    mesh.node(*label).ok_or_else(|| Error::Invalid(format!("foundation `{name}` node {label} has no structural degrees of freedom")))
                })
                .collect::<Result<_>>()?;
            let xyz: Vec<_> = labels
                .iter()
                .map(|label| model.label_xyz(*label).ok_or_else(|| Error::Invalid(format!("foundation `{name}` node {label} is absent"))))
                .collect::<Result<_>>()?;
            // Validate the reference geometry even when all coefficients are zero.
            fieldcraft_fem::foundation::matrix(&xyz, f.normal, f.tangential).map_err(|e| Error::Invalid(format!("foundation `{name}`: {e}")))?;
            out.push(Face {
                name: name.clone(),
                nodes,
                xyz,
                normal: f.normal,
                tangential: f.tangential,
                normal_damping: f.normal_damping,
                tangential_damping: f.tangential_damping,
                loss_factor: f.loss_factor,
            });
        }
    }
    Ok(out)
}
impl Face {
    pub fn add(&self, mesh: &Mesh, matrix: &mut Csr, normal: f64, tangential: f64) -> Result<()> {
        let data = fieldcraft_fem::foundation::matrix(&self.xyz, normal, tangential)
            .map_err(|e| Error::Invalid(format!("foundation `{}`: {e}", self.name)))?;
        add_block(matrix, &mesh.layout, &self.nodes, DofSet::TRANSLATIONS, &data)
    }
}

/// Exact dK/dk and dH/dk on the existing graph, including at k=0.
pub(crate) fn derivative(
    mesh: &Mesh,
    graph: &crate::sparse::NodeGraph,
    parameter: &fieldcraft_model::design::Parameter,
) -> Result<Option<(Csr, Csr)>> {
    use fieldcraft_model::design::{FoundationComponent, Parameter};
    let Parameter::Foundation { name, component } = parameter else {
        return Ok(None);
    };
    let mut k = graph.zero_matrix(&mesh.layout);
    let mut h = graph.zero_matrix(&mesh.layout);
    let (normal, tangential) = match component {
        FoundationComponent::Normal => (1., 0.),
        FoundationComponent::Tangential => (0., 1.),
    };
    for f in mesh.foundations.iter().filter(|f| &f.name == name) {
        f.add(mesh, &mut k, normal, tangential)?;
        f.add(mesh, &mut h, normal * f.loss_factor, tangential * f.loss_factor)?;
    }
    k.mirror_lower();
    h.mirror_lower();
    Ok(Some((k, h)))
}
