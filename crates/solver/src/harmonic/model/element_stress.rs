//! Element-local elastic recovery from solved phasors. Extrapolation uses the
//! same FE integration rule as assembly; contributions are never averaged across elements.
use super::*;
use std::collections::BTreeSet;

#[derive(Debug, Clone, Serialize)]
pub struct StressLocation {
    pub element: u32,
    pub position: &'static str,
    /// One-based local point number, not a vendor integration-point numbering claim.
    pub point: usize,
    pub node: Option<u32>,
    pub xyz: [f64; 3],
    pub natural_coordinates: Option<[f64; 3]>,
}

/// A bounded recovery plan for solid elements of the current model.
pub struct ElementStressRecovery {
    elements: Vec<system::Elem>,
    indices: Vec<Vec<usize>>,
    pub locations: Vec<StressLocation>,
}

impl ElementStressRecovery {
    pub fn new(model: &Model, h: &HarmonicResult, labels: &[u32], limits: &Limits) -> Result<Self> {
        if labels.is_empty() || labels.len() > 1024 || labels.iter().collect::<BTreeSet<_>>().len() != labels.len() {
            return Err(Error::Invalid("element stress requires 1..1024 unique element labels".into()));
        }
        let resolved = model.resolve()?;
        let mesh = Mesh::new(model, &resolved, limits.max_dofs)?;
        let mut elements = Vec::new();
        let mut indices = Vec::new();
        let mut locations = Vec::new();
        for label in labels {
            let e = mesh
                .elements
                .iter()
                .find(|e| e.label == *label && e.theory.is_none())
                .ok_or_else(|| Error::NotSupported(format!("element {label}: element-local harmonic stress currently requires a solid")))?;
            let kernel = e.kernel()?;
            let mut ix = Vec::new();
            for (point, &n) in e.nodes.iter().enumerate() {
                let node = *mesh.labels.get(n).ok_or_else(|| Error::Internal("stress node label".into()))?;
                ix.push(h.nodes.binary_search(&node).map_err(|_| Error::Invalid(format!("node {node} has no harmonic phasor")))?);
                locations.push(StressLocation {
                    element: *label,
                    position: "element_nodes",
                    point: point + 1,
                    node: Some(node),
                    xyz: mesh.xyz[n],
                    natural_coordinates: None,
                });
            }
            let system::Kernel::Solid(solid) = kernel else {
                return Err(Error::NotSupported("solid stress recovery only".into()));
            };
            for (point, q) in solid.integration_points().iter().enumerate() {
                locations.push(StressLocation {
                    element: *label,
                    position: "integration_points",
                    point: point + 1,
                    node: None,
                    xyz: solid.position(q.xi),
                    natural_coordinates: Some(q.xi),
                });
            }
            elements.push(e.clone());
            indices.push(ix);
        }
        if locations.len().checked_mul(h.frames.len()).is_none_or(|n| n > 1_000_000) {
            return Err(Error::Limit("element stress exceeds 1000000 location-frame samples; select fewer elements".into()));
        }
        Ok(Self { elements, indices, locations })
    }

    /// Global Cartesian components [xx, yy, zz, xy, xz, yz], in Pa.
    pub fn recover(&self, displacement: &[[f64; 3]]) -> Result<Vec<[f64; 6]>> {
        let mut out = Vec::new();
        for (e, ix) in self.elements.iter().zip(&self.indices) {
            fieldcraft_model::progress::check()?;
            let ue = ix
                .iter()
                .map(|&i| displacement.get(i).copied().ok_or_else(|| Error::Invalid("incomplete harmonic displacement field".into())))
                .collect::<Result<Vec<_>>>()?;
            let values = e.stresses_from(&ue, None, None)?.iter().map(|p| p.stress).collect::<Vec<_>>();
            let kernel = e.kernel()?;
            let nodal = extrapolation(kernel.ty(), kernel.rule()).map_err(|err| Error::element(e.label, err))?.apply(&values);
            out.extend(nodal);
            out.extend(values);
        }
        Ok(out)
    }
}
