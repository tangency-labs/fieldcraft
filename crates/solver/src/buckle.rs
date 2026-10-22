//! Linear buckling: the load factors λ at which a step's loads, applied as a linear preload,
//! make the structure lose stiffness, (K + λ K_σ) φ = 0, and the mode shapes φ (K.-J. Bathe,
//! *Finite Element Procedures*, 2nd ed., 2014, §6.8.2).
//!
//! **Preload.** The step is solved as a linear static step (its loads and prescribed motions);
//! the stresses at every element's stiffness points give the initial-stress stiffness
//! K_σ = ∫ (∇N_a · σ ∇N_b) I dV ([`fieldcraft_fem::finite`]). The eigenproblem holds every
//! boundary condition of the step at zero (the buckling motion is a perturbation of the
//! preloaded state). Pressures act as dead loads: the load stiffness of a follower pressure is
//! not part of K_σ.
//!
//! **Eigenproblem** by the block Lanczos of [`crate::eigen`], in the K inner product (K is
//! positive definite for a supported structure, while K_σ is indefinite), after R. G. Grimes,
//! J. G. Lewis, H. D. Simon, *A shifted block Lanczos algorithm for solving sparse symmetric
//! generalized eigenproblems*, SIAM J. Matrix Anal. Appl. 15(1):228–272, 1994, §3 (the buckling
//! transformations):
//! - **without a shift**, the operator K⁻¹(−K_σ) has eigenvalues θ = 1/λ: the lowest positive
//!   load factors are its largest positive θ, well separated, and K is factored by Cholesky.
//!   Converged negative θ are the load factors of the reversed load (reported, not returned as
//!   modes);
//! - **with a shift σ ≠ 0**, the buckling spectral transformation (K + σ K_σ)⁻¹ K has
//!   θ = λ/(λ − σ), largest for the load factors nearest above σ; K + σ K_σ is indefinite once σ
//!   passes the first load factor and is factored by LDLᵀ ([`crate::ldlt`]), whose inertia counts
//!   the load factors between 0 and σ (Sylvester's law), reported so none is missed silently.
//!
//! Modes are scaled to a largest component of 1 (positive). Solids and 2D elements only:
//! trusses, shells and beams have no initial-stress stiffness yet and are refused by name.

use std::collections::BTreeMap;

use fieldcraft_model::{Model, Procedure};
use serde::Serialize;

use crate::analysis::{
    BLOCK, FrequencyResult, Mode, Setup, Stats, TOLERANCE, assemble, budget, factor_error, factor_symmetric, lanczos_basis, pattern, reduction,
    split, static_on,
};
use crate::eigen::{Operator, Wanted, lanczos};
use crate::ldlt::{Kind, Symmetric};
use crate::sparse::{Csr, add_block};
use crate::system::Mesh;
use crate::{Error, Limits, Result};

/// What a buckling eigenproblem found besides its modes.
#[derive(Debug, Clone, Default, PartialEq, Serialize, serde::Deserialize)]
pub struct Buckling {
    /// Converged negative load factors (the load reversed buckles at |λ|), nearest to zero
    /// first; found without a shift only.
    pub negative: Vec<f64>,
    /// The shift, and the number of load factors between 0 and it (the inertia of K + σ K_σ).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub shift: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub below_shift: Option<usize>,
}

/// K⁻¹(−K_σ) in the K inner product: θ = 1/λ.
struct Plain<'a> {
    k: &'a Csr,
    /// −K_σ.
    g: &'a Csr,
    factor: &'a Symmetric,
}

impl Operator for Plain<'_> {
    fn n(&self) -> usize {
        self.k.n
    }
    fn mass(&self, x: &[f64], y: &mut [f64]) {
        self.k.mul(x, y);
    }
    fn solve(&self, x: &mut [f64]) -> Result<()> {
        self.factor.solve(x)
    }
    fn operand(&self, x: &[f64], _kx: &[f64]) -> Vec<f64> {
        let mut y = vec![0.0; self.g.n];
        self.g.mul(x, &mut y);
        y
    }
    fn eigenvalue(&self, _shift: f64, theta: f64) -> f64 {
        1.0 / theta
    }
}

/// (K + σ K_σ)⁻¹ K in the K inner product: θ = λ/(λ − σ).
struct Shifted<'a> {
    k: &'a Csr,
    /// K + σ K_σ and its factor, for iterative refinement.
    shifted: &'a Csr,
    factor: &'a Symmetric,
}

impl Operator for Shifted<'_> {
    fn n(&self) -> usize {
        self.k.n
    }
    fn mass(&self, x: &[f64], y: &mut [f64]) {
        self.k.mul(x, y);
    }
    fn solve(&self, x: &mut [f64]) -> Result<()> {
        self.factor.solve_refined(self.shifted, x, 2)
    }
    fn eigenvalue(&self, shift: f64, theta: f64) -> f64 {
        shift * theta / (theta - 1.0)
    }
}

/// Solve a buckling step (the step named `step`, or the model's only buckling step).
pub fn solve_buckle(model: &Model, step: Option<&str>, limits: &Limits) -> Result<FrequencyResult> {
    // A job cancelled before it starts stops here, before any setup work.
    fieldcraft_model::progress::check()?;
    let step = model.select_step(step, |p| matches!(p, Procedure::Buckle { .. }), "buckling")?;
    let Procedure::Buckle { modes, shift } = step.procedure else {
        return Err(Error::Internal("buckling step".into()));
    };
    step.procedure.validate()?;
    let resolved = model.resolve()?;
    crate::system::refuse_contact(model, "buckling")?;
    let mesh = Mesh::new(model, &resolved, limits.max_dofs)?;
    if let Some(t) = mesh.trusses.first() {
        return Err(Error::NotSupported(format!(
            "trusses in buckling steps (element {}): their initial-stress stiffness is not implemented",
            t.label
        )));
    }
    if let Some(s) = mesh.structural.first() {
        return Err(Error::NotSupported(format!(
            "shells and beams in buckling steps (element {}): their initial-stress stiffness is not implemented",
            s.label
        )));
    }
    // The preload: the step solved as a linear static step.
    let setup = Setup { model, resolved: &resolved, mesh: &mesh, step };
    let preload = static_on(&setup, limits)?;
    let mut warnings: Vec<String> = preload.warnings.clone();
    let stresses: BTreeMap<u32, Vec<[f64; 6]>> = preload.elements.iter().map(|e| (e.label, e.points.iter().map(|p| p.stress).collect())).collect();

    let (graph, order) = pattern(&mesh);
    let red = reduction(&setup, false, &order, budget(limits, &mesh, &graph, 2), &mut Vec::new())?;
    let neq = red.neq();
    if neq == 0 {
        return Err(Error::Invalid("every degree of freedom is held; there is nothing to buckle".into()));
    }
    let k = assemble(&mesh, &graph, None)?.k;
    let mut kg = graph.zero_matrix(&mesh.layout);
    let matrices = fieldcraft_parallel::map_grained(mesh.elements.len(), crate::analysis::GRAIN, |i| -> Result<_> {
        let e = mesh.elements.get(i).ok_or_else(|| Error::Internal("element index".into()))?;
        let kernel = e.kernel()?;
        let stress = stresses.get(&e.label).ok_or_else(|| Error::Internal(format!("element {} has no preload stress", e.label)))?;
        Ok((kernel.dofs().node, kernel.geometric_stiffness(stress).map_err(|err| Error::element(e.label, err))?.data))
    });
    for (e, m) in mesh.elements.iter().zip(matrices) {
        let (kinds, data) = m?;
        add_block(&mut kg, &mesh.layout, &e.nodes, kinds, &data)?;
    }
    let [kff, gff] = red.reduce([&k, &kg])?;
    drop((k, kg));
    if gff.val.iter().all(|v| *v == 0.0) {
        return Err(Error::Invalid(format!("step `{}` puts no stress into the model, so nothing can buckle (give it loads)", step.name)));
    }
    let max_basis = lanczos_basis(limits, neq);
    let count = usize::try_from(modes).unwrap_or(usize::MAX).min(neq);
    let sigma = shift.unwrap_or(0.0);
    let (pairs, entries, buckling) = if sigma == 0.0 {
        let factor = factor_symmetric(&kff, limits, Kind::Definite).map_err(|e| factor_error(e, &mesh, &red))?;
        let minus = Csr { val: gff.val.iter().map(|v| -v).collect(), ..gff.clone() };
        let op = Plain { k: &kff, g: &minus, factor: &factor };
        let pairs = lanczos(&op, 0.0, Wanted { count, min: 0.0, max: f64::INFINITY }, BLOCK, max_basis, TOLERANCE)?;
        let buckling = Buckling { negative: pairs.negative.clone(), shift: None, below_shift: None };
        (pairs, factor.entries(), buckling)
    } else {
        let shifted = kff.add_scaled(sigma, &gff)?;
        let factor = factor_symmetric(&shifted, limits, Kind::Indefinite).map_err(|e| match e {
            crate::cholesky::FactorError::Singular { .. } => {
                Error::Singular(format!("the shift {sigma} is a buckling load factor (K + σ K_σ is singular); move it a little"))
            }
            e => factor_error(e, &mesh, &red),
        })?;
        let below = factor.inertia().negative;
        if below > 0 {
            warnings.push(format!("{below} buckling load factor(s) lie between 0 and the shift {sigma}; the modes returned are those above it"));
        }
        let op = Shifted { k: &kff, shifted: &shifted, factor: &factor };
        let (min, max) = if sigma > 0.0 { (sigma, f64::INFINITY) } else { (f64::NEG_INFINITY, f64::INFINITY) };
        let pairs = lanczos(&op, sigma, Wanted { count, min, max }, BLOCK, max_basis, TOLERANCE)?;
        (pairs, factor.entries(), Buckling { negative: Vec::new(), shift: Some(sigma), below_shift: Some(below) })
    };
    if let Some(first) = buckling.negative.first() {
        warnings.push(format!("the reversed load buckles too: negative load factor {first:.6e} (and {} more found)", buckling.negative.len() - 1));
    }
    let mut out = Vec::with_capacity(pairs.values.len());
    for (lambda, phi) in pairs.values.iter().zip(&pairs.vectors) {
        // Scale to a largest component of 1 (Lanczos already made it positive).
        let big = phi.iter().fold(0.0f64, |m, v| m.max(v.abs()));
        let scaled: Vec<f64> = phi.iter().map(|v| if big > 0.0 { v / big } else { *v }).collect();
        let (shape, rotation, point_shape) = split(&red.expand(&scaled, false), &mesh);
        out.push(Mode { eigenvalue: *lambda, frequency: 0.0, shape, rotation, participation: [0.0; 3], effective_mass: [0.0; 3], point_shape });
    }
    if out.len() < count {
        warnings.push(format!("found {} of the {modes} buckling modes asked for", out.len()));
    }
    Ok(FrequencyResult {
        step: step.name.clone(),
        nodes: mesh.labels.clone(),
        reference_points: mesh.points.iter().map(|p| (p.name.clone(), p.label)).collect(),
        modes: out,
        total_mass: 0.0,
        stats: Stats { dofs: mesh.dofs(), equations: neq, factor_entries: entries, lanczos_vectors: pairs.basis },
        warnings,
        buckling: Some(buckling),
    })
}
