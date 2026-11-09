//! Same linear FE operators as harmonic response, before factorization; native Tᵀ A T reduction.
use super::*;
use fieldcraft_model::steps::MatrixKind;

pub struct MatrixExport {
    pub step: String,
    /// Independent coordinates in equation order: native model node label and 1-based DOF.
    pub degrees: Vec<(u32, u8)>,
    pub matrices: BTreeMap<MatrixKind, Csr>,
    pub warnings: Vec<String>,
}

pub fn export_matrices(model: &Model, name: Option<&str>, limits: &Limits) -> Result<MatrixExport> {
    fieldcraft_model::progress::check()?;
    let resolved = model.resolve()?;
    let step = model.select_step(name, |p| matches!(p, Procedure::MatrixGeneration { .. }), "matrix generation")?;
    let Procedure::MatrixGeneration { matrices } = &step.procedure else {
        return Err(Error::Internal("matrix step".into()));
    };
    system::refuse_contact(model, "matrix generation")?;
    if model
        .steps()
        .iter()
        .take_while(|s| s.name != step.name)
        .any(|s| matches!(s.procedure, Procedure::Static {} | Procedure::NonlinearStatic(_) | Procedure::ExplicitDynamic(_)))
    {
        return Err(Error::NotSupported("matrix generation after a general load history needs an explicit tangent/prestress mapping".into()));
    }
    if !step.loads.is_empty() || fieldcraft_model::initial::left_out(model, step).is_some() {
        return Err(Error::NotSupported("matrix export with loads or initial fields".into()));
    }
    if !system::inelastic_warnings(model).is_empty() {
        return Err(Error::NotSupported("matrix export requires linear elastic material behavior".into()));
    }
    for (name, section) in model.sections() {
        if let SectionKind::CompositeShell { layup } = &section.kind
            && layup.plies.iter().any(|p| model.materials().get(&p.material).and_then(|m| m.damping).is_some_and(|d| d != Damping::default()))
        {
            return Err(Error::NotSupported(format!("section `{name}`: per-ply matrix damping")));
        }
    }
    let mesh = Mesh::new(model, &resolved, limits.max_dofs)?;
    let mut warnings: Vec<_> = resolved.warnings.iter().map(ToString::to_string).chain(mesh.warnings.iter().cloned()).collect();
    let (graph, order) = analysis::pattern(&mesh);
    let setup = Setup { model, resolved: &resolved, mesh: &mesh, step };
    let red = analysis::reduction(&setup, true, &order, analysis::budget(limits, &mesh, &graph, 4), &mut warnings)?;
    if red.neq() == 0 || red.largest_prescribed() != 0.0 {
        return Err(Error::NotSupported("matrix export needs free coordinates and homogeneous prescribed motion".into()));
    }
    let assembled = analysis::assemble(&mesh, &graph, matrices.contains(&MatrixKind::Mass).then_some(MassKind::Consistent))?;
    let (c, g) = if matrices.contains(&MatrixKind::ViscousDamping) || matrices.contains(&MatrixKind::StructuralDamping) {
        damping_matrices(&mesh, &graph, &material_damping(model, &resolved), MassKind::Consistent)?
    } else {
        (None, None)
    };
    let mut full = BTreeMap::new();
    if matrices.contains(&MatrixKind::Stiffness) {
        full.insert(MatrixKind::Stiffness, assembled.k);
    }
    if matrices.contains(&MatrixKind::Mass) {
        full.insert(MatrixKind::Mass, assembled.m.ok_or_else(|| Error::Internal("missing matrix mass".into()))?);
    }
    for (kind, value) in [(MatrixKind::ViscousDamping, c), (MatrixKind::StructuralDamping, g)] {
        if matrices.contains(&kind) {
            full.insert(kind, value.unwrap_or_else(|| graph.zero_matrix(&mesh.layout)));
        }
    }
    let mut reduced = BTreeMap::new();
    let mut entries = 0usize;
    for (kind, mut a) in full {
        a.mirror_lower();
        let [mut r] = red.reduce([&a])?;
        r.mirror_lower();
        entries = entries.saturating_add(r.val.len());
        if entries > limits.max_factor_entries {
            return Err(Error::Limit("exported matrix storage exceeds entry budget".into()));
        }
        if r.val.iter().any(|v| !v.is_finite()) {
            return Err(Error::Invalid("nonfinite matrix coefficient".into()));
        }
        reduced.insert(kind, r);
    }
    let degrees = red
        .eq_dof
        .iter()
        .map(|&dof| {
            let (slot, kind) = mesh.layout.locate(dof).ok_or_else(|| Error::Internal("missing matrix coordinate".into()))?;
            if kind.direction().is_none() {
                return Err(Error::NotSupported("nonstructural matrix coordinate".into()));
            }
            let label = mesh
                .labels
                .get(slot)
                .copied()
                .or_else(|| slot.checked_sub(mesh.labels.len()).and_then(|i| mesh.points.get(i)).map(|p| p.label))
                .ok_or_else(|| Error::Internal("missing matrix node label".into()))?;
            Ok((label, kind.number()))
        })
        .collect::<Result<Vec<_>>>()?;
    Ok(MatrixExport { step: step.name.clone(), degrees, matrices: reduced, warnings })
}
