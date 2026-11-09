//! Exported operators must reconstruct frozen harmonic results without a new sweep.
use fieldcraft_engine::Session;
use fieldcraft_model::{MassKind, steps::MatrixKind};
use fieldcraft_solver::{
    Limits,
    complex::Vector,
    harmonic::{Matrices, Modal},
};
use serde_json::{Value, json};

#[test]
fn exported_shell_mass_reconstructs_frozen_response() {
    let recipe: Value = serde_json::from_str(include_str!("../../../tests/oracle/mitc8-navier-harmonic/recipe.json")).unwrap();
    let mut s = Session::default();
    for c in recipe.as_array().unwrap() {
        let command = c["command"].as_str().unwrap();
        if !command.starts_with("solve.") {
            s.execute(command, &c["params"]).unwrap();
        }
    }
    let original = recipe.as_array().unwrap().iter().find(|c| c["command"] == "model.step.create").unwrap();
    s.execute("model.step.create", &json!({"name":"Operators", "procedure":{"type":"matrix_generation", "matrices":["stiffness","mass","viscous_damping","structural_damping"]}, "boundary_conditions":original["params"]["boundary_conditions"]})).unwrap();
    let limits = Limits::default();
    for (name, mass, expected) in [
        ("consistent", MassKind::Consistent, [2.0559758370409166e-6, -2.5384524850928547e-8, 2.0559761816897116e-6, -2.5384524851292126e-8]),
        ("lumped", MassKind::Lumped, [2.0748815279850348e-6, -2.6051784086164216e-8, 2.075179723238668e-6, -2.6051784582294895e-8]),
    ] {
        let operators = fieldcraft_solver::export_matrices_with_mass(&s.model, Some("Operators"), mass, &limits).unwrap();
        let k = &operators.matrices[&MatrixKind::Stiffness];
        let m = &operators.matrices[&MatrixKind::Mass];
        let c = &operators.matrices[&MatrixKind::ViscousDamping];
        let g = &operators.matrices[&MatrixKind::StructuralDamping];
        assert!(c.val.iter().chain(&g.val).all(|v| *v == 0.0));
        let mesh = fieldcraft_solver::system::Mesh::new(&s.model, &s.model.resolve().unwrap(), limits.max_dofs).unwrap();
        let step = s.model.steps().iter().find(|p| p.name == "consistent_direct").unwrap();
        let full = fieldcraft_solver::system::loads(&s.model, &s.model.resolve().unwrap(), &mesh, step, &vec![true; mesh.dofs()]).unwrap();
        let force = Vector {
            real: operators
                .degrees
                .iter()
                .map(|&(node, dof)| full[mesh.dof(mesh.node(node).unwrap(), fieldcraft_fem::DofKind::from_number(dof).unwrap()).unwrap()])
                .collect(),
            imaginary: vec![0.0; k.n],
        };
        let matrices = Matrices { stiffness: k, mass: m, viscous: Some(c), loss: Some(g) };
        let modal = Modal::new(&matrices, 12, &limits).unwrap();
        let damping = fieldcraft_model::harmonic::Damping { alpha: 2.0, ..Default::default() };
        let omega = 60.0 * std::f64::consts::PI;
        let direct = matrices.direct(omega, damping, &force, &limits).unwrap();
        let reduced = modal.solve(omega, damping, &[], &force, &limits).unwrap();
        let probe = operators.degrees.iter().position(|x| *x == (33, 3)).unwrap();
        let values = [
            direct.displacement.real[probe],
            direct.displacement.imaginary[probe],
            reduced.displacement.real[probe],
            reduced.displacement.imaginary[probe],
        ];
        for (a, b) in values.into_iter().zip(expected) {
            assert!((a - b).abs() < b.abs() * 1e-9, "{name}: {a} vs {b}");
        }
        let mut projection = Vec::new();
        for shape in modal.shapes() {
            let mut mv = vec![0.0; k.n];
            m.mul(shape, &mut mv);
            let norm = shape.iter().zip(&mv).map(|(a, b)| a * b).sum::<f64>();
            assert!((norm - 1.0).abs() < 1e-10);
            projection.push(json!({"mass_norm":norm,"load":shape.iter().zip(&force.real).map(|(a,b)|a*b).sum::<f64>(),"probe":shape[probe]}));
        }
        let metadata = s.execute("solve.matrices", &json!({"step":"Operators","mass":name})).unwrap();
        assert_eq!(metadata["mass"], name);
        if let Ok(dir) = std::env::var("FIELDCRAFT_TEST_ARTIFACTS") {
            let path = std::path::Path::new(&dir).join(format!("{name}.json"));
            s.execute("solve.matrices", &json!({"step":"Operators","mass":name,"path":path})).unwrap();
            let record = json!({"degrees":operators.degrees,"force_real":force.real,"force_imaginary":force.imaginary,"probe_node":33,"probe_dof":3,"probe_equation":probe,"alpha":2.0,"beta":0.0,"loss_factor":0.0,"modal_damping_ratios":[],"eigenvalues":modal.eigenvalues(),"mode_shapes_equation_order":modal.shapes(),"projections":projection,"frequency_hz":30,"reconstructed":values,"frozen":expected});
            std::fs::write(std::path::Path::new(&dir).join(format!("{name}-projection.json")), serde_json::to_vec_pretty(&record).unwrap()).unwrap();
        }
    }
    assert!(s.execute("solve.matrices", &json!({"step":"Operators","mass":"unknown"})).is_err());
}

#[test]
fn material_rayleigh_export_uses_selected_mass() {
    let recipe: Value = serde_json::from_str(include_str!("../../../tests/oracle/solid-shell-pure-bending/recipe.json")).unwrap();
    let mut s = Session::default();
    for c in recipe.as_array().unwrap() {
        let cmd = c["command"].as_str().unwrap();
        if cmd.starts_with("solve.") || cmd == "model.step.create" {
            continue;
        }
        let mut params = c["params"].clone();
        if cmd == "model.material.create" {
            params["damping"] = json!({"alpha":0.7,"beta":0.0});
        }
        s.execute(cmd, &params).unwrap();
    }
    s.execute("model.step.create",&json!({"name":"Operators","procedure":{"type":"matrix_generation","matrices":["stiffness","mass","viscous_damping"]},"boundary_conditions":["clamp"]})).unwrap();
    for mass in [MassKind::Consistent, MassKind::Lumped] {
        let a = fieldcraft_solver::export_matrices_with_mass(&s.model, Some("Operators"), mass, &Limits::default()).unwrap();
        let k = &a.matrices[&MatrixKind::Stiffness];
        let m = &a.matrices[&MatrixKind::Mass];
        let c = &a.matrices[&MatrixKind::ViscousDamping];
        for ((&k, &m), &c) in k.val.iter().zip(&m.val).zip(&c.val) {
            let expected = 0.7 * m + 0.0 * k;
            assert!((c - expected).abs() <= 1e-12 * (1.0 + expected.abs()));
        }
    }
}
