//! Cross-feature regression on a two-cell guided elastic-foundation bar.
//! Central response differences use h/p=1e-4; tolerance 2e-5 relative +1e-9
//! absolute allows truncation/eigenvector normalization without hiding zero columns.
#[path = "../../solver/tests/common/mod.rs"]
mod common;
use common::*;
use fieldcraft_engine::Session;
use fieldcraft_model::{
    design::{DampingParameter, FoundationComponent, MaterialConstant, Parameter},
    harmonic::{Damping, Frequencies, Harmonic, Method},
    *,
};
use fieldcraft_solver::{FrequencyOptions, Limits, harmonic::solve_harmonic, sensitivity::*, solve_frequency};
use serde_json::json;

fn system() -> (Session, u32, u32) {
    let mut g = grid([2, 1, 1], ElementKind::Hex8, |p| p);
    g.define(iso(1000., 0.), Some(2.), full(), Procedure::Static {});
    g.face_set("bed", |p| near(p[0], 0.));
    g.face_set("tip", |p| near(p[0], 1.));
    g.model
        .add_foundation(
            "Bed",
            Foundation { set: "bed".into(), normal: 100., tangential: 100., normal_damping: 0., tangential_damping: 0., loss_factor: 0. },
        )
        .unwrap();
    g.bc("guide", BoundaryCondition::Displacement { target: Target::Set("all".into()), values: [None, Some(0.), Some(0.), None, None, None] });
    g.load("force", Load::Force { target: Target::Set("tip".into()), vector: [1., 0., 0.], distribution: Distribution::Total });
    let mut f = Step::new("F", Procedure::Frequency { modes: 1, min_hz: None, max_hz: None });
    f.boundary_conditions.insert("guide".into());
    g.model.add_step(f, None).unwrap();
    let mut h = Step::new(
        "H",
        Procedure::Harmonic(Harmonic {
            frequencies: Frequencies::List(vec![0.2, 0.6, 1.]),
            method: Method::Direct,
            damping: Damping { alpha: 0.2, beta: 0.003, loss_factor: 0.02 },
            load_factors: Default::default(),
        }),
    );
    h.boundary_conditions.insert("guide".into());
    h.loads.insert("force".into());
    g.model.add_step(h, None).unwrap();
    let root = g.nodes_where(|p| near(p[0], 0.))[0];
    let tip = g.nodes_where(|p| near(p[0], 1.))[0];
    let mut s = Session::default();
    s.model = g.model;
    (s, root, tip)
}
fn bed() -> Parameter {
    Parameter::Foundation { name: "Bed".into(), component: FoundationComponent::Normal }
}
fn damping() -> Parameter {
    Parameter::Damping { step: "H".into(), component: DampingParameter::Alpha }
}
fn close(a: f64, b: f64) {
    assert!((a - b).abs() <= 1e-9 + 2e-5 * b.abs(), "{a} != {b}");
}
fn modal(m: &Model, root: u32, tip: u32) -> (f64, f64, Vec<f64>) {
    let r = solve_frequency(m, Some("F"), &FrequencyOptions::default(), &Limits::default()).unwrap();
    let v = [root, tip].map(|n| r.modes[0].shape[r.nodes.binary_search(&n).unwrap()][0]);
    (r.modes[0].frequency, v[1] / v[0], v.to_vec())
}
#[test]
fn foundation_mac_frequency_and_nelson_component_match_differences() {
    let (mut s, root, tip) = system();
    let (f, ratio, shape) = modal(&s.model, root, tip);
    let r=s.execute("solve.inverse.joint",&json!({"parameters":[
        {"parameter":bed(),"lower_si":50,"upper_si":200,"scale_si":100},
        {"parameter":Parameter::Material{name:"Mat".into(),constant:MaterialConstant::Young},"lower_si":500,"upper_si":2000,"scale_si":1000}],
        "observations":[{"response":{"kind":"matched_frequency","group":"g","measured_mode":0},"value_si":f,"standard_deviation_si":0.001},{"response":{"kind":"mode_shape","group":"g","measured_mode":0,"sensor":1},"value_si":ratio,"standard_deviation_si":0.001}],
        "modal_observations":[{"id":"g","step":"F","sensors":[{"node":root,"dof":1},{"node":tip,"dof":1}],"measured_modes":[shape],"min_mac":0.99,"min_assignment_margin":0.001,"normalization":{"kind":"reference","reference_sensor":0}}]})).unwrap();
    assert!(r["modal_pairings"][0]["pairing"]["pairs"][0]["mac"].as_f64().unwrap() > 0.999999);
    let mut plus = s.model.clone();
    let mut minus = s.model.clone();
    bed().set(&mut plus, 100.01).unwrap();
    bed().set(&mut minus, 99.99).unwrap();
    let a = modal(&plus, root, tip);
    let b = modal(&minus, root, tip);
    for (row, fd) in [(0, (a.0 - b.0) / 0.02), (1, (a.1 - b.1) / 0.02)] {
        assert!(fd.abs() > 1e-8);
        close(r["jacobian_si"][row][0].as_f64().unwrap(), fd);
    }
}
#[test]
fn damping_and_foundation_share_a_harmonic_derivative_batch() {
    let (mut s, _, tip) = system();
    for method in [Method::Direct, Method::Modal { modes: 12, damping_ratios: vec![], residual_vectors: 0 }] {
        let mut step = s.model.step("H").unwrap().clone();
        if let Procedure::Harmonic(h) = &mut step.procedure {
            h.method = method;
        }
        s.model.replace_step(step).unwrap();
        let parameters = [damping(), bed()];
        let increments = [0.00002, 0.01];
        let requests = parameters
            .iter()
            .zip(increments)
            .map(|(p, h)| HarmonicSensitivityRequest {
                step: Some("H".into()),
                parameter: p.clone(),
                increment: h,
                observation: HarmonicObservation { node: tip, dof: 1 },
            })
            .collect::<Vec<_>>();
        let batch = harmonic_sensitivity_batch(&s.model, &requests, &Limits::default()).unwrap();
        assert_eq!(batch.results.len(), 2);
        for (col, p) in parameters.iter().enumerate() {
            let value = p.value(&s.model).unwrap();
            let mut plus = s.model.clone();
            let mut minus = s.model.clone();
            p.set(&mut plus, value + increments[col]).unwrap();
            p.set(&mut minus, value - increments[col]).unwrap();
            let a = solve_harmonic(&plus, Some("H"), &Limits::default()).unwrap();
            let b = solve_harmonic(&minus, Some("H"), &Limits::default()).unwrap();
            let i = a.nodes.binary_search(&tip).unwrap();
            for (j, (a, b)) in a.frames.iter().zip(&b.frames).enumerate() {
                for (k, fd) in [
                    (0, (a.real[i][0] - b.real[i][0]) / (2. * increments[col])),
                    (1, (a.imaginary[i][0] - b.imaginary[i][0]) / (2. * increments[col])),
                ] {
                    close(batch.results[col].frames[j].direct[k], fd);
                }
            }
        }
    }
}
#[test]
fn damping_static_observation_is_zero_in_joint_fit() {
    let (mut s, _, tip) = system();
    let stat = fieldcraft_solver::solve_static(&s.model, Some("Step"), &Limits::default()).unwrap();
    let u = stat.displacement[stat.nodes.binary_search(&tip).unwrap()][0];
    let h = solve_harmonic(&s.model, Some("H"), &Limits::default()).unwrap();
    let i = h.nodes.binary_search(&tip).unwrap();
    let r=s.execute("solve.inverse.joint",&json!({"parameters":[{"parameter":damping(),"lower_si":0.1,"upper_si":0.4,"scale_si":0.2},{"parameter":bed(),"lower_si":50,"upper_si":200,"scale_si":100}],"observations":[{"response":{"kind":"static","step":"Step","node":tip,"dof":1},"value_si":u,"standard_deviation_si":1e-5},{"response":{"kind":"harmonic_imaginary","step":"H","node":tip,"dof":1,"frame":1},"value_si":h.frames[1].imaginary[i][0],"standard_deviation_si":1e-5}]})).unwrap();
    assert_eq!(r["jacobian_si"][0][0].as_f64().unwrap(), 0.);
}
#[test]
fn foundation_force_random_matches_frame_integral_and_base_is_explicit() {
    let (mut s, _, tip) = system();
    let h = solve_harmonic(&s.model, Some("H"), &Limits::default()).unwrap();
    let i = h.nodes.binary_search(&tip).unwrap();
    s.execute("solve.harmonic", &json!({"step":"H","threads":1})).unwrap();
    let r = s.execute("solve.random", &json!({"step":"H","input_psd":[[0.2,2.],[1.,4.]]})).unwrap();
    let psd = h
        .frames
        .iter()
        .enumerate()
        .map(|(j, f)| {
            let input = 2. + (f.frequency - 0.2) * 2. / 0.8;
            let value = (0..3).map(|k| f.real[i][k] * f.real[i][k] + f.imaginary[i][k] * f.imaginary[i][k]).sum::<f64>() * input;
            close(r["node_displacement_psd"][j][i].as_f64().unwrap(), value);
            value
        })
        .collect::<Vec<_>>();
    let variance = h.frames.windows(2).zip(psd.windows(2)).map(|(f, p)| (f[1].frequency - f[0].frequency) * (p[0] + p[1]) / 2.).sum::<f64>();
    close(r["node_displacement_rms"][i].as_f64().unwrap(), variance.sqrt());
    let random=serde_json::from_value(json!({"frequencies":[0.2,1.],"excitations":[{"type":"base_acceleration","dof":1,"spectrum":[{"frequency":0.2,"psd":2.},{"frequency":1.,"psd":4.}]}],"damping":{"alpha":0.2,"beta":0.,"loss_factor":0.}})).unwrap();
    s.model.add_step(Step::new("Base", Procedure::RandomVibration(random)), None).unwrap();
    let err = s.execute("solve.random", &json!({"step":"Base"})).unwrap_err().to_string();
    assert!(err.contains("foundation") && err.contains("base"), "{err}");
}
