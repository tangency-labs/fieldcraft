//! Fixed-seed, prespecified 200-replicate local joint-ellipsoid coverage check.
//! Consistent-mass 100-element rod. Independent Gaussian errors are injected in
//! static tip displacement and first frequency; no covariance from noisy ratios.
use fieldcraft_model::{
    design::{MaterialConstant, Parameter},
    *,
};
use fieldcraft_solver::{Limits, sensitivity::*, solve_frequency, solve_static};
fn rod() -> Model {
    let mut m = Model::default();
    let cells = 100;
    for i in 0..=cells {
        m.set_node(i + 1, [i as f64 / cells as f64, 0.0, 0.0]).unwrap();
    }
    for i in 0..cells {
        m.set_element(i + 1, ElementKind::Rod2, vec![i + 1, i + 2]).unwrap();
    }
    for (name, nodes, elements) in
        [("all", (1..=cells + 1).collect(), (1..=cells).collect()), ("root", vec![1], vec![]), ("tip", vec![cells + 1], vec![])]
    {
        m.add_set(name, Set { nodes, elements, ..Default::default() }).unwrap();
    }
    m.add_material("bar", Material { elastic: Some(Elastic::Isotropic { young: 1000.0, poisson: 0.25 }), density: Some(2.0), ..Default::default() })
        .unwrap();
    m.add_section(
        "bar",
        Section {
            set: "all".into(),
            material: "bar".into(),
            kind: SectionKind::Rod { area: 1.0 },
            formulation: Default::default(),
            orientation: None,
        },
    )
    .unwrap();
    m.add_boundary_condition("root", BoundaryCondition::Fixed { target: Target::Set("root".into()) }).unwrap();
    m.add_boundary_condition(
        "guide",
        BoundaryCondition::Displacement { target: Target::Set("all".into()), values: [None, Some(0.0), Some(0.0), None, None, None] },
    )
    .unwrap();
    m.add_load("force", Load::Force { target: Target::Set("tip".into()), vector: [1.0, 0.0, 0.0], distribution: Distribution::Total }).unwrap();
    for (name, procedure) in [("S", Procedure::Static {}), ("F", Procedure::Frequency { modes: 1, min_hz: None, max_hz: None })] {
        let mut s = Step::new(name, procedure);
        s.boundary_conditions.extend(["root".into(), "guide".into()]);
        if name == "S" {
            s.loads.insert("force".into());
        }
        m.add_step(s, None).unwrap();
    }
    m
}
#[test]
#[ignore = "heavy: 200 prespecified FE-backed noise replicates"]
fn joint_ellipsoid_coverage_200() {
    let m = rod();
    let limits = Limits::default();
    let parameters = vec![
        JointParameter {
            parameter: Parameter::Material { name: "bar".into(), constant: MaterialConstant::Young },
            lower_si: 800.0,
            upper_si: 1300.0,
            scale_si: 1000.0,
        },
        JointParameter {
            parameter: Parameter::Material { name: "bar".into(), constant: MaterialConstant::Density },
            lower_si: 1.5,
            upper_si: 2.5,
            scale_si: 2.0,
        },
    ];
    let truth = [1100.0, 2.1];
    let mut target = m.clone();
    for (p, x) in parameters.iter().zip(truth) {
        p.parameter.set(&mut target, x).unwrap();
    }
    let u = solve_static(&target, Some("S"), &limits).unwrap().displacement[100][0];
    let f = solve_frequency(&target, Some("F"), &Default::default(), &limits).unwrap().modes[0].frequency;
    let mut fit = JointFit {
        parameters,
        observations: vec![
            JointObservation {
                response: JointResponse::Static { step: "S".into(), node: 101, dof: 1 },
                value_si: u,
                standard_deviation_si: u * 0.002,
            },
            JointObservation { response: JointResponse::Frequency { step: "F".into(), mode: 0 }, value_si: f, standard_deviation_si: f * 0.002 },
        ],
        regularization: 0.0,
        max_iterations: 32,
        relative_increment: 1e-5,
    };
    // LCG uniforms, fixed before observing any coverage result; Box-Muller normals.
    let mut seed = 0x2511_0010_2026_1107_u64;
    let mut uniform = || {
        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        ((seed >> 11) as f64 + 0.5) / 9007199254740992.0
    };
    let n = 200;
    let mut covered = 0;
    for _ in 0..n {
        let radius = (-2.0 * uniform().ln()).sqrt();
        let angle = std::f64::consts::TAU * uniform();
        fit.observations[0].value_si = u + fit.observations[0].standard_deviation_si * radius * angle.cos();
        fit.observations[1].value_si = f + fit.observations[1].standard_deviation_si * radius * angle.sin();
        let r = inverse_joint_uniform_material(&m, &fit, &limits).unwrap();
        let c = r.covariance_si.as_ref().unwrap();
        let d = [r.estimate_si[0] - truth[0], r.estimate_si[1] - truth[1]];
        let determinant = c[0][0] * c[1][1] - c[0][1] * c[1][0];
        let distance = (c[1][1] * d[0] * d[0] - 2.0 * c[0][1] * d[0] * d[1] + c[0][0] * d[1] * d[1]) / determinant;
        if distance <= r.joint_chi_square {
            covered += 1;
        }
    }
    // Exact binomial central quantiles, recurrence from k=0. No fitted acceptance band.
    let mut probability = 0.05_f64.powi(n);
    let mut cdf = 0.0;
    let mut lower = None;
    let mut upper = None;
    for k in 0..=n {
        cdf += probability;
        if lower.is_none() && cdf >= 0.025 {
            lower = Some(k);
        }
        if upper.is_none() && cdf >= 0.975 {
            upper = Some(k);
        }
        if k < n {
            probability *= ((n - k) as f64 / (k + 1) as f64) * 19.0;
        }
    }
    let lower = lower.unwrap();
    let upper = upper.unwrap();
    println!("coverage seed=0x2511001020261107 n={n} covered={covered} binomial95=[{lower},{upper}] cells=100 sigma_relative=0.002");
    assert!((lower..=upper).contains(&covered));
}
