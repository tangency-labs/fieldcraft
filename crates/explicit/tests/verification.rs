//! Analytic verification of the explicit solver (docs/solvers/literature.md, A10): elastic waves
//! in a bar (front at c·t, stress σ = ρ c v), a bar striking a rigid wall (contact for 2L/c,
//! rebound at the impact speed), the energy balance and bit-for-bit repeatability.
//!
//! The bar has Poisson's ratio 0, so a uniaxial stress wave travels at the rod speed c = √(E/ρ)
//! without lateral inertia effects (d'Alembert's solution of the one-dimensional wave equation).

use fieldcraft_explicit::{ExplicitResult, Limits, solve_explicit};
use fieldcraft_model::*;

const YOUNG: f64 = 1.0e4;
const DENSITY: f64 = 1.0;
const LENGTH: f64 = 10.0;
const WIDTH: f64 = 0.1;

fn wave_speed() -> f64 {
    (YOUNG / DENSITY).sqrt()
}

/// A bar along x from 0 to LENGTH, WIDTH × WIDTH, with `cells` cells along it: one hex8 per
/// cell, or six tet4 (Kuhn's split of the cube), with ν = 0; its sets, and an explicit step.
struct Bar {
    model: Model,
    kind: ElementKind,
    cells: usize,
}

impl Bar {
    fn new(kind: ElementKind, cells: usize, step: ExplicitDynamic) -> Bar {
        let mut m = Model::default();
        let h = LENGTH / cells as f64;
        let label = |i: usize, j: usize, k: usize| (1 + i + (cells + 1) * (j + 2 * k)) as u32;
        for k in 0..2 {
            for j in 0..2 {
                for i in 0..=cells {
                    m.set_node(label(i, j, k), [i as f64 * h, j as f64 * WIDTH, k as f64 * WIDTH]).unwrap();
                }
            }
        }
        let corners = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]];
        let mut next = 1;
        for i in 0..cells {
            let at = |c: [usize; 3]| label(i + c[0], c[1], c[2]);
            if kind == ElementKind::Hex8 {
                m.set_element(next, kind, corners.iter().map(|c| at(*c)).collect()).unwrap();
                next += 1;
                continue;
            }
            for perm in [[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]] {
                let mut v = [[0usize; 3]; 4];
                for s in 0..3 {
                    v[s + 1] = v[s];
                    v[s + 1][perm[s]] = 1;
                }
                let x = v.map(|c| [c[0] as f64 * h, c[1] as f64 * WIDTH, c[2] as f64 * WIDTH]);
                let e = |a: usize| [0, 1, 2].map(|k| x[a][k] - x[0][k]);
                let (a, b, c) = (e(1), e(2), e(3));
                let vol = a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0]);
                if vol < 0.0 {
                    v.swap(1, 2);
                }
                m.set_element(next, kind, v.iter().map(|c| at(*c)).collect()).unwrap();
                next += 1;
            }
        }
        let all: Vec<u32> = m.elements().keys().copied().collect();
        m.add_set("all", Set { elements: all, ..Set::default() }).unwrap();
        let end_nodes = |i: usize| (0..4).map(|n| label(i, n % 2, n / 2)).collect::<Vec<_>>();
        m.add_set("start", Set { nodes: end_nodes(0), ..Set::default() }).unwrap();
        m.add_set("end", Set { nodes: end_nodes(cells), ..Set::default() }).unwrap();
        let mut faces = Vec::new();
        for (l, e) in m.elements() {
            for (f, local) in e.kind.faces().iter().enumerate() {
                if local.iter().all(|n| m.node(e.nodes[*n]).unwrap().xyz[0] == LENGTH) {
                    faces.push(ElementFace { element: *l, face: f as u8 + 1 });
                }
            }
        }
        m.add_set("end_face", Set { faces, ..Set::default() }).unwrap();
        m.add_material(
            "Mat",
            Material { elastic: Some(Elastic::Isotropic { young: YOUNG, poisson: 0.0 }), density: Some(DENSITY), ..Material::default() },
        )
        .unwrap();
        let integration = if kind == ElementKind::Hex8 { Integration::Reduced } else { Integration::Full };
        let formulation = Formulation { integration, ..Formulation::default() };
        m.add_section("bar", Section { orientation: None, set: "all".into(), material: "Mat".into(), kind: SectionKind::Solid {}, formulation })
            .unwrap();
        m.add_step(Step::new("Dynamic", Procedure::ExplicitDynamic(step)), None).unwrap();
        Bar { model: m, kind, cells }
    }

    fn per_cell(&self) -> usize {
        if self.kind == ElementKind::Hex8 { 1 } else { 6 }
    }

    /// σ_xx per cell (the mean over its elements) and the cell's centre.
    fn axial_stress(&self, r: &ExplicitResult, frame: usize) -> Vec<(f64, f64)> {
        let f = &r.frames[frame];
        let h = LENGTH / self.cells as f64;
        (0..self.cells)
            .map(|i| {
                let s: f64 = f.stress[i * self.per_cell()..(i + 1) * self.per_cell()].iter().map(|s| s[0]).sum();
                ((i as f64 + 0.5) * h, s / self.per_cell() as f64)
            })
            .collect()
    }

    /// The axial velocity of the nodes on the bar's axis line y = z = 0, by x.
    fn axial_velocity(&self, r: &ExplicitResult, frame: usize) -> Vec<(f64, f64)> {
        let f = &r.frames[frame];
        let h = LENGTH / self.cells as f64;
        (0..=self.cells).map(|i| (i as f64 * h, f.velocity[r.nodes.binary_search(&(1 + i as u32)).unwrap()][0])).collect()
    }

    fn activate(&mut self, bc: Option<(&str, BoundaryCondition)>, load: Option<(&str, Load)>) {
        if let Some((name, bc)) = bc {
            self.model.add_boundary_condition(name, bc).unwrap();
            self.model.set_boundary_condition_active("Dynamic", name, true).unwrap();
        }
        if let Some((name, load)) = load {
            self.model.add_load(name, load).unwrap();
            self.model.set_load_active("Dynamic", name, true).unwrap();
        }
    }
}

/// Where a profile first falls below half of `level` going along x (linear interpolation).
fn front(profile: &[(f64, f64)], level: f64) -> f64 {
    for w in profile.windows(2) {
        let (a, b) = (w[0].1 / level, w[1].1 / level);
        if a >= 0.5 && b < 0.5 {
            return w[0].0 + (a - 0.5) / (a - b) * (w[1].0 - w[0].0);
        }
    }
    f64::NAN
}

fn mean_between(profile: &[(f64, f64)], from: f64, to: f64) -> f64 {
    let inside: Vec<f64> = profile.iter().filter(|(x, _)| *x >= from && *x <= to).map(|(_, v)| *v).collect();
    inside.iter().sum::<f64>() / inside.len() as f64
}

/// A pressure step on the free end of a bar sends a compressive wave in −x: front at L − c t,
/// stress −p, particle velocity −p/(ρc), so σ = ρ c v.
#[test]
fn pressure_wave_in_a_hex_bar() {
    let c = wave_speed();
    let t = 0.6 * LENGTH / c;
    let mut step = ExplicitDynamic::new(t);
    step.output_interval = Some(t / 4.0);
    let mut bar = Bar::new(ElementKind::Hex8, 200, step);
    let p = 1.0;
    bar.activate(None, Some(("push", Load::Pressure { set: "end_face".into(), magnitude: p })));
    let r = solve_explicit(&bar.model, None, &Limits::default()).unwrap();
    assert!(r.completed, "{:?}", r.warnings);
    let last = r.frames.len() - 1;
    assert!((r.frames[last].time - t).abs() < 1e-12 * t);
    // Measured from the loaded end, the profiles are those of a wave travelling in +x.
    let flip = |v: Vec<(f64, f64)>| -> Vec<(f64, f64)> { v.into_iter().rev().map(|(x, s)| (LENGTH - x, s)).collect() };
    let stress = flip(bar.axial_stress(&r, last));
    let velocity = flip(bar.axial_velocity(&r, last));
    let h = LENGTH / 200.0;
    let x_front = front(&stress, -p);
    let x_vfront = front(&velocity, -p / (DENSITY * c));
    let sigma = mean_between(&stress, 0.1 * c * t, 0.8 * c * t);
    let v = mean_between(&velocity, 0.1 * c * t, 0.8 * c * t);
    eprintln!(
        "hex8 pressure wave: front {x_front:.4} and {x_vfront:.4} vs c t = {:.4} (h = {h}); σ {sigma:.5} vs −p; v {v:.6} vs {:.6}; ρ c v = {:.5}; energy error {:.2e}; Δt {:.3e}, {} increments",
        c * t,
        -p / (DENSITY * c),
        DENSITY * c * v,
        r.energy.relative_error(r.initial_kinetic),
        r.initial_increment,
        r.increments
    );
    assert!((x_front - c * t).abs() < 2.0 * h, "stress front {x_front} vs {}", c * t);
    assert!((x_vfront - c * t).abs() < 2.0 * h, "velocity front {x_vfront} vs {}", c * t);
    assert!((sigma + p).abs() < 0.01 * p, "σ = {sigma}");
    assert!((DENSITY * c * v - sigma).abs() < 0.01 * p, "ρ c v = {} vs σ = {sigma}", DENSITY * c * v);
    // Work of the pressure = kinetic + strain energy (+ the viscous dissipation at the front).
    assert!(r.energy.external > 0.0 && r.energy.relative_error(r.initial_kinetic) < 0.01, "{:?}", r.energy);
    assert!(r.warnings.iter().all(|w| w.contains("unit system")), "{:?}", r.warnings);
    // The exact energy input is p A (c t)(p / E) = p² A c t / E; half is kinetic, half strain.
    let input = p * p * WIDTH * WIDTH * c * t / YOUNG;
    eprintln!("hex8 pressure wave: work {:.6e} vs p² A c t / E = {input:.6e}", r.energy.external);
    assert!((r.energy.external - input).abs() < 0.005 * input, "{} vs {input}", r.energy.external);
}

/// A velocity v prescribed at one end of a tetrahedral bar (nodal-pressure tet4): front at c t,
/// stress −ρ c v behind it.
#[test]
fn velocity_wave_in_a_tet_bar() {
    let c = wave_speed();
    let t = 0.6 * LENGTH / c;
    let v = 0.01;
    let mut step = ExplicitDynamic::new(t);
    step.output_interval = Some(t / 2.0);
    let mut bar = Bar::new(ElementKind::Tet4, 200, step);
    let push = BoundaryCondition::Displacement { target: Target::Set("start".into()), values: [Some(v * t), None, None, None, None, None] };
    bar.activate(Some(("push", push)), None);
    let r = solve_explicit(&bar.model, None, &Limits::default()).unwrap();
    let last = r.frames.len() - 1;
    let stress = bar.axial_stress(&r, last);
    let h = LENGTH / 200.0;
    let sigma_exact = -DENSITY * c * v;
    let x_front = front(&stress, sigma_exact);
    let sigma = mean_between(&stress, 0.1 * c * t, 0.8 * c * t);
    eprintln!(
        "tet4 velocity wave: front {x_front:.4} vs c t = {:.4} (h = {h}); σ {sigma:.6} vs −ρcv = {sigma_exact:.6}; energy error {:.2e}; Δt {:.3e}",
        c * t,
        r.energy.relative_error(r.initial_kinetic),
        r.initial_increment
    );
    assert!((x_front - c * t).abs() < 3.0 * h, "front {x_front} vs {}", c * t);
    assert!((sigma - sigma_exact).abs() < 0.02 * sigma_exact.abs(), "σ = {sigma}");
    assert!(r.energy.relative_error(r.initial_kinetic) < 0.01, "{:?}", r.energy);
    // The boundary condition does the work: ρ c v² A t (kinetic and strain energy).
    let input = DENSITY * c * v * v * WIDTH * WIDTH * t;
    eprintln!("tet4 velocity wave: work {:.6e} vs ρ c v² A t = {input:.6e}", r.energy.constraint);
    assert!((r.energy.constraint - input).abs() < 0.01 * input, "{} vs {input}", r.energy.constraint);
}

/// A bar moving at v₀ strikes a rigid wall: it stays in contact for 2L/c, then leaves at v₀ in
/// the opposite direction with its kinetic energy (elastic impact). The comparison is with the
/// elastic solution, so the linear bulk viscosity (which dissipates a few percent of the energy
/// at the wave fronts) is off; the kinematic wall stops the first layer of nodes dead, an O(h)
/// loss that the constraint work accounts for.
#[test]
fn bar_rebounds_from_a_rigid_wall() {
    let c = wave_speed();
    let v0 = 0.01;
    let end = 3.0 * LENGTH / c;
    for (kind, linear) in [(ElementKind::Hex8, 0.0), (ElementKind::Tet4, 0.0), (ElementKind::Hex8, 0.06)] {
        let mut step = ExplicitDynamic::new(end);
        step.history_interval = Some(end / 3000.0);
        step.output_interval = Some(end);
        step.initial_velocities = vec![InitialVelocity { set: "all_nodes".into(), velocity: [-v0, 0.0, 0.0] }];
        step.rigid_walls = vec![RigidWall { set: "start".into(), point: [0.0; 3], normal: [1.0, 0.0, 0.0], friction: Friction::NONE }];
        step.bulk_viscosity.linear = linear;
        let mut bar = Bar::new(kind, 200, ExplicitDynamic::new(end));
        let nodes: Vec<u32> = bar.model.nodes().keys().copied().collect();
        bar.model.add_set("all_nodes", Set { nodes, ..Set::default() }).unwrap();
        bar.model.replace_step(Step::new("Dynamic", Procedure::ExplicitDynamic(step))).unwrap();
        let r = solve_explicit(&bar.model, Some("Dynamic"), &Limits::default()).unwrap();
        let hist = &r.history;
        let peak = DENSITY * c * v0 * WIDTH * WIDTH;
        let force: Vec<f64> = hist.constraint_force.iter().map(|f| f[0]).collect();
        // Contact lasts until the wall force drops for good.
        let release = hist.time.iter().zip(&force).filter(|(_, f)| **f > 0.01 * peak).map(|(t, _)| *t).fold(0.0, f64::max);
        let rebound = hist.momentum.last().unwrap()[0] / r.total_mass;
        let error = r.energy.relative_error(r.initial_kinetic);
        let kinetic = r.energy.kinetic / r.initial_kinetic;
        let mid = hist.time.iter().position(|t| *t >= LENGTH / c).unwrap();
        eprintln!(
            "{} impact (linear q {linear}): contact {release:.5} vs 2L/c = {:.5}; wall force at L/c {:.4e} vs ρcv₀A = {peak:.4e}; rebound {rebound:.6} vs {v0}; final KE/KE₀ {kinetic:.4}; energy balance {error:.2e}; {:?}",
            kind.name(),
            2.0 * LENGTH / c,
            force[mid],
            r.energy
        );
        assert!(r.warnings.iter().all(|w| w.contains("unit system")), "{:?}", r.warnings);
        assert!(error < 0.01, "{} {:?}", kind.name(), r.energy);
        // The first layer of nodes (half an element's mass) stops dead on the wall.
        let stop = 0.5 * (0.5 * DENSITY * WIDTH * WIDTH * LENGTH / 200.0) * v0 * v0;
        assert!((r.energy.constraint + stop).abs() < 1e-6 * stop, "{} {:?}", kind.name(), r.energy);
        if linear > 0.0 {
            assert!(r.energy.viscous > 0.01 * r.initial_kinetic, "{:?}", r.energy);
            continue;
        }
        assert!((release - 2.0 * LENGTH / c).abs() < 0.01 * 2.0 * LENGTH / c, "{} contact {release}", kind.name());
        assert!((rebound - v0).abs() < 0.01 * v0, "{} rebound {rebound}", kind.name());
        assert!(kinetic > 0.99, "{} kinetic {kinetic}", kind.name());
        // While in contact the plateau force is ρ c v₀ A.
        assert!((force[mid] - peak).abs() < 0.03 * peak, "{} force {} vs {peak}", kind.name(), force[mid]);
    }
}

/// The same model gives the same numbers, bit for bit.
#[test]
fn runs_are_bit_identical() {
    let end = 0.3 * LENGTH / wave_speed();
    let mut step = ExplicitDynamic::new(end);
    step.hourglass = Hourglass { kind: HourglassKind::Stiffness, coefficient: 0.1 };
    let mut bar = Bar::new(ElementKind::Hex8, 40, step);
    bar.activate(
        None,
        Some(("push", Load::Force { target: Target::Set("end".into()), vector: [0.0, 0.02, -0.01], distribution: Distribution::Total })),
    );
    let a = solve_explicit(&bar.model, None, &Limits::default()).unwrap();
    let b = solve_explicit(&bar.model, None, &Limits::default()).unwrap();
    assert_eq!(a, b);
    assert!(a.increments > 10);
}

/// Hostile and unsupported models are errors, never panics; the caps hold.
#[test]
fn caps_and_unsupported_models_are_errors() {
    let end = 0.5 * LENGTH / wave_speed();
    let bar = Bar::new(ElementKind::Hex8, 20, ExplicitDynamic::new(end));
    let tight = |f: fn(&mut Limits)| {
        let mut l = Limits::default();
        f(&mut l);
        l
    };
    let e = solve_explicit(&bar.model, None, &tight(|l| l.max_increments = 3)).unwrap_err().to_string();
    assert!(e.contains("increments") && e.contains("element"), "{e}");
    assert!(solve_explicit(&bar.model, None, &tight(|l| l.max_dofs = 10)).is_err());
    assert!(solve_explicit(&bar.model, None, &tight(|l| l.max_output_values = 10)).is_err());
    assert!(solve_explicit(&bar.model, Some("Nope"), &Limits::default()).is_err());

    // A step capped below what it needs stops early, says so and keeps what it computed.
    let mut capped = ExplicitDynamic::new(end);
    capped.max_increments = 1_000_000;
    let mut short = Bar::new(ElementKind::Hex8, 20, capped);
    let r = solve_explicit(&short.model, None, &Limits::default()).unwrap();
    assert!(r.completed && r.increments > 5);
    let needed = r.increments;
    let Procedure::ExplicitDynamic(mut s) = short.model.step("Dynamic").unwrap().procedure.clone() else { unreachable!() };
    s.max_increments = needed;
    s.dt_scale = 0.89;
    short.model.replace_step(Step::new("Dynamic", Procedure::ExplicitDynamic(s))).unwrap();
    let r = solve_explicit(&short.model, None, &Limits::default());
    // Either refused up front (the estimate exceeds the cap) or stopped at the cap with a warning.
    match r {
        Ok(r) => assert!(!r.completed && r.warnings.iter().any(|w| w.contains("cap")) && r.time < end, "{r:?}"),
        Err(e) => assert!(e.to_string().contains("cap"), "{e}"),
    }

    // Fully integrated hex8, orthotropic tets and missing densities are refused by name.
    let mut full = Bar::new(ElementKind::Hex8, 4, ExplicitDynamic::new(end));
    full.model
        .replace_section(
            "bar",
            Section {
                orientation: None,
                set: "all".into(),
                material: "Mat".into(),
                kind: SectionKind::Solid {},
                formulation: Formulation::default(),
            },
        )
        .unwrap();
    assert!(solve_explicit(&full.model, None, &Limits::default()).unwrap_err().to_string().contains("not supported yet"));
    let mut ortho = Bar::new(ElementKind::Tet4, 4, ExplicitDynamic::new(end));
    let constants = Elastic::EngineeringConstants { young: [YOUNG; 3], poisson: [0.0; 3], shear: [YOUNG / 2.0; 3] };
    ortho.model.replace_material("Mat", Material { elastic: Some(constants), density: Some(DENSITY), ..Material::default() }).unwrap();
    assert!(solve_explicit(&ortho.model, None, &Limits::default()).unwrap_err().to_string().contains("orthotropic"));
    let mut light = Bar::new(ElementKind::Tet4, 4, ExplicitDynamic::new(end));
    light
        .model
        .replace_material("Mat", Material { elastic: Some(Elastic::Isotropic { young: YOUNG, poisson: 0.0 }), density: None, ..Material::default() })
        .unwrap();
    assert!(solve_explicit(&light.model, None, &Limits::default()).is_err());
    // A node behind its wall at the start.
    let mut behind = ExplicitDynamic::new(end);
    behind.rigid_walls = vec![RigidWall { set: "end".into(), point: [LENGTH + 1.0, 0.0, 0.0], normal: [-1.0, 0.0, 0.0], friction: Friction::NONE }];
    let mut wall = Bar::new(ElementKind::Hex8, 4, ExplicitDynamic::new(end));
    wall.model.replace_step(Step::new("Dynamic", Procedure::ExplicitDynamic(behind.clone()))).unwrap();
    assert!(solve_explicit(&wall.model, None, &Limits::default()).is_ok(), "the wall is beyond the end node, which is on its allowed side");
    behind.rigid_walls[0].point = [LENGTH - 1.0, 0.0, 0.0];
    wall.model.replace_step(Step::new("Dynamic", Procedure::ExplicitDynamic(behind))).unwrap();
    assert!(solve_explicit(&wall.model, None, &Limits::default()).unwrap_err().to_string().contains("behind"));
}

#[test]
fn explicit_refuses_a_named_material_orientation_until_it_rotates_materials() {
    let mut bar = Bar::new(ElementKind::Hex8, 1, ExplicitDynamic::new(0.001));
    let mut section = bar.model.section("bar").unwrap().clone();
    section.orientation = Some(SectionOrientation {
        frame: MaterialOrientation::Rectangular { axis1: [1.0, 1.0, 0.0], axis2: [0.0, 0.0, 1.0] },
        shell_reference_axis: 1,
    });
    bar.model.replace_section("bar", section).unwrap();
    let error = solve_explicit(&bar.model, None, &Limits::default()).unwrap_err().to_string();
    assert!(error.contains("section `bar`") && error.contains("material orientation") && error.contains("explicit dynamics"), "{error}");
}
