//! Large deformation and plasticity in explicit dynamics against closed-form solutions
//! (docs/solvers/literature.md, A4, A6, A10):
//! - simple shear of a hypoelastic material (J. K. Dienes, *On the analysis of rotation and
//!   stress rate in deforming bodies*, Acta Mech. 32:217–232, 1979): the Jaumann rate gives the
//!   oscillating σ12 = G sin γ, σ11 = −σ22 = G(1 − cos γ), the Green–Naghdi rate the monotonic
//!   solution in β = arctan(γ/2);
//! - a stressed element turned rigidly by 90° in many increments: the stress turns with it and
//!   keeps its magnitude (objectivity), where the small-strain update does not;
//! - uniaxial stress under large stretch with J2 plasticity, loaded and unloaded: Cauchy stress
//!   against the logarithmic strain follows the bilinear flow curve and unloads elastically;
//! - the Taylor bar (G. I. Taylor, Proc. R. Soc. A 194:289–299, 1948): the yield stress Taylor's
//!   one-dimensional estimate recovers from the computed final shape.
//!
//! Run with `--nocapture` to see the errors.

use fieldcraft_explicit::elements::{Motion, State, Workspace, internal_forces};
use fieldcraft_explicit::system::System;
use fieldcraft_explicit::{ExplicitResult, Failure, Limits, solve_explicit};
use fieldcraft_model::*;

/// A block [0, size] of `n` cells: one hex8 per cell, or six tet4 (Kuhn's split).
fn block(n: [usize; 3], size: [f64; 3], kind: ElementKind) -> Model {
    let mut m = Model::default();
    let label = |i: usize, j: usize, k: usize| (1 + i + (n[0] + 1) * (j + (n[1] + 1) * k)) as u32;
    for k in 0..=n[2] {
        for j in 0..=n[1] {
            for i in 0..=n[0] {
                let x = [i as f64 * size[0] / n[0] as f64, j as f64 * size[1] / n[1] as f64, k as f64 * size[2] / n[2] as f64];
                m.set_node(label(i, j, k), x).unwrap();
            }
        }
    }
    let corners = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]];
    let mut next = 1;
    for k in 0..n[2] {
        for j in 0..n[1] {
            for i in 0..n[0] {
                let at = |c: [usize; 3]| label(i + c[0], j + c[1], k + c[2]);
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
                    let p = v.map(|c| [0, 1, 2].map(|d| c[d] as f64 * size[d] / n[d] as f64));
                    let e = |a: usize| [0, 1, 2].map(|d| p[a][d] - p[0][d]);
                    let (a, b, c) = (e(1), e(2), e(3));
                    let vol = a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0]);
                    if vol < 0.0 {
                        v.swap(1, 2);
                    }
                    m.set_element(next, kind, v.iter().map(|c| at(*c)).collect()).unwrap();
                    next += 1;
                }
            }
        }
    }
    let all: Vec<u32> = m.elements().keys().copied().collect();
    m.add_set("all", Set { elements: all, ..Set::default() }).unwrap();
    m
}

fn node_set(m: &mut Model, name: &str, f: impl Fn([f64; 3]) -> bool) {
    let nodes: Vec<u32> = m.nodes().iter().filter(|(_, n)| f(n.xyz)).map(|(l, _)| *l).collect();
    m.add_set(name, Set { nodes, ..Set::default() }).unwrap();
}

fn define(m: &mut Model, kind: ElementKind, material: Material, step: ExplicitDynamic) {
    m.add_material("Mat", material).unwrap();
    let integration = if kind == ElementKind::Hex8 { Integration::Reduced } else { Integration::Full };
    let formulation = Formulation { integration, ..Formulation::default() };
    m.add_section("solid", Section { orientation: None, set: "all".into(), material: "Mat".into(), kind: SectionKind::Solid {}, formulation })
        .unwrap();
    m.add_step(Step::new("Dynamic", Procedure::ExplicitDynamic(step)), None).unwrap();
}

fn bc(m: &mut Model, name: &str, bc: BoundaryCondition) {
    m.add_boundary_condition(name, bc).unwrap();
    m.set_boundary_condition_active("Dynamic", name, true).unwrap();
}

fn elastic(young: f64, poisson: f64, density: f64) -> Material {
    Material { elastic: Some(Elastic::Isotropic { young, poisson }), density: Some(density), ..Material::default() }
}

/// Dienes' closed forms: (σ11, σ12) / G at shear γ.
fn dienes(rate: StressRate, gamma: f64) -> (f64, f64) {
    match rate {
        StressRate::Jaumann => (1.0 - gamma.cos(), gamma.sin()),
        StressRate::GreenNaghdi => {
            let beta = (gamma / 2.0).atan();
            let lc = beta.cos().ln();
            let s11 = 4.0 * ((2.0 * beta).cos() * lc + beta * (2.0 * beta).sin() - beta.sin().powi(2));
            let s12 = 2.0 * (2.0 * beta).cos() * (2.0 * beta - 2.0 * (2.0 * beta).tan() * lc - beta.tan());
            (s11, s12)
        }
    }
}

/// A unit cube sheared to γ = 6 by prescribing every node: x = X + γ(t) Y e_x.
#[test]
fn dienes_simple_shear_with_both_rates() {
    let (young, poisson) = (1.0, 0.3);
    let g = young / (2.0 * (1.0 + poisson));
    let gamma_end = 6.0;
    for kind in [ElementKind::Hex8, ElementKind::Tet4] {
        for rate in [StressRate::Jaumann, StressRate::GreenNaghdi] {
            let mut m = block([1, 1, 1], [1.0; 3], kind);
            let end = 4000.0;
            let mut step = ExplicitDynamic::new(end);
            step.stress_rate = rate;
            step.output_interval = Some(end / 60.0);
            step.max_increments = 100_000;
            define(&mut m, kind, elastic(young, poisson, 1.0), step);
            node_set(&mut m, "bottom", |x| x[1] == 0.0);
            node_set(&mut m, "top", |x| x[1] == 1.0);
            bc(&mut m, "hold", BoundaryCondition::Fixed { target: Target::Set("bottom".into()) });
            bc(
                &mut m,
                "shear",
                BoundaryCondition::Displacement {
                    target: Target::Set("top".into()),
                    values: [Some(gamma_end), Some(0.0), Some(0.0), None, None, None],
                },
            );
            let r = solve_explicit(&m, None, &Limits::default()).unwrap();
            assert!(r.completed);
            let mut worst: f64 = 0.0;
            for f in &r.frames {
                let gamma = gamma_end * f.time / end;
                let (s11, s12) = dienes(rate, gamma);
                for s in &f.stress {
                    let e = [(s[0] / g - s11).abs(), (s[1] / g + s11).abs(), (s[3] / g - s12).abs(), s[2].abs() / g, s[4].abs() / g, s[5].abs() / g];
                    worst = worst.max(e.into_iter().fold(0.0, f64::max));
                }
            }
            let last = r.frames.last().unwrap();
            eprintln!(
                "Dienes {} {rate:?}: σ/G at γ = 6: ({:.5}, {:.5}) vs ({:.5}, {:.5}); largest error over the run {worst:.2e} G; {} increments",
                kind.name(),
                last.stress[0][0] / g,
                last.stress[0][3] / g,
                dienes(rate, gamma_end).0,
                dienes(rate, gamma_end).1,
                r.increments
            );
            assert!(worst < 2e-3, "{} {rate:?}: {worst}", kind.name());
        }
    }
}

/// Drive one element through a sequence of nodal displacements with the solver's own element
/// update, returning its stress after each.
struct Driver {
    sys: System,
    state: State,
    ws: Workspace,
    u: Vec<f64>,
    time: f64,
}

impl Driver {
    fn new(model: &Model) -> Driver {
        let resolved = model.resolve().unwrap();
        let step = model.step("Dynamic").unwrap();
        let Procedure::ExplicitDynamic(settings) = &step.procedure else { panic!() };
        let sys = System::new(model, &resolved, step, settings, &Limits { max_dofs: 1000, ..Limits::default() }).unwrap();
        let (state, ws) = (State::new(&sys), Workspace::new(&sys));
        let u = vec![0.0; sys.dofs()];
        Driver { sys, state, ws, u, time: 0.0 }
    }

    /// Move the nodes to `x(X)` over an increment of length `dt`.
    fn go(&mut self, to: impl Fn([f64; 3]) -> [f64; 3], dt: f64) -> [f64; 6] {
        self.try_go(to, dt).unwrap()
    }

    fn try_go(&mut self, to: impl Fn([f64; 3]) -> [f64; 3], dt: f64) -> fieldcraft_explicit::Result<[f64; 6]> {
        let mut du = vec![0.0; self.u.len()];
        for (i, x0) in self.sys.xyz.iter().enumerate() {
            let x = to(*x0);
            for k in 0..3 {
                du[3 * i + k] = x[k] - x0[k] - self.u[3 * i + k];
                self.u[3 * i + k] = x[k] - x0[k];
            }
        }
        self.time += dt;
        let v: Vec<f64> = du.iter().map(|d| d / dt).collect();
        let bulk = BulkViscosity { linear: 0.0, quadratic: 0.0 };
        let m = Motion { u: &self.u, du: &du, v: &v, dt, time: self.time, bulk: &bulk };
        internal_forces(&self.sys, &m, &mut self.state, &mut self.ws)?;
        Ok(self.state.elements[0].point.stress)
    }
}

fn rotation_z(angle: f64, x: [f64; 3]) -> [f64; 3] {
    let (s, c) = angle.sin_cos();
    [c * x[0] - s * x[1], s * x[0] + c * x[1], x[2]]
}

/// Stretch an element 1% along x, then turn it rigidly by 90° about z in 50 increments: the
/// stress turns with it (σ_yy takes σ_xx's value), exactly for both rates, and the internal
/// forces stay those of the turned stress. The small-strain update keeps σ_xx and adds a
/// spurious stress from the rotation's "strain".
#[test]
fn rigid_rotation_turns_the_stress_with_the_body() {
    for (rate, nonlinear) in [(StressRate::Jaumann, true), (StressRate::GreenNaghdi, true), (StressRate::Jaumann, false)] {
        let mut m = block([1, 1, 1], [1.0; 3], ElementKind::Hex8);
        let mut step = ExplicitDynamic::new(1.0);
        step.stress_rate = rate;
        step.nonlinear_geometry = nonlinear;
        define(&mut m, ElementKind::Hex8, elastic(1000.0, 0.3, 1.0), step);
        let mut d = Driver::new(&m);
        let stretch = |x: [f64; 3]| [1.01 * x[0], x[1], x[2]];
        let before = d.go(stretch, 1.0);
        let n = 50;
        let mut after = before;
        for k in 1..=n {
            let angle = std::f64::consts::FRAC_PI_2 * k as f64 / n as f64;
            after = d.go(|x| rotation_z(angle, stretch(x)), 1.0);
        }
        // R σ Rᵀ for 90° about z: xx ↔ yy, xy → −xy, xz → yz, yz → −xz.
        let want = [before[1], before[0], before[2], -before[3], -before[5], before[4]];
        let error = after.iter().zip(want).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max) / before[0].abs();
        eprintln!("rigid rotation {rate:?} nonlinear {nonlinear}: σ before {before:.4?}, after {after:.4?}; error {error:.2e} of σ_xx");
        if nonlinear {
            assert!(error < 1e-12, "{rate:?}: {error}");
            // Forces of the turned stress balance on the turned element: ∑ f = 0, ∑ x × f = 0.
            let f = &d.ws.stress_force;
            let total: f64 = (0..3).map(|k| f.iter().skip(k).step_by(3).sum::<f64>().abs()).sum();
            assert!(total < 1e-10 * before[0].abs());
        } else {
            assert!(error > 0.5, "the small-strain update is not objective: {error}");
        }
    }
}

/// Uniaxial stress to a logarithmic strain of 0.3 and back by 0.006 (inside the elastic range) in a single hexahedron with
/// bilinear isotropic hardening: the lateral stretch is iterated each increment so the lateral
/// stresses vanish. The Cauchy stress follows σ = σ_y + H ε_p with ln λ = σ/E + ε_p (for a
/// hypoelastic–plastic law under coaxial stretching the integral of D is the logarithmic
/// strain), then unloads with slope E.
#[test]
fn j2_uniaxial_stress_strain_curve_and_unloading_at_large_strain() {
    let (young, poisson, y0, h) = (200_000.0, 0.3, 250.0, 2000.0);
    for kind in [ElementKind::Hex8, ElementKind::Tet4] {
        let mut m = block([1, 1, 1], [1.0; 3], kind);
        let plastic = Plastic::Bilinear { yield_stress: y0, hardening_modulus: h, hardening: Hardening::Isotropic };
        define(&mut m, kind, Material { plastic: Some(plastic), ..elastic(young, poisson, 7.8e-9) }, ExplicitDynamic::new(1.0));
        let mut d = Driver::new(&m);
        let mut path: Vec<f64> = (1..=600).map(|k| 0.3 * k as f64 / 600.0).collect();
        path.extend((1..=200).map(|k| 0.3 - 0.006 * k as f64 / 200.0));
        let (mut lateral, mut worst, mut peak) = (0.0f64, 0.0f64, (0.0, 0.0, 0.0));
        let mut last = [0.0; 6];
        for (step, axial) in path.iter().enumerate() {
            // Secant iteration on the lateral log strain for σ_yy = 0 (σ_zz follows by symmetry).
            let saved = (d.state.clone(), d.u.clone(), d.time);
            let trial = |lat: f64, d: &mut Driver| {
                (d.state, d.u, d.time) = saved.clone();
                let (a, l) = (axial.exp(), lat.exp());
                d.go(|x| [a * x[0], l * x[1], l * x[2]], 1e-6)
            };
            let (mut a, mut b) = (lateral, lateral - 1e-4);
            let (mut fa, mut fb) = (trial(a, &mut d)[1], trial(b, &mut d)[1]);
            for _ in 0..60 {
                if fb.abs() < 1e-9 || fb == fa {
                    break;
                }
                let c = b - fb * (b - a) / (fb - fa);
                (a, fa) = (b, fb);
                b = c;
                fb = trial(b, &mut d)[1];
            }
            lateral = b;
            let s = trial(b, &mut d);
            let p = d.state.elements.iter().map(|e| e.point.plastic_strain).fold(0.0, f64::max);
            if step < 600 {
                // Loading: σ = y0 + H ε_p and ε = σ/E + ε_p, so σ = y0 + E_t (ε − y0/E).
                let et = young * h / (young + h);
                let exact = if *axial * young <= y0 { young * axial } else { y0 + et * (axial - y0 / young) };
                worst = worst.max((s[0] - exact).abs() / exact);
                peak = (*axial, s[0], p);
            } else {
                // Unloading: elastic from the peak.
                let exact = peak.1 - young * (peak.0 - axial);
                worst = worst.max((s[0] - exact).abs() / peak.1);
                assert!((p - peak.2).abs() < 1e-12, "no plastic flow while unloading");
            }
            last = s;
        }
        eprintln!(
            "J2 uniaxial {}: peak σ {:.4} at ln λ = 0.3 (exact {:.4}); after unloading to 0.294 σ = {:.4}; largest relative error {worst:.2e}",
            kind.name(),
            peak.1,
            y0 + young * h / (young + h) * (0.3 - y0 / young),
            last[0]
        );
        assert!(worst < 1e-6, "{}: {worst}", kind.name());
    }
}

/// The Taylor bar: a copper-like bar of square section hits a rigid wall at 190 m/s (units mm,
/// t, s, MPa) and mushrooms. Taylor's one-dimensional analysis of a rigid–perfectly plastic bar
/// (Taylor 1948) relates the final length L₁ and the length X of the undeformed rear part to the
/// yield stress: σ_y = ρ U² (L − X) / (2 (L − L₁) ln(L / X)). Applied to the computed shape it
/// should give back the yield stress the material has, within the accuracy of the estimate
/// itself (which neglects elastic waves and radial inertia): 25%.
#[test]
fn taylor_bar_recovers_its_yield_stress() {
    let (rho, young, poisson, yield_stress) = (8.93e-9, 117_000.0, 0.35, 400.0);
    let (length, width, speed) = (25.4, 7.6, 190_000.0);
    let end = 80e-6;
    let mut step = ExplicitDynamic::new(end);
    step.output_interval = Some(end / 4.0);
    step.history_interval = Some(end / 400.0);
    step.initial_velocities = vec![InitialVelocity { set: "nodes".into(), velocity: [0.0, 0.0, -speed] }];
    step.rigid_walls = vec![RigidWall { set: "nodes".into(), point: [0.0; 3], normal: [0.0, 0.0, 1.0], friction: Friction::NONE }];
    step.hourglass = Hourglass { kind: HourglassKind::Stiffness, coefficient: 0.1 };
    let mut m = block([4, 4, 30], [width, width, length], ElementKind::Hex8);
    node_set(&mut m, "nodes", |_| true);
    let plastic = Plastic::Bilinear { yield_stress, hardening_modulus: 0.0, hardening: Hardening::Isotropic };
    define(&mut m, ElementKind::Hex8, Material { plastic: Some(plastic), ..elastic(young, poisson, rho) }, step);
    let r = solve_explicit(&m, None, &Limits::default()).unwrap();
    let f = r.frames.last().unwrap();
    // Final length: the top face's mean height; undeformed part: the axis above the last
    // element with plastic strain, along the centre column of elements.
    let z = |label: u32| m.node(label).unwrap().xyz[2] + f.displacement[r.nodes.binary_search(&label).unwrap()][2];
    let top: Vec<u32> = m.nodes().iter().filter(|(_, n)| n.xyz[2] == length).map(|(l, _)| *l).collect();
    let l1 = top.iter().map(|l| z(*l)).sum::<f64>() / top.len() as f64;
    let column: Vec<usize> = (0..30).map(|k| 4 * 4 * k + 4 + 1).collect();
    // The rear part counts as undeformed from the first element (from the impact end) whose
    // plastic strain is below 1%, small against the mushroom's.
    let deformed = column.iter().position(|&i| f.plastic_strain[i] < 0.01).unwrap_or(30);
    let x = length * (30 - deformed) as f64 / 30.0;
    let estimate = rho * speed * speed * (length - x) / (2.0 * (length - l1) * (length / x).ln());
    let kinetic = r.energy.kinetic / r.initial_kinetic;
    eprintln!(
        "Taylor bar: L₁/L = {:.4}, X/L = {:.4}; Taylor's estimate σ_y = {estimate:.1} vs {yield_stress} ({:+.1}%); final KE/KE₀ {kinetic:.4}; plastic work {:.3e} of {:.3e}; balance {:.2e}; {} increments",
        l1 / length,
        x / length,
        100.0 * (estimate / yield_stress - 1.0),
        r.energy.plastic,
        r.initial_kinetic,
        r.energy.relative_error(r.initial_kinetic),
        r.increments
    );
    assert!(r.completed && l1 < 0.95 * length && l1 > 0.5 * length);
    assert!((estimate / yield_stress - 1.0).abs() < 0.25, "{estimate}");
    assert!(r.energy.relative_error(r.initial_kinetic) < 0.01, "{:?}", r.energy);
    assert!(r.energy.plastic > 0.5 * r.initial_kinetic, "most of the impact energy goes into plastic work: {:?}", r.energy);
}

fn taylor_with(material: Material) -> (Model, ExplicitResult) {
    let end = 40e-6;
    let mut step = ExplicitDynamic::new(end);
    step.output_interval = Some(end / 2.0);
    step.initial_velocities = vec![InitialVelocity { set: "nodes".into(), velocity: [0.0, 0.0, -300_000.0] }];
    step.rigid_walls = vec![RigidWall { set: "nodes".into(), point: [0.0; 3], normal: [0.0, 0.0, 1.0], friction: Friction::NONE }];
    step.hourglass = Hourglass { kind: HourglassKind::Stiffness, coefficient: 0.1 };
    let mut m = block([2, 2, 12], [7.6, 7.6, 25.4], ElementKind::Hex8);
    node_set(&mut m, "nodes", |_| true);
    define(&mut m, ElementKind::Hex8, material, step);
    let r = solve_explicit(&m, None, &Limits::default()).unwrap();
    (m, r)
}

/// Johnson–Cook strength, adiabatic heating, damage and erosion in a fast impact: elements at
/// the impact face reach the erosion strain (compression keeps their damage low, as the
/// fracture strain grows with confinement), are removed and reported; the run stays balanced
/// and repeatable bit for bit.
#[test]
fn johnson_cook_impact_erodes_elements() {
    // OFHC copper constants as published by Johnson & Cook (1983, 1985), units mm, t, s, MPa, K.
    let plastic = Plastic::JohnsonCook {
        a: 90.0,
        b: 292.0,
        n: 0.31,
        c: 0.025,
        m: 1.09,
        reference_strain_rate: 1.0,
        room_temperature: 294.0,
        melt_temperature: 1356.0,
        specific_heat: Some(3.83e8),
        taylor_quinney: 0.9,
    };
    let damage = Damage::JohnsonCook { d1: 0.54, d2: 4.89, d3: -3.03, d4: 0.014, d5: 1.12 };
    let material = Material {
        plastic: Some(plastic),
        damage: Some(damage),
        erosion: Some(Erosion { plastic_strain: 0.8 }),
        ..elastic(124_000.0, 0.34, 8.96e-9)
    };
    let (_, r) = taylor_with(material.clone());
    let f = r.frames.last().unwrap();
    let hottest = f.temperature.iter().copied().fold(0.0, f64::max);
    eprintln!(
        "Johnson–Cook impact: {} elements eroded ({:?}…), eroded energy {:.3e}; hottest element {hottest:.1} K; largest damage {:.3}; balance {:.2e}",
        r.eroded.len(),
        r.eroded.first(),
        r.energy.eroded,
        f.damage.iter().copied().fold(0.0, f64::max),
        r.energy.relative_error(r.initial_kinetic)
    );
    assert!(!r.eroded.is_empty() && r.eroded.iter().all(|e| matches!(e.failure, Failure::Damage | Failure::PlasticStrain)));
    assert_eq!(f.failed.iter().filter(|x| **x).count(), r.eroded.len());
    assert!(hottest > 294.0 + 10.0);
    assert!(r.energy.eroded > 0.0 && r.warnings.iter().any(|w| w.contains("eroded")));
    assert!(r.energy.relative_error(r.initial_kinetic) < 0.02, "{:?}", r.energy);
    let (_, again) = taylor_with(material);
    assert!(r == again, "bit-identical runs");
}

/// Damage with a constant fracture strain D1 (D2 = D4 = D5 = 0): an element stretched in
/// tension fails when its plastic strain reaches D1, and from then on carries no stress.
#[test]
fn damage_fails_an_element_at_its_fracture_strain() {
    let mut m = block([1, 1, 1], [1.0; 3], ElementKind::Hex8);
    let plastic = Plastic::Bilinear { yield_stress: 250.0, hardening_modulus: 1000.0, hardening: Hardening::Isotropic };
    let damage = Damage::JohnsonCook { d1: 0.2, d2: 0.0, d3: 0.0, d4: 0.0, d5: 0.0 };
    define(
        &mut m,
        ElementKind::Hex8,
        Material { plastic: Some(plastic), damage: Some(damage), ..elastic(200_000.0, 0.3, 7.8e-9) },
        ExplicitDynamic::new(1.0),
    );
    let mut d = Driver::new(&m);
    let mut failed_at = None;
    for k in 1..=1000 {
        let stretch = 1.0 + 0.5 * k as f64 / 1000.0;
        let before = d.state.elements[0].point.plastic_strain;
        d.go(|x| [stretch * x[0], x[1] / stretch.sqrt(), x[2] / stretch.sqrt()], 1e-6);
        let e = &d.state.elements[0];
        if e.failed.is_some() {
            failed_at = Some((before, e.point.plastic_strain, e.failed));
            break;
        }
    }
    let (before, at, why) = failed_at.expect("the element fails");
    eprintln!("damage: failed between plastic strain {before:.5} and {at:.5} (D1 = 0.2)");
    assert_eq!(why, Some(Failure::Damage));
    assert!(before < 0.2 && at >= 0.2 - 1e-12, "{before} {at}");
    // Failed: no forces from it.
    d.go(|x| [1.6 * x[0], x[1], x[2]], 1e-6);
    assert!(d.ws.stress_force.iter().all(|f| *f == 0.0));
}

/// Inelastic models that cannot run are refused by name, never panics.
#[test]
fn inelastic_hostile_models_are_errors() {
    // Plasticity with orthotropic elasticity is refused by the model itself.
    let mut m = block([1, 1, 1], [1.0; 3], ElementKind::Hex8);
    let ortho = Elastic::EngineeringConstants { young: [1.0; 3], poisson: [0.0; 3], shear: [0.5; 3] };
    let plastic = Plastic::Bilinear { yield_stress: 1.0, hardening_modulus: 0.0, hardening: Hardening::Isotropic };
    assert!(m.add_material("Bad", Material { elastic: Some(ortho), density: Some(1.0), plastic: Some(plastic), ..Material::default() }).is_err());
    // An element turned inside out stops the run with an error naming it.
    define(&mut m, ElementKind::Hex8, elastic(1.0, 0.3, 1.0), ExplicitDynamic::new(1.0));
    let mut d = Driver::new(&m);
    d.go(|x| [x[0], x[1], 0.5 * x[2]], 1.0);
    let e = d.try_go(|x| [x[0], x[1], 0.5 - x[2]], 1.0).unwrap_err().to_string();
    assert!(e.contains("element 1") && e.contains("inside out"), "{e}");
    // A crushed element shrinks the increment until the cap stops the run, with a warning.
    let mut m = block([1, 1, 1], [1.0; 3], ElementKind::Hex8);
    let mut step = ExplicitDynamic::new(1.0);
    step.max_increments = 20_000;
    define(&mut m, ElementKind::Hex8, elastic(1.0, 0.3, 1.0), step);
    node_set(&mut m, "top", |x| x[2] == 1.0);
    node_set(&mut m, "bottom", |x| x[2] == 0.0);
    bc(&mut m, "hold", BoundaryCondition::Fixed { target: Target::Set("bottom".into()) });
    bc(
        &mut m,
        "crush",
        BoundaryCondition::Displacement { target: Target::Set("top".into()), values: [Some(0.0), Some(0.0), Some(-2.0), None, None, None] },
    );
    match solve_explicit(&m, None, &Limits::default()) {
        Ok(r) => assert!(!r.completed && r.warnings.iter().any(|w| w.contains("cap"))),
        Err(e) => assert!(e.to_string().contains("element 1"), "{e}"),
    }
}
