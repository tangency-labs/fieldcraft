//! Analytic verification of the structural elements: shell patch tests (membrane and bending,
//! quad4 and tri3, in a tilted plane), the simply supported plate under uniform pressure against
//! Navier's series (deflection and first frequency), cantilever beams against Timoshenko and
//! Euler–Bernoulli theory (deflection and frequencies), a portal frame against the
//! slope-deflection hand calculation, and mixed meshes of solids with shells and rods.
//!
//! Every expected value is a closed form or a converged series; nothing comes from another
//! program.

use std::collections::BTreeMap;

use fieldcraft_model::{
    BeamProfile, BoundaryCondition, Distribution, Elastic, ElementFace, ElementKind, Formulation, Load, Material, Model, Procedure, Section,
    SectionKind, Set, Step, Target,
};
use fieldcraft_solver::{FrequencyOptions, Limits, MassKind, StaticResult, solve_frequency, solve_static};

/// A model under construction: nodes and elements numbered from 1, one material, one step.
struct Build {
    model: Model,
    xyz: BTreeMap<u32, [f64; 3]>,
    elements: u32,
}

impl Build {
    fn new(elastic: Elastic, density: f64, procedure: Procedure) -> Build {
        let mut model = Model::default();
        model.add_material("Mat", Material { elastic: Some(elastic), density: Some(density), ..Material::default() }).unwrap();
        model.add_step(Step::new("Step", procedure), None).unwrap();
        Build { model, xyz: BTreeMap::new(), elements: 0 }
    }

    fn node(&mut self, label: u32, xyz: [f64; 3]) {
        self.model.set_node(label, xyz).unwrap();
        self.xyz.insert(label, xyz);
    }

    fn element(&mut self, kind: ElementKind, nodes: Vec<u32>) -> u32 {
        self.elements += 1;
        self.model.set_element(self.elements, kind, nodes).unwrap();
        self.elements
    }

    fn section(&mut self, name: &str, elements: Vec<u32>, kind: SectionKind) {
        self.model.add_set(name, Set { elements, ..Set::default() }).unwrap();
        self.model
            .add_section(name, Section { orientation: None, set: name.into(), material: "Mat".into(), kind, formulation: Formulation::default() })
            .unwrap();
    }

    fn nodes_where(&self, f: impl Fn([f64; 3]) -> bool) -> Vec<u32> {
        self.xyz.iter().filter(|(_, x)| f(**x)).map(|(l, _)| *l).collect()
    }

    /// A boundary condition on the nodes where `f` holds: `values` per direction (None = free).
    fn hold(&mut self, name: &str, values: [Option<f64>; 6], f: impl Fn([f64; 3]) -> bool) {
        let nodes = self.nodes_where(f);
        self.hold_nodes(name, values, nodes);
    }

    fn hold_nodes(&mut self, name: &str, values: [Option<f64>; 6], nodes: Vec<u32>) {
        self.model.add_set(name, Set { nodes, ..Set::default() }).unwrap();
        self.model.add_boundary_condition(name, BoundaryCondition::Displacement { target: Target::Set(name.into()), values }).unwrap();
        self.model.set_boundary_condition_active("Step", name, true).unwrap();
    }

    fn load(&mut self, name: &str, load: Load) {
        self.model.add_load(name, load).unwrap();
        self.model.set_load_active("Step", name, true).unwrap();
    }

    fn force(&mut self, name: &str, nodes: Vec<u32>, vector: [f64; 3]) {
        self.model.add_set(name, Set { nodes, ..Set::default() }).unwrap();
        self.load(name, Load::Force { target: Target::Set(name.into()), vector, distribution: Distribution::Total });
    }

    /// A structured surface mesh of n[0] × n[1] cells on the map of [0, 1]², quadrilaterals or
    /// two triangles per cell, starting at node label `base`; returns the elements.
    fn surface(&mut self, n: [usize; 2], kind: ElementKind, base: u32, map: impl Fn(f64, f64) -> [f64; 3]) -> Vec<u32> {
        let label = |i: usize, j: usize| base + (i + (n[0] + 1) * j) as u32;
        for j in 0..=n[1] {
            for i in 0..=n[0] {
                self.node(label(i, j), map(i as f64 / n[0] as f64, j as f64 / n[1] as f64));
            }
        }
        let mut out = Vec::new();
        for j in 0..n[1] {
            for i in 0..n[0] {
                let (a, b, c, d) = (label(i, j), label(i + 1, j), label(i + 1, j + 1), label(i, j + 1));
                if kind == ElementKind::Quad4 {
                    out.push(self.element(kind, vec![a, b, c, d]));
                } else {
                    out.push(self.element(kind, vec![a, b, c]));
                    out.push(self.element(kind, vec![a, c, d]));
                }
            }
        }
        out
    }

    fn solve(&self) -> StaticResult {
        solve_static(&self.model, None, &Limits::default()).unwrap()
    }
}

fn at(r: &StaticResult, label: u32) -> usize {
    r.nodes.binary_search(&label).unwrap()
}

fn rel(a: f64, b: f64) -> f64 {
    (a - b).abs() / b.abs()
}

const FREE: Option<f64> = None;

// ---------------------------------------------------------------------------------------------
// Shell patch tests
// ---------------------------------------------------------------------------------------------

/// The distorted patch of five quadrilaterals in a 0.24 × 0.12 rectangle (the irregular patch of
/// R. H. MacNeal, R. L. Harder, *A proposed standard set of problems to test finite element
/// accuracy*, Finite Elements in Analysis and Design 1(1):3–20, 1985), in local plane
/// coordinates; the four inner nodes are 5–8.
const PATCH: [[f64; 2]; 8] = [[0.0, 0.0], [0.24, 0.0], [0.24, 0.12], [0.0, 0.12], [0.04, 0.02], [0.18, 0.03], [0.16, 0.08], [0.08, 0.08]];
/// Thick enough that membrane and bending stiffness are within a few orders of magnitude, so the
/// comparison measures the elements and not the round-off of an ill-conditioned system.
const PATCH_THICKNESS: f64 = 0.02;
const PATCH_QUADS: [[u32; 4]; 5] = [[1, 2, 6, 5], [2, 3, 7, 6], [3, 4, 8, 7], [4, 1, 5, 8], [5, 6, 7, 8]];

/// A tilted plane: local (x, y) → global, with its axes.
fn tilted(p: [f64; 2]) -> [f64; 3] {
    let (e1, e2) = axes();
    [0, 1, 2].map(|k| 1.0 + p[0] * e1[k] + p[1] * e2[k])
}

fn axes() -> ([f64; 3], [f64; 3]) {
    let e1 = [0.6, 0.8, 0.0];
    let (c, s) = (0.5f64.sqrt(), 0.5f64.sqrt());
    // e2 = cos · (−0.8, 0.6, 0) + sin · z.
    let e2 = [-0.8 * c, 0.6 * c, s];
    (e1, e2)
}

fn normal() -> [f64; 3] {
    let (a, b) = axes();
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

/// Local displacement (u, v, w) and rotations (θx, θy, θz) to global 6-vector.
fn to_global(local: [f64; 6]) -> [f64; 6] {
    let (e1, e2) = axes();
    let e3 = normal();
    let g = |a: [f64; 3]| [0, 1, 2].map(|k| a[0] * e1[k] + a[1] * e2[k] + a[2] * e3[k]);
    let t = g([local[0], local[1], local[2]]);
    let r = g([local[3], local[4], local[5]]);
    [t[0], t[1], t[2], r[0], r[1], r[2]]
}

/// The patch with a prescribed exact field on its boundary nodes 1–4: returns the worst nodal
/// error of the inner nodes (relative to the field's size) and the solution.
fn patch(kind: ElementKind, field: impl Fn([f64; 2]) -> [f64; 6]) -> (f64, StaticResult, Build) {
    let mut b = Build::new(Elastic::Isotropic { young: 1e6, poisson: 0.25 }, 1.0, Procedure::Static {});
    for (i, p) in PATCH.iter().enumerate() {
        b.node(i as u32 + 1, tilted(*p));
    }
    let mut elements = Vec::new();
    for q in PATCH_QUADS {
        if kind == ElementKind::Quad4 {
            elements.push(b.element(kind, q.to_vec()));
        } else {
            elements.push(b.element(kind, vec![q[0], q[1], q[2]]));
            elements.push(b.element(kind, vec![q[0], q[2], q[3]]));
        }
    }
    b.section("shell", elements, SectionKind::Shell { thickness: PATCH_THICKNESS, offset: 0.0, points: 5 });
    for (i, p) in PATCH.iter().enumerate().take(4) {
        let g = to_global(field(*p));
        b.hold_nodes(&format!("edge{i}"), g.map(Some), vec![i as u32 + 1]);
    }
    let r = b.solve();
    let size = PATCH.iter().flat_map(|p| field(*p)).fold(0.0f64, |m, v| m.max(v.abs()));
    let mut worst = 0.0f64;
    for (i, p) in PATCH.iter().enumerate().skip(4) {
        let want = to_global(field(*p));
        let k = at(&r, i as u32 + 1);
        let got = [r.displacement[k], r.rotation[k]].concat();
        for c in 0..6 {
            worst = worst.max((got[c] - want[c]).abs() / size);
        }
    }
    (worst, r, b)
}

#[test]
fn shell_membrane_patch_test() {
    // u = 10⁻³ (x + y/2), v = 10⁻³ (y + x/2): constant strain, no in-plane rotation.
    let field = |p: [f64; 2]| [1e-3 * (p[0] + p[1] / 2.0), 1e-3 * (p[1] + p[0] / 2.0), 0.0, 0.0, 0.0, 0.0];
    for kind in [ElementKind::Quad4, ElementKind::Tri3] {
        let (err, r, _) = patch(kind, field);
        // σ in local axes: E/(1 − ν²) [[1, ν], [ν, 1]] ε, τ = G γ.
        let (e, nu) = (1e6, 0.25);
        let c = e / (1.0 - nu * nu);
        let (sx, sy, txy) = (c * (1e-3 + nu * 1e-3), c * (1e-3 + nu * 1e-3), e / (2.0 * (1.0 + nu)) * 1e-3);
        let mut stress_err = 0.0f64;
        for s in &r.shells {
            let [a1, a2, _] = s.stress.axes;
            for layer in [&s.stress.bottom, &s.stress.middle, &s.stress.top] {
                for g in layer {
                    // Back to the patch's local axes: σ_ab = a·S·b.
                    let t = [[g[0], g[3], g[4]], [g[3], g[1], g[5]], [g[4], g[5], g[2]]];
                    let (e1, e2) = axes();
                    let p = |a: [f64; 3], b: [f64; 3]| (0..3).map(|i| (0..3).map(|j| a[i] * t[i][j] * b[j]).sum::<f64>()).sum::<f64>();
                    stress_err = stress_err.max(((p(e1, e1) - sx).abs()).max((p(e2, e2) - sy).abs()).max((p(e1, e2) - txy).abs()) / sx);
                }
            }
            assert!(a1.iter().chain(&a2).all(|v| v.is_finite()));
        }
        eprintln!("membrane patch {kind:?}: displacement error {err:.1e}, stress error {stress_err:.1e}");
        assert!(err < 1e-9 && stress_err < 1e-9, "{kind:?}: {err:e} {stress_err:e}");
    }
}

#[test]
fn shell_bending_patch_test() {
    // w = 10⁻³ (x² + x y + y²)/2 with θx = w,y and θy = −w,x: constant curvature, no shear.
    let k = 1e-3;
    let field = |p: [f64; 2]| {
        let [x, y] = p;
        [0.0, 0.0, k * (x * x + x * y + y * y) / 2.0, k * (x / 2.0 + y), -k * (x + y / 2.0), 0.0]
    };
    for kind in [ElementKind::Quad4, ElementKind::Tri3] {
        let (err, r, _) = patch(kind, field);
        // κx = θy,x = −k, κy = −θx,y = −k, κxy = θy,y − θx,x = −k.
        let (e, nu, t) = (1e6, 0.25, PATCH_THICKNESS);
        let d = e * t * t * t / (12.0 * (1.0 - nu * nu));
        let want = [d * (-k - nu * k), d * (-k - nu * k), d * (1.0 - nu) / 2.0 * -k];
        let mut moment_err = 0.0f64;
        for s in &r.shells {
            // Local element axes may differ from the patch's in-plane axes by a rotation in the
            // plane; the invariants Mx + My and Mx My − Mxy² do not.
            for f in &s.stress.forces {
                let (tr, det) = (f[3] + f[4], f[3] * f[4] - f[5] * f[5]);
                let (wtr, wdet) = (want[0] + want[1], want[0] * want[1] - want[2] * want[2]);
                moment_err = moment_err.max(rel(tr, wtr)).max(rel(det, wdet)).max(f[6].abs().max(f[7].abs()) / wtr.abs());
            }
        }
        eprintln!("bending patch {kind:?}: displacement error {err:.1e}, moment error {moment_err:.1e}");
        assert!(err < 1e-9 && moment_err < 1e-8, "{kind:?}: {err:e} {moment_err:e}");
    }
}

// ---------------------------------------------------------------------------------------------
// Simply supported plate (Navier)
// ---------------------------------------------------------------------------------------------

/// Navier's series for the centre deflection of a simply supported square plate under uniform
/// pressure: w = 16 q a⁴/(π⁶ D) Σ_{m, n odd} (−1)^((m+n)/2 − 1) / (m n (m² + n²)²)
/// (S. Timoshenko, S. Woinowsky-Krieger, *Theory of Plates and Shells*, 2nd ed., McGraw-Hill
/// 1959, §30).
fn navier(q: f64, a: f64, d: f64) -> f64 {
    let mut s = 0.0;
    for m in (1..400).step_by(2) {
        for n in (1..400).step_by(2) {
            let (mf, nf) = (f64::from(m), f64::from(n));
            let sign = if ((m + n) / 2 - 1) % 2 == 0 { 1.0 } else { -1.0 };
            s += sign / (mf * nf * (mf * mf + nf * nf).powi(2));
        }
    }
    16.0 * q * a.powi(4) / (std::f64::consts::PI.powi(6) * d) * s
}

/// A quarter of a simply supported square plate (side a, edges x = 0 and y = 0 held, symmetry
/// on x = y = a/2), n × n cells, under pressure q or for its modes.
fn plate(kind: ElementKind, n: usize, procedure: Procedure) -> Build {
    let (a, t, e, nu, rho) = (1.0, 0.01, 1e7, 0.3, 1.0);
    let mut b = Build::new(Elastic::Isotropic { young: e, poisson: nu }, rho, procedure);
    let elements = b.surface([n, n], kind, 1, |u, v| [u * a / 2.0, v * a / 2.0, 0.0]);
    b.section("plate", elements.clone(), SectionKind::Shell { thickness: t, offset: 0.0, points: 5 });
    let on = |x: f64, c: f64| (x - c).abs() < 1e-9;
    // Hard simple support: w and the twist along the edge held.
    b.hold("x0", [FREE, FREE, Some(0.0), Some(0.0), FREE, FREE], |p| on(p[0], 0.0));
    b.hold("y0", [FREE, FREE, Some(0.0), FREE, Some(0.0), FREE], |p| on(p[1], 0.0));
    b.hold("xs", [Some(0.0), FREE, FREE, FREE, Some(0.0), Some(0.0)], |p| on(p[0], a / 2.0));
    b.hold("ys", [FREE, Some(0.0), FREE, Some(0.0), FREE, Some(0.0)], |p| on(p[1], a / 2.0));
    if matches!(b.model.steps()[0].procedure, Procedure::Static {}) {
        let faces = elements.iter().map(|e| ElementFace { element: *e, face: 1 }).collect();
        b.model.add_set("top", Set { faces, ..Set::default() }).unwrap();
        b.load("q", Load::Pressure { set: "top".into(), magnitude: 1.0 });
    }
    b
}

#[test]
fn simply_supported_plate_deflection_and_frequency() {
    let (a, t, e, nu, rho) = (1.0, 0.01, 1e7, 0.3, 1.0);
    let d = e * t * t * t / (12.0 * (1.0 - nu * nu));
    let w = navier(1.0, a, d);
    let f11 = std::f64::consts::PI * 2.0 / (a * a) * (d / (rho * t)).sqrt() / 2.0;
    for kind in [ElementKind::Quad4, ElementKind::Tri3] {
        let mut errors = Vec::new();
        for n in [4, 8, 16] {
            let b = plate(kind, n, Procedure::Static {});
            let r = b.solve();
            let centre = b.nodes_where(|p| (p[0] - 0.5).abs() < 1e-9 && (p[1] - 0.5).abs() < 1e-9)[0];
            let got = -r.displacement[at(&r, centre)][2];
            assert!((r.reaction_total[2] + r.applied_total[2]).abs() < 1e-9 * r.applied_total[2].abs());
            let b = plate(kind, n, Procedure::Frequency { modes: 1, min_hz: None, max_hz: None });
            let fr = solve_frequency(&b.model, None, &FrequencyOptions::default(), &Limits::default()).unwrap();
            errors.push((n, rel(got, w), rel(fr.modes[0].frequency, f11)));
        }
        eprintln!("plate {kind:?} vs Navier w = {w:.6e}, f11 = {f11:.6}: (cells, deflection error, frequency error) {errors:.5?}");
        let (_, dw, df) = errors[2];
        assert!(dw < 0.01 && df < 0.01, "{kind:?} {errors:?}");
        assert!(errors[0].1 > errors[2].1, "converges: {errors:?}");
    }
}

// ---------------------------------------------------------------------------------------------
// Beams
// ---------------------------------------------------------------------------------------------

/// A straight cantilever along a skew direction, `n` beam elements, clamped at node 1.
fn cantilever(n: usize, length: f64, profile: BeamProfile, procedure: Procedure) -> (Build, [f64; 3], [f64; 3]) {
    let dir = [2.0 / 3.0, 2.0 / 3.0, 1.0 / 3.0];
    let orientation = [0.0, 0.0, 1.0];
    let mut b = Build::new(Elastic::Isotropic { young: 210e3, poisson: 0.3 }, 7.85e-9, procedure);
    for i in 0..=n {
        let s = length * i as f64 / n as f64;
        b.node(i as u32 + 1, dir.map(|c| c * s));
    }
    let elements: Vec<u32> = (1..=n as u32).map(|i| b.element(ElementKind::Beam2, vec![i, i + 1])).collect();
    b.section("beam", elements, SectionKind::Beam { profile, orientation });
    b.hold_nodes("clamp", [Some(0.0); 6], vec![1]);
    // Local y: the orientation made normal to the axis; local z = x × y.
    let dot = dir[2];
    let y = [0, 1, 2].map(|k| orientation[k] - dot * dir[k]);
    let ny = (y[0] * y[0] + y[1] * y[1] + y[2] * y[2]).sqrt();
    let y = y.map(|c| c / ny);
    let z = [dir[1] * y[2] - dir[2] * y[1], dir[2] * y[0] - dir[0] * y[2], dir[0] * y[1] - dir[1] * y[0]];
    (b, y, z)
}

#[test]
fn cantilever_beam_deflection_is_exact() {
    // A stocky I beam, so shear deformation matters: δ = P L³/(3 E I) + P L/(G A_s) in each plane.
    let profile = BeamProfile::I { height: 200.0, width: 100.0, web: 6.0, flange: 10.0 };
    let length: f64 = 1000.0;
    let p = fieldcraft_fem::beam::properties(&profile, 0.3);
    let (e, g) = (210e3, 210e3 / 2.6);
    for n in [1, 4] {
        for (along_y, name) in [(true, "weak"), (false, "strong")] {
            let (mut b, y, z) = cantilever(n, length, profile, Procedure::Static {});
            let d = if along_y { y } else { z };
            b.force("tip", vec![n as u32 + 1], d.map(|c| 1000.0 * c));
            let r = b.solve();
            let u = r.displacement[at(&r, n as u32 + 1)];
            let got = u[0] * d[0] + u[1] * d[1] + u[2] * d[2];
            let (i, a_s) = if along_y { (p.iz, p.shear_y) } else { (p.iy, p.shear_z) };
            let want = 1000.0 * length.powi(3) / (3.0 * e * i) + 1000.0 * length / (g * a_s);
            eprintln!(
                "cantilever {n} element(s), {name} axis: tip {got:.6e} vs {want:.6e} (shear share {:.1}%)",
                100.0 * 1000.0 * length / (g * a_s) / want
            );
            assert!(rel(got, want) < 1e-9, "{got} {want}");
            // The clamp's moment and the root section force.
            let m = r.reaction_moment[at(&r, 1)];
            let mm = (m[0] * m[0] + m[1] * m[1] + m[2] * m[2]).sqrt();
            assert!(rel(mm, 1000.0 * length) < 1e-9, "{m:?}");
            let root = &r.beams[0].forces.ends[0];
            let (shear, moment) = if along_y { (root[1], root[5]) } else { (root[2], root[4]) };
            assert!(rel(shear.abs(), 1000.0) < 1e-9 && rel(moment.abs(), 1000.0 * length) < 1e-9, "{root:?}");
        }
    }
}

#[test]
fn cantilever_beam_frequencies() {
    // A slender pipe, L/r = 1000 (shear and rotary inertia below 10⁻⁴ in the third mode):
    // Euler–Bernoulli ω_k = (β_k L)² √(E I / (ρ A L⁴)).
    let profile = BeamProfile::Pipe { radius: 10.0, thickness: 1.0 };
    let length: f64 = 10000.0;
    let p = fieldcraft_fem::beam::properties(&profile, 0.3);
    let base = (210e3 * p.iy / (7.85e-9 * p.area * length.powi(4))).sqrt() / (2.0 * std::f64::consts::PI);
    let beta = [1.875_104_068_711_961, 4.694_091_132_974_175, 7.854_757_438_237_613];
    for mass in [MassKind::Consistent, MassKind::Lumped] {
        let (b, _, _) = cantilever(40, length, profile, Procedure::Frequency { modes: 6, min_hz: None, max_hz: None });
        let r = solve_frequency(&b.model, None, &FrequencyOptions { mass, shift: None }, &Limits::default()).unwrap();
        // The pipe bends alike in both planes: modes come in pairs.
        let errors: Vec<f64> = beta.iter().enumerate().map(|(k, bl)| r.modes[2 * k].frequency / (bl * bl * base) - 1.0).collect();
        eprintln!("cantilever frequencies ({mass:?}, 40 elements) vs Euler–Bernoulli: {errors:.5?}");
        let mass_total = 7.85e-9 * p.area * length;
        assert!(rel(r.total_mass, mass_total) < 1e-12, "{} {mass_total}", r.total_mass);
        for (k, e) in errors.iter().enumerate() {
            // Lumped mass converges from below, a few 10⁻³ low at 40 elements.
            let tolerance = if mass == MassKind::Consistent { 5e-4 } else { 5e-3 };
            assert!(e.abs() < tolerance, "mode {k}: {e}");
            assert!(rel(r.modes[2 * k].frequency, r.modes[2 * k + 1].frequency) < 1e-6);
        }
    }
}

#[test]
fn portal_frame_matches_the_slope_deflection_solution() {
    // Columns of height h, a beam of span l, fixed bases, a horizontal load H at the top.
    // Slope-deflection with Timoshenko members (Przemieniecki 1968, §5.6): with
    // a = E I_c / (h (1 + Φ_c)) and b = E I_b / (l (1 + Φ_b)), Φ = 12 E I / (G A_s L²), the end
    // moments are M_ij = k [(4 + Φ) θ_i + (2 − Φ) θ_j − 6 Δ/L]; antisymmetric sway gives the joint
    // rotation θ = 6 a Δ / (h (a (4 + Φ_c) + 6 b)), the storey shear
    // a (6 θ − 12 Δ/h) = −H h / 2, and the base moment M = a ((2 − Φ_c) θ − 6 Δ/h). The columns'
    // axial deformation is left out of the hand calculation: their shortening and lengthening
    // under the overturning moment tilts the girder by about 2 H h² / (l² E A), about 2.6 · 10⁻³
    // of θ here; sway and base moment feel it far less.
    let (h, l, load) = (8000.0, 10000.0, 1000.0);
    let column = BeamProfile::Rectangle { width: 100.0, height: 100.0 };
    let girder = BeamProfile::Rectangle { width: 100.0, height: 200.0 };
    let e = 210e3;
    let mut b = Build::new(Elastic::Isotropic { young: e, poisson: 0.3 }, 7.85e-9, Procedure::Static {});
    let n = 8;
    // Nodes: left column 1..=n+1 (base to top), beam to the right column top, right column down.
    let mut path = Vec::new();
    for i in 0..=n {
        path.push([0.0, 0.0, h * i as f64 / n as f64]);
    }
    for i in 1..=n {
        path.push([l * i as f64 / n as f64, 0.0, h]);
    }
    for i in 1..=n {
        path.push([l, 0.0, h - h * i as f64 / n as f64]);
    }
    for (i, x) in path.iter().enumerate() {
        b.node(i as u32 + 1, *x);
    }
    let count = path.len() as u32;
    let mut columns = Vec::new();
    let mut beam = Vec::new();
    for i in 1..count {
        let el = b.element(ElementKind::Beam2, vec![i, i + 1]);
        if (n as u32 + 1..2 * n as u32 + 1).contains(&i) { beam.push(el) } else { columns.push(el) }
    }
    // Both members bend in the x–z plane: their profile height (local z) must lie in it, so
    // the orientation (local y, the width) is along global y.
    b.section("columns", columns, SectionKind::Beam { profile: column, orientation: [0.0, 1.0, 0.0] });
    b.section("girder", beam, SectionKind::Beam { profile: girder, orientation: [0.0, 1.0, 0.0] });
    b.hold_nodes("bases", [Some(0.0); 6], vec![1, count]);
    b.force("push", vec![n as u32 + 1], [load, 0.0, 0.0]);
    let r = b.solve();
    let g = e / 2.6;
    let (pc, pb) = (fieldcraft_fem::beam::properties(&column, 0.3), fieldcraft_fem::beam::properties(&girder, 0.3));
    let (phi_c, phi_b) = (12.0 * e * pc.iy / (g * pc.shear_z * h * h), 12.0 * e * pb.iy / (g * pb.shear_z * l * l));
    let (a, bb) = (e * pc.iy / (h * (1.0 + phi_c)), e * pb.iy / (l * (1.0 + phi_b)));
    let c = 6.0 * a / (h * (a * (4.0 + phi_c) + 6.0 * bb));
    let sway = -load * h / 2.0 / (a * (6.0 * c - 12.0 / h));
    let theta = c * sway;
    let base = a * ((2.0 - phi_c) * theta - 6.0 * sway / h);
    let top = r.displacement[at(&r, n as u32 + 1)][0];
    let m = r.reaction_moment[at(&r, 1)][1];
    let rot = r.rotation[at(&r, n as u32 + 1)][1];
    eprintln!(
        "portal frame: sway {top:.6e} vs {sway:.6e} ({:.1e}), joint rotation {rot:.6e} vs {theta:.6e} ({:.1e}), base moment {m:.6e} vs {:.6e} ({:.1e})",
        rel(top, sway),
        rel(rot.abs(), theta.abs()),
        base.abs(),
        rel(m.abs(), base.abs())
    );
    assert!(rel(top, sway) < 5e-4 && rel(rot.abs(), theta.abs()) < 3e-3 && rel(m.abs(), base.abs()) < 5e-4);
    assert!((r.reaction_total[0] + load).abs() < 1e-9 * load);
}

// ---------------------------------------------------------------------------------------------
// Mixed meshes
// ---------------------------------------------------------------------------------------------

/// A hex8 column 1 × 1 × 2 (two elements) clamped at z = 0, with four rods of area 1/4 from its
/// top corners up to four loaded nodes: the elongation is the sum of both parts' (ν = 0).
#[test]
fn rods_continue_a_solid_column() {
    let e = 1000.0;
    let mut b = Build::new(Elastic::Isotropic { young: e, poisson: 0.0 }, 1.0, Procedure::Static {});
    let corners = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
    for z in 0..3 {
        for (i, c) in corners.iter().enumerate() {
            b.node(1 + i as u32 + 4 * z, [c[0], c[1], f64::from(z)]);
        }
    }
    for (i, c) in corners.iter().enumerate() {
        b.node(13 + i as u32, [c[0], c[1], 5.0]);
    }
    let hexes: Vec<u32> = (0..2).map(|k| b.element(ElementKind::Hex8, (1..=8).map(|n| n + 4 * k).collect())).collect();
    let rods: Vec<u32> = (0..4).map(|i| b.element(ElementKind::Rod2, vec![9 + i, 13 + i])).collect();
    b.section("column", hexes, SectionKind::Solid {});
    b.section("rods", rods, SectionKind::Rod { area: 0.25 });
    b.hold_nodes("base", [Some(0.0), Some(0.0), Some(0.0), FREE, FREE, FREE], vec![1, 2, 3, 4]);
    // Rods carry no lateral stiffness: hold their ends sideways.
    b.hold_nodes("guides", [Some(0.0), Some(0.0), FREE, FREE, FREE, FREE], vec![13, 14, 15, 16]);
    b.force("pull", vec![13, 14, 15, 16], [0.0, 0.0, 4.0]);
    let r = b.solve();
    let want = 4.0 * 2.0 / (e * 1.0) + 4.0 * 3.0 / (e * 1.0);
    for n in 13..=16 {
        let got = r.displacement[at(&r, n)][2];
        assert!(rel(got, want) < 1e-12, "{got} {want}");
    }
    let rods: Vec<_> = r.elements.iter().filter(|e| e.kind == ElementKind::Rod2).collect();
    assert_eq!(rods.len(), 4);
    assert!(rods.iter().flat_map(|e| &e.nodal).all(|s| rel(fieldcraft_fem::von_mises(s), 4.0) < 1e-12), "{rods:?}");
    // Solid nodes have no rotations.
    assert!(r.rotation.iter().all(|x| *x == [0.0; 3]));
}

/// A hex8 slab with a shell skin on its top face sharing the nodes, stretched by prescribed end
/// displacements: the reaction is E (A_solid + A_shell) ε (ν = 0, uniform strain in both).
#[test]
fn shell_skin_on_a_solid_slab() {
    let (e, length, width, depth, skin) = (1000.0, 4.0, 2.0, 1.0, 0.1);
    let mut b = Build::new(Elastic::Isotropic { young: e, poisson: 0.0 }, 1.0, Procedure::Static {});
    let (nx, ny) = (4usize, 2usize);
    let label = |i: usize, j: usize, k: usize| (1 + i + (nx + 1) * (j + (ny + 1) * k)) as u32;
    for k in 0..2 {
        for j in 0..=ny {
            for i in 0..=nx {
                b.node(label(i, j, k), [length * i as f64 / nx as f64, width * j as f64 / ny as f64, depth * k as f64]);
            }
        }
    }
    let mut hexes = Vec::new();
    let mut quads = Vec::new();
    for j in 0..ny {
        for i in 0..nx {
            let q = |k| vec![label(i, j, k), label(i + 1, j, k), label(i + 1, j + 1, k), label(i, j + 1, k)];
            hexes.push(b.element(ElementKind::Hex8, [q(0), q(1)].concat()));
            quads.push(b.element(ElementKind::Quad4, q(1)));
        }
    }
    b.section("slab", hexes, SectionKind::Solid {});
    b.section("skin", quads, SectionKind::Shell { thickness: skin, offset: 0.0, points: 5 });
    let strain = 1e-3;
    b.hold("left", [Some(0.0), Some(0.0), Some(0.0), FREE, FREE, FREE], |p| p[0] == 0.0);
    b.hold("right", [Some(strain * length), Some(0.0), Some(0.0), FREE, FREE, FREE], |p| p[0] == length);
    let r = b.solve();
    let want = e * (width * depth + width * skin) * strain;
    let pull: f64 = b.nodes_where(|p| p[0] == length).iter().map(|n| r.reaction[at(&r, *n)][0]).sum();
    eprintln!("solid with shell skin: end reaction {pull:.9e} vs {want:.9e}");
    assert!(rel(pull, want) < 1e-10, "{pull}");
    let top = r.shells[0].stress.middle[0][0];
    assert!(rel(top, e * strain) < 1e-10, "{top}");
}

// ---------------------------------------------------------------------------------------------
// Errors, not crashes
// ---------------------------------------------------------------------------------------------

#[test]
fn structural_models_that_cannot_be_solved_are_errors() {
    // A beam whose orientation is along it.
    let (mut b, _, _) = cantilever(2, 100.0, BeamProfile::Circle { radius: 1.0 }, Procedure::Static {});
    b.model
        .replace_section(
            "beam",
            Section {
                orientation: None,
                set: "beam".into(),
                material: "Mat".into(),
                kind: SectionKind::Beam { profile: BeamProfile::Circle { radius: 1.0 }, orientation: [2.0, 2.0, 1.0] },
                formulation: Formulation::default(),
            },
        )
        .unwrap();
    b.force("tip", vec![3], [1.0, 0.0, 0.0]);
    let err = solve_static(&b.model, None, &Limits::default()).unwrap_err();
    assert!(err.to_string().contains("parallel"), "{err}");

    // Degenerate first- and second-order shells name the element.
    for (kind, nodes, want) in [(ElementKind::Quad8, 8, "normal"), (ElementKind::Quad4, 4, "degenerate")] {
        let mut b = Build::new(Elastic::Isotropic { young: 1.0, poisson: 0.0 }, 1.0, Procedure::Static {});
        for i in 0..nodes {
            b.node(i as u32 + 1, [f64::from(i as u8 % 2), 0.0, 0.0]);
        }
        let el = b.element(kind, (1..=nodes as u32).collect());
        b.section("s", vec![el], SectionKind::Shell { thickness: 1.0, offset: 0.0, points: 5 });
        let err = solve_static(&b.model, None, &Limits::default()).unwrap_err();
        assert!(err.to_string().contains(want), "{kind:?}: {err}");
    }

    // A moment on the nodes of a rod: they have no rotations.
    let mut b = Build::new(Elastic::Isotropic { young: 1.0, poisson: 0.0 }, 1.0, Procedure::Static {});
    b.node(1, [0.0; 3]);
    b.node(2, [1.0, 0.0, 0.0]);
    let el = b.element(ElementKind::Rod2, vec![1, 2]);
    b.section("r", vec![el], SectionKind::Rod { area: 1.0 });
    b.hold_nodes("all", [Some(0.0); 6], vec![1]);
    b.model.add_set("end", Set { nodes: vec![2], ..Set::default() }).unwrap();
    b.load("twist", Load::Moment { target: Target::Set("end".into()), vector: [1.0, 0.0, 0.0], distribution: Distribution::Total });
    let err = solve_static(&b.model, None, &Limits::default()).unwrap_err();
    assert!(err.to_string().contains("rotations"), "{err}");
}

fn quadratic_shell_model(procedure: Procedure) -> Build {
    let mut b = Build::new(Elastic::Isotropic { young: 1200.0, poisson: 0.0 }, 1.0, procedure);
    for (i, p) in fieldcraft_fem::basis::Basis::Quad8.nodes().iter().enumerate() {
        b.node(i as u32 + 1, [p[0] + 1.0, 0.5 * (p[1] + 1.0), 0.0]);
    }
    let el = b.element(ElementKind::Quad8, (1..=8).collect());
    b.section("panel", vec![el], SectionKind::Shell { thickness: 0.1, offset: 0.0, points: 5 });
    b.hold_nodes("clamp", [Some(0.0); 6], vec![1, 4, 8]);
    for (node, w) in [(2, 1.0 / 6.0), (3, 1.0 / 6.0), (6, 2.0 / 3.0)] {
        let name = format!("m{node}");
        b.model.add_set(&name, Set { nodes: vec![node], ..Set::default() }).unwrap();
        b.load(&name, Load::Moment { target: Target::Set(name.clone()), vector: [0.0, 1e-4 * w, 0.0], distribution: Distribution::Total });
    }
    b
}

#[test]
fn mitc8_solves_pure_bending_and_refuses_unimplemented_procedures() {
    let mut b = quadratic_shell_model(Procedure::Static {});
    let r = b.solve();
    let tip = at(&r, 6);
    assert!((r.displacement[tip][2] / -0.002 - 1.0).abs() < 1e-7, "{:?}", r.displacement[tip]);
    assert_eq!(r.shells[0].stress.top.len(), 8);
    for s in &r.shells[0].stress.top {
        assert!((s[0] - 0.06).abs() < 1e-8, "{s:?}");
    }
    let mut section = b.model.sections()["panel"].clone();
    section.formulation.integration = fieldcraft_model::Integration::Reduced;
    b.model.replace_section("panel", section).unwrap();
    assert!(solve_static(&b.model, None, &Limits::default()).unwrap_err().to_string().contains("S8R"));
    let b = quadratic_shell_model(Procedure::Frequency { modes: 1, min_hz: None, max_hz: None });
    assert!(solve_frequency(&b.model, None, &FrequencyOptions::default(), &Limits::default()).unwrap_err().to_string().contains("MITC8 mass"));
    let b = quadratic_shell_model(Procedure::Buckle { modes: 1, shift: None });
    assert!(fieldcraft_solver::solve_buckle(&b.model, None, &Limits::default()).unwrap_err().to_string().contains("MITC8"));
}

#[test]
fn mitc8_shared_directors_and_global_rotations_are_compatible() {
    use fieldcraft_fem::{
        basis::Basis,
        frame::{cross, dot},
    };
    use fieldcraft_solver::system::{Mesh, StructuralKernel};
    let build = |permutation: [usize; 8], curved: bool| {
        let mut b = Build::new(Elastic::Isotropic { young: 1200.0, poisson: 0.0 }, 1.0, Procedure::Static {});
        let mut labels = BTreeMap::new();
        for j in 0..3 {
            for (i, x) in [0.0f64, 0.5, 1.0, 1.75, 2.5].into_iter().enumerate() {
                if j == 1 && i % 2 == 1 {
                    continue;
                }
                let label = labels.len() as u32 + 1;
                labels.insert((i, j), label);
                b.node(label, [x, j as f64 / 2.0, if curved { 0.05 * x.powi(3) } else { 0.0 }]);
            }
        }
        let mut els = Vec::new();
        for i in [0, 2] {
            let n = [(i, 0), (i + 2, 0), (i + 2, 2), (i, 2), (i + 1, 0), (i + 2, 1), (i + 1, 2), (i, 1)].map(|p| labels[&p]);
            let order = if i == 0 { [0, 1, 2, 3, 4, 5, 6, 7] } else { permutation };
            els.push(b.element(ElementKind::Quad8, order.map(|j| n[j]).to_vec()));
        }
        b.section("panel", els, SectionKind::Shell { thickness: 0.1, offset: 0.0, points: 3 });
        let resolved = b.model.resolve().unwrap();
        Mesh::new(&b.model, &resolved, 1000).unwrap()
    };
    let mut reference = None;
    for permutation in [[0, 1, 2, 3, 4, 5, 6, 7], [1, 2, 3, 0, 5, 6, 7, 4], [0, 3, 2, 1, 7, 6, 5, 4]] {
        let mesh = build(permutation, true);
        let mut k = fieldcraft_fem::dense::Matrix::zeros(6 * mesh.xyz.len());
        for el in &mesh.structural {
            let StructuralKernel::Mitc8(shell) = el.kernel().unwrap() else { panic!() };
            let ke = shell.stiffness().unwrap();
            for i in 0..48 {
                for j in 0..48 {
                    k.add(6 * el.nodes[i / 6] + i % 6, 6 * el.nodes[j / 6] + j % 6, ke.get(i, j));
                }
            }
            let omega = [0.3, -0.2, 0.7];
            let u: Vec<_> = el
                .xyz
                .iter()
                .flat_map(|x| {
                    let v = cross(omega, *x);
                    [v[0], v[1], v[2], omega[0], omega[1], omega[2]]
                })
                .collect();
            for row in shell.strain_operator([0.31, -0.27, 0.8], true).unwrap() {
                assert!(row.iter().zip(&u).map(|(a, b)| a * b).sum::<f64>().abs() < 2e-12);
            }
        }
        if let Some(old) = &reference {
            let old: &fieldcraft_fem::dense::Matrix = old;
            let scale = old.data.iter().map(|v| v.abs()).fold(0.0, f64::max);
            assert!(k.data.iter().zip(&old.data).all(|(a, b)| (a - b).abs() < 1e-11 * scale));
        } else {
            reference = Some(k);
        }
        let a = &mesh.structural[0];
        let b = &mesh.structural[1];
        for (i, n) in a.nodes.iter().enumerate() {
            if let Some(j) = b.nodes.iter().position(|m| m == n) {
                let na = a.directors.unwrap()[i];
                let nb = b.directors.unwrap()[j];
                assert!((dot(na, nb).abs() - 1.0).abs() < 1e-14);
            }
        }
        // Independently interpolate each side's thickened common edge with arbitrary
        // shared global rotations. Equality is checked away from the nodes as well.
        let edge = |el: &fieldcraft_solver::system::Structural, y: f64| {
            let local = el.xyz.iter().enumerate().filter(|(_, x)| x[0] == 1.0).map(|(i, _)| i).collect::<Vec<_>>();
            let mut p = [0.0; 3];
            for &i in &local {
                let yn = el.xyz[i][1];
                let w = if yn == 0.0 {
                    2.0 * (y - 0.5) * (y - 1.0)
                } else if yn == 1.0 {
                    2.0 * y * (y - 0.5)
                } else {
                    4.0 * y * (1.0 - y)
                };
                let x = el.xyz[i];
                let n = el.directors.unwrap()[i];
                let sign = if n[2] > 0.0 { 1.0 } else { -1.0 };
                let rotation = [x[1] * 0.1, 0.2 * x[0], -0.3 * x[1]];
                let v = cross(rotation, n);
                for j in 0..3 {
                    p[j] += w * (x[j] + 0.04 * sign * n[j] + 0.04 * sign * v[j]);
                }
            }
            p
        };
        for y in [0.0, 0.23, 0.5, 0.81, 1.0] {
            let x = edge(a, y);
            let z = edge(b, y);
            assert!(x.iter().zip(z).all(|(a, b)| (a - b).abs() < 1e-14));
        }
    }
    // Two-element pure bending patch energy, independent Eh^3/12 integral.
    let mesh = build([0, 1, 2, 3, 4, 5, 6, 7], false);
    let mut got = 0.0;
    for el in &mesh.structural {
        let u: Vec<_> = el.xyz.iter().flat_map(|x| [0.0, 0.0, -0.5 * x[0] * x[0], 0.0, x[0], 0.0]).collect();
        let k = el.kernel().unwrap().stiffness().unwrap();
        got += 0.5 * k.mul_vec(&u).iter().zip(&u).map(|(a, b)| a * b).sum::<f64>();
    }
    assert!((got / (0.5 * 2.5 * 1200.0 * 0.1f64.powi(3) / 12.0) - 1.0).abs() < 1e-10);
    // Coincident midside coordinates with different labels are NOT a common
    // quadratic edge: refuse rather than silently leave the rotation field split.
    let mut bad = Build::new(Elastic::Isotropic { young: 1200.0, poisson: 0.0 }, 1.0, Procedure::Static {});
    for (&label, &x) in mesh.labels.iter().zip(&mesh.xyz) {
        bad.node(label, x);
    }
    let mut ids = Vec::new();
    for (i, el) in mesh.structural.iter().enumerate() {
        let mut labels: Vec<_> = el.nodes.iter().map(|&n| mesh.labels[n]).collect();
        if i == 1 {
            let j = el.xyz.iter().position(|p| p[0] == 1.0 && p[1] == 0.5).unwrap();
            bad.node(99, el.xyz[j]);
            labels[j] = 99;
        }
        ids.push(bad.element(ElementKind::Quad8, labels));
    }
    bad.section("panel", ids, SectionKind::Shell { thickness: 0.1, offset: 0.0, points: 3 });
    let resolved = bad.model.resolve().unwrap();
    assert!(Mesh::new(&bad.model, &resolved, 1000).unwrap_err().to_string().contains("nonconforming"));
    assert_eq!(Basis::Quad8.nodes().len(), 8);
}
