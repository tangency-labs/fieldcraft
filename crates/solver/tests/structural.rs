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
    BeamProfile, BeamTheory, BoundaryCondition, Direction, Distribution, Elastic, ElementFace, ElementKind, Formulation, Load, Material, Model,
    PointMass, Procedure, Section, SectionKind, Set, Spring, Step, Target,
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
        self.model.add_section(name, Section { set: name.into(), material: "Mat".into(), kind, formulation: Formulation::default() }).unwrap();
    }

    /// Makes the beam section `name` shear-rigid (Euler–Bernoulli).
    fn euler_bernoulli(&mut self, name: &str) {
        let mut section = self.model.section(name).unwrap().clone();
        section.formulation.beam_theory = BeamTheory::EulerBernoulli;
        self.model.replace_section(name, section).unwrap();
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
    let p = profile.properties(0.3);
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
    let p = profile.properties(0.3);
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
    let (pc, pb) = (column.properties(0.3), girder.properties(0.3));
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

#[test]
fn euler_bernoulli_cantilever_deflection_is_exact() {
    // The stocky I beam of `cantilever_beam_deflection_is_exact` with a shear-rigid section:
    // δ = P L³/(3 E I) in each plane, at any number of elements (Przemieniecki 1968, §5.6 with Φ = 0).
    let profile = BeamProfile::I { height: 200.0, width: 100.0, web: 6.0, flange: 10.0 };
    let length: f64 = 1000.0;
    let p = profile.properties(0.3);
    let (e, g) = (210e3, 210e3 / 2.6);
    for n in [1, 4] {
        for (along_y, name) in [(true, "weak"), (false, "strong")] {
            let (mut b, y, z) = cantilever(n, length, profile, Procedure::Static {});
            b.euler_bernoulli("beam");
            let d = if along_y { y } else { z };
            b.force("tip", vec![n as u32 + 1], d.map(|c| 1000.0 * c));
            let r = b.solve();
            let u = r.displacement[at(&r, n as u32 + 1)];
            let got = u[0] * d[0] + u[1] * d[1] + u[2] * d[2];
            let (i, a_s) = if along_y { (p.iz, p.shear_y) } else { (p.iy, p.shear_z) };
            let want = 1000.0 * length.powi(3) / (3.0 * e * i);
            eprintln!(
                "Euler–Bernoulli cantilever {n} element(s), {name} axis: tip {got:.6e} vs {want:.6e} (the shear term left out is {:.1}% of it)",
                100.0 * 1000.0 * length / (g * a_s) / want
            );
            assert!(rel(got, want) < 1e-9, "{got} {want}");
            let m = r.reaction_moment[at(&r, 1)];
            let mm = (m[0] * m[0] + m[1] * m[1] + m[2] * m[2]).sqrt();
            assert!(rel(mm, 1000.0 * length) < 1e-9, "{m:?}");
            let root = &r.beams[0].forces.ends[0];
            let (shear, moment) = if along_y { (root[1], root[5]) } else { (root[2], root[4]) };
            assert!(rel(shear.abs(), 1000.0) < 1e-9 && rel(moment.abs(), 1000.0 * length) < 1e-9, "{root:?}");
        }
    }
}

/// The portal frame of `portal_frame_matches_the_slope_deflection_solution` (fixed bases, `n`
/// elements per member) with half the horizontal load on each joint, so it sways antisymmetrically:
/// the sway, the joints' rotation about Y, the left joint's lift and the left base's moment about Y.
fn sway_frame(n: usize, shear_rigid: bool) -> [f64; 4] {
    let (h, l, load) = (SWAY_FRAME.0, SWAY_FRAME.1, SWAY_FRAME.2);
    let mut b = Build::new(Elastic::Isotropic { young: SWAY_FRAME_YOUNG, poisson: 0.3 }, 7.85e-9, Procedure::Static {});
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
    let (left, right) = (n as u32 + 1, 2 * n as u32 + 1);
    let mut columns = Vec::new();
    let mut beam = Vec::new();
    for i in 1..count {
        let el = b.element(ElementKind::Beam2, vec![i, i + 1]);
        if (left..right).contains(&i) { beam.push(el) } else { columns.push(el) }
    }
    // Both members bend in the x–z plane: the orientation (local y, the width) is along global y.
    b.section("columns", columns, SectionKind::Beam { profile: SWAY_COLUMN, orientation: [0.0, 1.0, 0.0] });
    b.section("girder", beam, SectionKind::Beam { profile: SWAY_GIRDER, orientation: [0.0, 1.0, 0.0] });
    if shear_rigid {
        b.euler_bernoulli("columns");
        b.euler_bernoulli("girder");
    }
    b.hold_nodes("bases", [Some(0.0); 6], vec![1, count]);
    b.force("push", vec![left, right], [load, 0.0, 0.0]);
    let r = b.solve();
    assert!((r.reaction_total[0] + load).abs() < 1e-9 * load);
    [r.displacement[at(&r, left)][0], r.rotation[at(&r, left)][1], r.displacement[at(&r, left)][2], r.reaction_moment[at(&r, 1)][1]]
}

/// Column height, girder span and total horizontal load of [`sway_frame`].
const SWAY_FRAME: (f64, f64, f64) = (8000.0, 10000.0, 1000.0);
const SWAY_FRAME_YOUNG: f64 = 210e3;
const SWAY_COLUMN: BeamProfile = BeamProfile::Rectangle { width: 100.0, height: 100.0 };
const SWAY_GIRDER: BeamProfile = BeamProfile::Rectangle { width: 100.0, height: 200.0 };

#[test]
fn euler_bernoulli_portal_frame_matches_the_slope_deflection_solution() {
    // Slope-deflection without shear (Przemieniecki 1968, §5.6 with Φ = 0; Cook et al. 2002, ch. 2):
    // end moments M_ij = (2 E I / L) (2 θ_i + θ_j − 3 ψ), ψ the chord rotation. Antisymmetric sway:
    // both joints move Δ and turn θ, and the columns stretch and shorten by ±v under the girder's
    // end shears, so the girder's chord turns by 2 v / l. With a = E I_c / h, b = E I_b / l and the
    // column's axial stiffness s = E A_c / h:
    //   column shear      a (6 θ' + 12 Δ / h) = H h / 2   (θ' = −θ, counter-clockwise),
    //   joint moment      (4 a + 6 b) θ' + 6 a Δ / h + 12 b v / l = 0,
    //   joint lift        s v = −(12 b / l) (θ' + 2 v / l),
    // which give v = κ l θ with κ = 12 b / (s l² + 24 b), the girder's restraint softened to
    // b' = b (1 − 2 κ), θ = H h / (4 (a + 6 b')), Δ = h (H h / 2 + 6 a θ) / (12 a) and the base
    // moment 2 a (3 Δ / h − θ). The elements are exact for end loads, so this holds to round-off
    // at any number of elements.
    let (h, l, load) = SWAY_FRAME;
    let (pc, pb) = (SWAY_COLUMN.properties(0.3), SWAY_GIRDER.properties(0.3));
    let (a, b, s) = (SWAY_FRAME_YOUNG * pc.iy / h, SWAY_FRAME_YOUNG * pb.iy / l, SWAY_FRAME_YOUNG * pc.area / h);
    let kappa = 12.0 * b / (s * l * l + 24.0 * b);
    let theta = load * h / (4.0 * (a + 6.0 * b * (1.0 - 2.0 * kappa)));
    let sway = h * (load * h / 2.0 + 6.0 * a * theta) / (12.0 * a);
    let want = [sway, theta, kappa * l * theta, -2.0 * a * (3.0 * sway / h - theta)];
    for n in [1, 3] {
        let got = sway_frame(n, true);
        let errors = [0, 1, 2, 3].map(|k| (got[k] - want[k]).abs() / want[k].abs());
        eprintln!("Euler–Bernoulli portal frame, {n} element(s) per member: sway, rotation, lift, base moment {got:?} vs {want:?} ({errors:?})");
        assert!(errors.iter().all(|e| *e < 1e-9), "{errors:?}");
    }
    // The Timoshenko frame is not this one: its sway is larger by the members' shear (about 5·10⁻⁴).
    let timoshenko = sway_frame(1, false);
    assert!(timoshenko[0] / want[0] - 1.0 > 1e-4, "{} {}", timoshenko[0], want[0]);
}

/// A cantilever along X that bends in the X–Z plane only (every other freedom of every node is
/// held), `n` elements, for a frequency step of three modes.
fn planar_cantilever(n: usize, length: f64, profile: BeamProfile, shear_rigid: bool) -> Build {
    let procedure = Procedure::Frequency { modes: 3, min_hz: None, max_hz: None };
    let mut b = Build::new(Elastic::Isotropic { young: 210e3, poisson: 0.3 }, 7.85e-9, procedure);
    for i in 0..=n {
        b.node(i as u32 + 1, [length * i as f64 / n as f64, 0.0, 0.0]);
    }
    let elements: Vec<u32> = (1..=n as u32).map(|i| b.element(ElementKind::Beam2, vec![i, i + 1])).collect();
    b.section("beam", elements, SectionKind::Beam { profile, orientation: [0.0, 1.0, 0.0] });
    if shear_rigid {
        b.euler_bernoulli("beam");
    }
    b.hold_nodes("clamp", [Some(0.0); 6], vec![1]);
    b.hold_nodes("plane", [Some(0.0), Some(0.0), FREE, Some(0.0), FREE, Some(0.0)], (2..=n as u32 + 1).collect());
    b
}

#[test]
fn euler_bernoulli_cantilever_frequencies() {
    // A stocky beam (L/h = 10), where Timoshenko's third frequency is a tenth lower: the
    // shear-rigid element has no shear and no rotary inertia in bending, so its frequencies are
    // Euler–Bernoulli's ω_k = (β_k L)² √(E I / (ρ A L⁴)) at any slenderness.
    let profile = BeamProfile::Rectangle { width: 50.0, height: 100.0 };
    let length: f64 = 1000.0;
    let p = profile.properties(0.3);
    let base = (210e3 * p.iy / (7.85e-9 * p.area * length.powi(4))).sqrt() / (2.0 * std::f64::consts::PI);
    let beta = [1.875_104_068_711_961, 4.694_091_132_974_175, 7.854_757_438_237_613];
    let errors = |n: usize, mass: MassKind, shear_rigid: bool| -> [f64; 3] {
        let b = planar_cantilever(n, length, profile, shear_rigid);
        let r = solve_frequency(&b.model, None, &FrequencyOptions { mass, shift: None }, &Limits::default()).unwrap();
        [0, 1, 2].map(|k| r.modes[k].frequency / (beta[k] * beta[k] * base) - 1.0)
    };
    // Consistent mass: from above, fourth order in the element length. The leading term is
    // (β_k L / n)⁴ / 1440 (measured: 6.9·10⁻⁴ to 6.94·10⁻⁴ of (β_k L / n)⁴ in all three modes).
    let (coarse, fine) = (errors(10, MassKind::Consistent, true), errors(20, MassKind::Consistent, true));
    eprintln!("Euler–Bernoulli cantilever frequencies, consistent mass: 10 elements {coarse:?}, 20 elements {fine:?}");
    for k in 0..3 {
        let leading = (beta[k] / 20.0).powi(4) / 1440.0;
        assert!(fine[k] > 0.0 && fine[k] < 1.1 * leading, "mode {k}: {} against {leading}", fine[k]);
        let rate = (coarse[k] / fine[k]).log2();
        assert!((3.9..4.1).contains(&rate), "mode {k}: order {rate}");
    }
    // Lumped mass: from below, second order (a few 10⁻³ low at 40 elements, as the Timoshenko beam's).
    let (coarse, fine) = (errors(20, MassKind::Lumped, true), errors(40, MassKind::Lumped, true));
    eprintln!("Euler–Bernoulli cantilever frequencies, lumped mass: 20 elements {coarse:?}, 40 elements {fine:?}");
    for k in 0..3 {
        assert!(fine[k] < 0.0 && fine[k].abs() < 5e-3, "mode {k}: {}", fine[k]);
        let rate = (coarse[k] / fine[k]).log2();
        assert!((1.9..2.1).contains(&rate), "mode {k}: order {rate}");
    }
    // The Timoshenko beam on the same mesh is not Euler–Bernoulli's: its third frequency is a tenth lower.
    let timoshenko = errors(20, MassKind::Consistent, false);
    assert!(timoshenko[2] < -0.05, "{timoshenko:?}");
}

// ---------------------------------------------------------------------------------------------
// Pipes under pressure
// ---------------------------------------------------------------------------------------------

/// Lamé's constant A = (p_i a² − p_o b²)/(b² − a²) of a pipe under pressure (S. P. Timoshenko,
/// J. N. Goodier, *Theory of Elasticity*, 3rd ed., 1970, §28): half of σ_r + σ_θ anywhere in the
/// wall, and the wall's axial stress under the closed-end force.
fn lame_mean(profile: BeamProfile, internal: f64, external: f64) -> f64 {
    let (a, b) = profile.bore().unwrap();
    (internal * a * a - external * b * b) / (b * b - a * a)
}

#[test]
fn a_closed_pipe_under_pressure_follows_lame() {
    // a = 10, b = 20, p_i = 100, p_o = 40: A = −20 and B = 8000, so σ_r = −100 and −40 and
    // σ_θ = 60 and 0 at the bore and outside. Free to stretch, the closed pipe carries σ_x = A
    // and lengthens by (1 − 2ν) A/E; held at both ends it keeps σ_x = 2νA and its anchors
    // carry the effective force, (1 − 2ν) A A_w along the axis at the first.
    let profile = BeamProfile::Pipe { radius: 20.0, thickness: 10.0 };
    let (length, e, nu): (f64, f64, f64) = (1000.0, 210e3, 0.3);
    let (internal, external) = (100.0, 40.0);
    let mean = lame_mean(profile, internal, external);
    assert!((mean + 20.0).abs() < 1e-12, "{mean}");
    let area = profile.area();
    let scale = (mean * area).abs();
    let dir = [2.0 / 3.0, 2.0 / 3.0, 1.0 / 3.0];
    let along = |v: [f64; 3]| v[0] * dir[0] + v[1] * dir[1] + v[2] * dir[2];
    for held in [false, true] {
        let (mut b, _, _) = cantilever(4, length, profile, Procedure::Static {});
        if held {
            b.hold_nodes("end", [Some(0.0); 6], vec![5]);
        }
        b.load("steam", Load::PipePressure { set: "beam".into(), internal, external });
        let r = b.solve();
        let (stretch, wall, anchor) = if held {
            (0.0, 2.0 * nu * mean * area, (1.0 - 2.0 * nu) * mean * area)
        } else {
            (length * (1.0 - 2.0 * nu) * mean / e, mean * area, 0.0)
        };
        let (tip, pull) = (along(r.displacement[at(&r, 5)]), along(r.reaction[at(&r, 1)]));
        eprintln!("closed pipe, held {held}: stretch {tip:.9e} vs {stretch:.9e}, anchor {pull:.9e} vs {anchor:.9e}");
        assert!((tip - stretch).abs() < 1e-9 * (length * mean / e).abs(), "{tip} {stretch}");
        assert!((pull - anchor).abs() < 1e-9 * scale, "{pull} {anchor}");
        assert_eq!(r.beams.len(), 4);
        for beam in &r.beams {
            for end in beam.forces.ends {
                assert!((end[0] - wall).abs() < 1e-9 * scale && end[1].abs() < 1e-9 * scale && end[2].abs() < 1e-9 * scale, "{end:?} {wall}");
            }
            let pipe = beam.pipe.expect("a pipe under pressure reports its wall stresses");
            let want = [-100.0, -40.0, 60.0, 0.0, wall / area, wall / area];
            for (got, want) in pipe.radial.iter().chain(&pipe.hoop).chain(&pipe.axial).zip(want) {
                assert!((got - want).abs() < 1e-9 * 100.0, "{pipe:?}");
            }
        }
    }
    // Two pressures on the same pipes add up; without one a pipe reports no wall stress.
    let (mut b, _, _) = cantilever(2, length, profile, Procedure::Static {});
    b.load("inside", Load::PipePressure { set: "beam".into(), internal, external: 0.0 });
    b.load("outside", Load::PipePressure { set: "beam".into(), internal: 0.0, external });
    let r = b.solve();
    assert!((r.beams[0].pipe.unwrap().hoop[0] - 60.0).abs() < 1e-9 * 100.0, "{:?}", r.beams[0].pipe);
    assert!((r.beams[1].forces.ends[1][0] - mean * area).abs() < 1e-9 * scale, "{:?}", r.beams[1].forces);
    let (mut b, _, _) = cantilever(2, length, profile, Procedure::Static {});
    b.force("tip", vec![3], [1.0, 0.0, 0.0]);
    assert!(b.solve().beams.iter().all(|beam| beam.pipe.is_none()));
    // A pressure on a beam that is no pipe is refused by the model's check.
    let (mut b, _, _) = cantilever(2, length, BeamProfile::Circle { radius: 20.0 }, Procedure::Static {});
    b.load("steam", Load::PipePressure { set: "beam".into(), internal, external });
    let err = solve_static(&b.model, None, &Limits::default()).unwrap_err();
    assert!(err.to_string().contains("pipe profile that has a bore"), "{err}");
}

#[test]
fn a_pressurised_pipe_loop_between_two_anchors() {
    // A pipe loop in the x–z plane: two legs of height h up from anchors a span w apart, joined
    // by a run, under internal pressure. Unrestrained, every member would lengthen by the
    // strain ε = (1 − 2ν) A/E (A Lamé's constant), a similar loop; the anchors hold the span,
    // so the run's growth ε w is taken up by the thrust X between the anchors, which acts at
    // the height ȳ of the loop's elastic centre (the centroid of its centre line; no moment and
    // no vertical force there, by the centroid and by symmetry). With bending, axial and shear
    // flexibility (Castigliano; Timoshenko members, Przemieniecki 1968, §5.6):
    // X (I_c/(E I) + w/(E A) + 2h/(G A_s)) = ε w, I_c = ∫ (y − ȳ)² ds over the centre line.
    // The legs carry no effective force, so the corners rise by ε h; the run is compressed by
    // X and bent by the constant moment X (h − ȳ), which lifts its middle by κ w²/8.
    let profile = BeamProfile::Pipe { radius: 50.0, thickness: 5.0 };
    let (h, w, internal): (f64, f64, f64) = (2000.0, 3000.0, 10.0);
    let (e, nu): (f64, f64) = (210e3, 0.3);
    let g = e / (2.0 * (1.0 + nu));
    let p = profile.properties(nu);
    let mean = lame_mean(profile, internal, 0.0);
    let (end_force, strain) = (mean * p.area, (1.0 - 2.0 * nu) * mean / e);
    let mut b = Build::new(Elastic::Isotropic { young: e, poisson: nu }, 7.85e-9, Procedure::Static {});
    let path = [[0.0, 0.0, 0.0], [0.0, 0.0, h / 2.0], [0.0, 0.0, h], [w / 2.0, 0.0, h], [w, 0.0, h], [w, 0.0, h / 2.0], [w, 0.0, 0.0]];
    for (i, x) in path.iter().enumerate() {
        b.node(i as u32 + 1, *x);
    }
    let elements: Vec<u32> = (1..7).map(|i| b.element(ElementKind::Beam2, vec![i, i + 1])).collect();
    b.section("loop", elements, SectionKind::Beam { profile, orientation: [0.0, 1.0, 0.0] });
    b.hold_nodes("anchors", [Some(0.0); 6], vec![1, 7]);
    b.load("steam", Load::PipePressure { set: "loop".into(), internal, external: 0.0 });
    let r = b.solve();
    let centre = (h * h + w * h) / (2.0 * h + w);
    let inertia = 2.0 * (h.powi(3) / 12.0 + h * (h / 2.0 - centre).powi(2)) + w * (h - centre).powi(2);
    let thrust = strain * w / (inertia / (e * p.iy) + w / (e * p.area) + 2.0 * h / (g * p.shear_z));
    let corner = w * (strain - thrust / (e * p.area)) / 2.0;
    let crown = strain * h + thrust * (h - centre) / (e * p.iy) * w * w / 8.0;
    let got = [
        ("thrust at the first anchor", r.reaction[at(&r, 1)][0], thrust),
        ("thrust at the second anchor", -r.reaction[at(&r, 7)][0], thrust),
        ("anchor moment", r.reaction_moment[at(&r, 1)][1].abs(), thrust * centre),
        ("corner 5 along x", r.displacement[at(&r, 5)][0], corner),
        ("corner 3 along x", -r.displacement[at(&r, 3)][0], corner),
        ("corner 3 along z", r.displacement[at(&r, 3)][2], strain * h),
        ("crown along z", r.displacement[at(&r, 4)][2], crown),
        ("wall force in the run", r.beams[2].forces.ends[1][0], end_force - thrust),
        ("wall force in a leg", r.beams[0].forces.ends[0][0], end_force),
    ];
    for (what, value, want) in got {
        eprintln!("pipe loop, {what}: {value:.12e} vs {want:.12e} ({:.1e})", rel(value, want));
        // The thrust is 1.1 · 10⁻³ of the pressure pair it balances, so its round-off is that much larger.
        assert!(rel(value, want) < 1e-8, "{what}: {value} {want}");
    }
    assert!(r.reaction_total.iter().all(|v| v.abs() < 1e-9 * end_force), "{:?}", r.reaction_total);
    let hoop = r.beams[3].pipe.unwrap().hoop;
    let (a2, b2) = (45.0 * 45.0, 50.0 * 50.0);
    assert!(rel(hoop[0], internal * (a2 + b2) / (b2 - a2)) < 1e-12 && rel(hoop[1], 2.0 * internal * a2 / (b2 - a2)) < 1e-12, "{hoop:?}");
}

// ---------------------------------------------------------------------------------------------
// Beams in 2D models
// ---------------------------------------------------------------------------------------------

/// A plane stress quad4 on four nodes (counter-clockwise) held by `values`: it makes the model a
/// 2D one (a model is 2D by its 2D sections) and, held at every node, carries nothing; the
/// frames below stand on it.
fn foundation(b: &mut Build, nodes: [u32; 4], values: [Option<f64>; 6]) {
    let plate = b.element(ElementKind::Quad4, nodes.to_vec());
    b.section("foundation", vec![plate], SectionKind::PlaneStress { thickness: 1.0 });
    b.hold_nodes("ground", values, nodes.to_vec());
}

#[test]
fn planar_cantilever_beam_is_exact() {
    // The cantilever of `cantilever_beam_deflection_is_exact` in a 2D model: three unknowns per
    // beam node, the profile's height in the plane (orientation along Z), so the plane bends
    // with I_y and the shear area along z. Tip shear: δ = P L³/(3 E I) + P L/(G A_s); tip
    // tension: P L/(E A); tip moment: θ = M L/(E I), δ = M L²/(2 E I).
    let profile = BeamProfile::I { height: 200.0, width: 100.0, web: 6.0, flange: 10.0 };
    let (length, load): (f64, f64) = (1000.0, 1000.0);
    let p = profile.properties(0.3);
    let (e, g) = (210e3, 210e3 / 2.6);
    // The axis and the in-plane direction a quarter turn counter-clockwise from it.
    let (dir, across) = ([0.6, 0.8, 0.0], [-0.8, 0.6, 0.0]);
    for n in [1usize, 4] {
        for case in ["shear", "tension", "moment"] {
            let mut b = Build::new(Elastic::Isotropic { young: e, poisson: 0.3 }, 7.85e-9, Procedure::Static {});
            for i in 0..=n {
                let s = length * i as f64 / n as f64;
                b.node(i as u32 + 1, dir.map(|c| c * s));
            }
            let elements: Vec<u32> = (1..=n as u32).map(|i| b.element(ElementKind::Beam2, vec![i, i + 1])).collect();
            b.section("beam", elements, SectionKind::Beam { profile, orientation: [0.0, 0.0, 1.0] });
            for (label, xy) in [(901, [-100.0, -100.0]), (902, [0.0, -100.0]), (903, [-100.0, 0.0])] {
                b.node(label, [xy[0], xy[1], 0.0]);
            }
            foundation(&mut b, [901, 902, 1, 903], [Some(0.0); 6]);
            let tip = n as u32 + 1;
            match case {
                "shear" => b.force("tip", vec![tip], across.map(|c| load * c)),
                "tension" => b.force("tip", vec![tip], dir.map(|c| load * c)),
                _ => {
                    b.model.add_set("tip", Set { nodes: vec![tip], ..Set::default() }).unwrap();
                    b.load("tip", Load::Moment { target: Target::Set("tip".into()), vector: [0.0, 0.0, load], distribution: Distribution::Total });
                }
            }
            let r = b.solve();
            assert_eq!(r.stats.dofs, 3 * (n + 1) + 2 * 3, "three unknowns per beam node, two per plate node");
            let (u, theta) = (r.displacement[at(&r, tip)], r.rotation[at(&r, tip)]);
            let along = |v: [f64; 3]| u[0] * v[0] + u[1] * v[1];
            let moment = r.reaction_moment[at(&r, 1)];
            let root = &r.beams[0].forces.ends[0];
            eprintln!("planar cantilever, {n} element(s), {case}: tip {u:?}, rotation {theta:?}, clamp moment {moment:?}, root {root:?}");
            match case {
                "shear" => {
                    let want = load * length.powi(3) / (3.0 * e * p.iy) + load * length / (g * p.shear_z);
                    assert!(rel(along(across), want) < 1e-9, "{} {want}", along(across));
                    assert!(rel(moment[2], -load * length) < 1e-9, "{moment:?}");
                    assert!(rel(root[2].abs(), load) < 1e-9 && rel(root[4].abs(), load * length) < 1e-9, "{root:?}");
                }
                "tension" => {
                    assert!(rel(along(dir), load * length / (e * p.area)) < 1e-9, "{}", along(dir));
                    assert!(rel(root[0], load) < 1e-9, "{root:?}");
                }
                _ => {
                    assert!(rel(theta[2], load * length / (e * p.iy)) < 1e-9, "{theta:?}");
                    assert!(rel(along(across), load * length * length / (2.0 * e * p.iy)) < 1e-9, "{}", along(across));
                    assert!(rel(moment[2], -load) < 1e-9, "{moment:?}");
                }
            }
            // Nothing leaves the plane.
            assert!(u[2] == 0.0 && theta[0] == 0.0 && theta[1] == 0.0 && moment[0] == 0.0 && moment[1] == 0.0, "{u:?} {theta:?} {moment:?}");
            assert!(root[1] == 0.0 && root[3] == 0.0 && root[5] == 0.0, "{root:?}");
        }
    }
}

#[test]
fn planar_portal_frame_matches_the_slope_deflection_solution() {
    // The frame of `portal_frame_matches_the_slope_deflection_solution` in the XY plane of a 2D
    // model (columns along Y, the girder along X), with the same hand calculation and
    // tolerances.
    let (h, l, load) = (8000.0, 10000.0, 1000.0);
    let column = BeamProfile::Rectangle { width: 100.0, height: 100.0 };
    let girder = BeamProfile::Rectangle { width: 100.0, height: 200.0 };
    let e = 210e3;
    let mut b = Build::new(Elastic::Isotropic { young: e, poisson: 0.3 }, 7.85e-9, Procedure::Static {});
    let n = 8;
    let mut path = Vec::new();
    for i in 0..=n {
        path.push([0.0, h * i as f64 / n as f64, 0.0]);
    }
    for i in 1..=n {
        path.push([l * i as f64 / n as f64, h, 0.0]);
    }
    for i in 1..=n {
        path.push([l, h - h * i as f64 / n as f64, 0.0]);
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
    b.section("columns", columns, SectionKind::Beam { profile: column, orientation: [0.0, 0.0, 1.0] });
    b.section("girder", beam, SectionKind::Beam { profile: girder, orientation: [0.0, 0.0, 1.0] });
    // The bases stand on the foundation's top corners.
    b.node(901, [0.0, -1000.0, 0.0]);
    b.node(902, [l, -1000.0, 0.0]);
    foundation(&mut b, [901, 902, count, 1], [Some(0.0); 6]);
    b.force("push", vec![n as u32 + 1], [load, 0.0, 0.0]);
    let r = b.solve();
    assert_eq!(r.stats.dofs, 3 * path.len() + 2 * 2);
    let g = e / 2.6;
    let (pc, pb) = (column.properties(0.3), girder.properties(0.3));
    let (phi_c, phi_b) = (12.0 * e * pc.iy / (g * pc.shear_z * h * h), 12.0 * e * pb.iy / (g * pb.shear_z * l * l));
    let (a, bb) = (e * pc.iy / (h * (1.0 + phi_c)), e * pb.iy / (l * (1.0 + phi_b)));
    let c = 6.0 * a / (h * (a * (4.0 + phi_c) + 6.0 * bb));
    let sway = -load * h / 2.0 / (a * (6.0 * c - 12.0 / h));
    let theta = c * sway;
    let base = a * ((2.0 - phi_c) * theta - 6.0 * sway / h);
    let top = r.displacement[at(&r, n as u32 + 1)][0];
    let m = r.reaction_moment[at(&r, 1)][2];
    let rot = r.rotation[at(&r, n as u32 + 1)][2];
    eprintln!(
        "planar portal frame: sway {top:.6e} vs {sway:.6e} ({:.1e}), joint rotation {rot:.6e} vs {theta:.6e} ({:.1e}), base moment {m:.6e} vs {:.6e} ({:.1e})",
        rel(top, sway),
        rel(rot.abs(), theta.abs()),
        base.abs(),
        rel(m.abs(), base.abs())
    );
    assert!(rel(top, sway) < 5e-4 && rel(rot.abs(), theta.abs()) < 3e-3 && rel(m.abs(), base.abs()) < 5e-4);
    assert!((r.reaction_total[0] + load).abs() < 1e-9 * load);
}

/// A plane stress plate with a beam along each loaded edge, sharing its nodes: in uniform
/// tension (ν = 0) the plate and the beams are springs in parallel, u = F L / (E t w + E A).
#[test]
fn plane_beams_stiffen_a_plane_stress_plate() {
    let (l, w, t, young, force) = (4.0, 1.0, 0.2, 1000.0, 30.0);
    let profile = BeamProfile::Rectangle { width: 0.1, height: 0.25 };
    let mut b = Build::new(Elastic::Isotropic { young, poisson: 0.0 }, 1.0, Procedure::Static {});
    for (label, xyz) in [(1, [0.0, 0.0, 0.0]), (2, [l, 0.0, 0.0]), (3, [l, w, 0.0]), (4, [0.0, w, 0.0])] {
        b.node(label, xyz);
    }
    let plate = b.element(ElementKind::Quad4, vec![1, 2, 3, 4]);
    b.section("plate", vec![plate], SectionKind::PlaneStress { thickness: t });
    let stiffeners = vec![b.element(ElementKind::Beam2, vec![1, 2]), b.element(ElementKind::Beam2, vec![4, 3])];
    b.section("stiffeners", stiffeners, SectionKind::Beam { profile, orientation: [0.0, 0.0, 1.0] });
    b.hold_nodes("left", [Some(0.0), FREE, FREE, FREE, FREE, FREE], vec![1, 4]);
    b.hold_nodes("corner", [FREE, Some(0.0), FREE, FREE, FREE, FREE], vec![1]);
    b.model.add_set("right", Set { nodes: vec![2, 3], ..Set::default() }).unwrap();
    b.load("pull", Load::Force { target: Target::Set("right".into()), vector: [force, 0.0, 0.0], distribution: Distribution::PerNode });
    let r = b.solve();
    let u = r.displacement[at(&r, 2)][0];
    let exact = 2.0 * force * l / (young * t * w + young * 2.0 * profile.area());
    let turn = r.rotation.iter().map(|x| x[2].abs()).fold(0.0, f64::max);
    eprintln!("plane beams on a plate: u {u:.12} vs {exact:.12}; largest rotation {turn:.3e}; {} unknowns", r.stats.dofs);
    assert!(rel(u, exact) < 1e-12 && turn < 1e-12 * exact / l, "{u} {exact} {turn}");
    assert_eq!(r.stats.dofs, 3 * 4, "Ux, Uy and Rz at every node");
}

/// The pin-ended double cross of NAFEMS test FV2 (*Selected Benchmarks for Natural Frequency
/// Analysis*, NAFEMS 1987): eight arms of length 5 m at 45° meeting at a rigid joint, each pinned
/// at its far end, square section 0.125 m, E = 200 GPa, ρ = 8000 kg/m³, four elements per arm.
#[test]
fn pin_ended_double_cross_in_plane_frequencies() {
    let (arm, side, young, density): (f64, f64, f64, f64) = (5.0, 0.125, 200e9, 8000.0);
    let profile = BeamProfile::Rectangle { width: side, height: side };
    let n = 4u32;
    let mut b = Build::new(Elastic::Isotropic { young, poisson: 0.3 }, density, Procedure::Frequency { modes: 16, min_hz: None, max_hz: None });
    b.node(1, [0.0; 3]);
    let (mut elements, mut ends) = (Vec::new(), Vec::new());
    for k in 0..8u32 {
        let angle = f64::from(k) * std::f64::consts::FRAC_PI_4;
        let mut previous = 1;
        for i in 1..=n {
            let label = 1 + k * n + i;
            let s = arm * f64::from(i) / f64::from(n);
            b.node(label, [s * angle.cos(), s * angle.sin(), 0.0]);
            elements.push(b.element(ElementKind::Beam2, vec![previous, label]));
            previous = label;
        }
        ends.push(previous);
    }
    b.section("arms", elements, SectionKind::Beam { profile, orientation: [0.0, 0.0, 1.0] });
    // The pins: every arm end held in X and Y, its rotation free. The plate between four of
    // them makes the model 2D and carries nothing.
    let plate = b.element(ElementKind::Quad4, vec![ends[0], ends[2], ends[4], ends[6]]);
    b.section("plate", vec![plate], SectionKind::PlaneStress { thickness: 0.01 });
    b.hold_nodes("pins", [Some(0.0), Some(0.0), FREE, FREE, FREE, FREE], ends.clone());
    let r = solve_frequency(&b.model, None, &FrequencyOptions { mass: MassKind::Consistent, shift: None }, &Limits::default()).unwrap();
    assert_eq!(r.stats.dofs, 3 * (1 + 8 * n as usize));
    let p = profile.properties(0.3);
    let g = young / 2.6;
    let pi = std::f64::consts::PI;
    let euler = |lambda: f64| lambda * lambda / (2.0 * pi * arm * arm) * (young * p.iy / (density * p.area)).sqrt();
    // A pinned–pinned Timoshenko beam's frequency of wave number k = jπ/L, the smaller root of
    // a ω⁴ − b ω² + c = 0.
    let timoshenko = |j: f64| {
        let k = j * pi / arm;
        let a = density * p.iy / (g * p.shear_z);
        let bq = 1.0 + p.iy / p.area * k * k + young * p.iy * k * k / (g * p.shear_z);
        let c = young * p.iy * k.powi(4) / (density * p.area);
        (2.0 * c / (bq + (bq * bq - 4.0 * a * c).sqrt())).sqrt() / (2.0 * pi)
    };
    // Mode 1: every arm a pinned–pinned beam (the joint turns, and by symmetry carries no moment).
    // Four elements per half wave: the cubic element with its consistent mass is (kh)⁴/1440 =
    // 2.6e-4 high at kh = π/4.
    let first = r.modes[0].frequency;
    eprintln!("double cross mode 1: {first:.6} Hz vs Timoshenko {:.6} ({:.1e})", timoshenko(1.0), rel(first, timoshenko(1.0)));
    assert!(rel(first, timoshenko(1.0)) < 1e-3, "{first}");
    // NAFEMS's targets (Hz) by mode, from Euler–Bernoulli theory (λ²/(2π L²) √(E I/ρ A) with
    // λ = π, 3.9266, 2π, 7.0686: pinned–pinned and clamped–pinned arms): shear and rotary inertia,
    // which the element has and the targets leave out, lower them by 0.1 % to 0.5 %.
    assert!(rel(euler(pi), 11.336) < 1e-4 && rel(euler(3.926_602_312_047_919), 17.709) < 1e-4);
    assert!(rel(euler(2.0 * pi), 45.345) < 1e-4 && rel(euler(7.068_582_745_628_732), 57.390) < 1e-4);
    assert_eq!(r.modes.len(), 16);
    for (i, m) in r.modes.iter().enumerate() {
        let target = [11.336, 17.709, 45.345, 57.390][usize::from(i > 0) + usize::from(i > 7) + usize::from(i > 8)];
        eprintln!("double cross mode {}: {:.6} Hz vs {target} ({:+.2e})", i + 1, m.frequency, m.frequency / target - 1.0);
        assert!(rel(m.frequency, target) < 0.01, "mode {}: {} vs {target}", i + 1, m.frequency);
    }
}

#[test]
fn planar_beam_models_that_cannot_be_solved_are_errors() {
    // A beam from a plate's corner; the plate's base held, and the corner's rotation.
    let build = |orientation: [f64; 3], z: f64, plate: SectionKind| {
        let mut b = Build::new(Elastic::Isotropic { young: 1.0, poisson: 0.0 }, 1.0, Procedure::Static {});
        for (label, xyz) in [(1, [0.0, 0.0, 0.0]), (2, [1.0, 0.0, 0.0]), (3, [0.5, 1.0, z]), (901, [0.0, -1.0, 0.0]), (902, [1.0, -1.0, 0.0])] {
            b.node(label, xyz);
        }
        let beam = b.element(ElementKind::Beam2, vec![1, 3]);
        b.section("beam", vec![beam], SectionKind::Beam { profile: BeamProfile::Circle { radius: 0.1 }, orientation });
        let quad = b.element(ElementKind::Quad4, vec![901, 902, 2, 1]);
        b.section("plate", vec![quad], plate);
        b.hold_nodes("base", [Some(0.0); 6], vec![901, 902]);
        b.hold_nodes("turn", [FREE, FREE, FREE, FREE, FREE, Some(0.0)], vec![1]);
        b.force("tip", vec![3], [1.0, 0.0, 0.0]);
        b
    };
    let stress = SectionKind::PlaneStress { thickness: 1.0 };
    let error = |b: &Build| solve_static(&b.model, None, &Limits::default()).unwrap_err().to_string();
    // Either way along Z is the same section.
    assert!(solve_static(&build([0.0, 0.0, -2.0], 0.0, stress).model, None, &Limits::default()).is_ok());
    for (b, want) in [
        (build([0.0, 1.0, 0.0], 0.0, stress), "along Z"),
        (build([0.0, 0.0, 1.0], 0.5, stress), "off the XY plane"),
        (build([0.0, 0.0, 1.0], 0.0, SectionKind::Axisymmetric { twist: false }), "cannot be mixed"),
    ] {
        let e = error(&b);
        assert!(e.contains(want), "{want}: {e}");
    }
    let mut b = build([0.0, 0.0, 1.0], 0.0, stress);
    b.model.add_set("end", Set { nodes: vec![3], ..Set::default() }).unwrap();
    b.load("twist", Load::Moment { target: Target::Set("end".into()), vector: [1.0, 0.0, 0.0], distribution: Distribution::Total });
    let e = error(&b);
    assert!(e.contains("about the Z axis"), "{e}");
    // The planar beam is a Timoshenko beam: an Euler–Bernoulli section is refused by name, and a
    // pipe pressure has no pipe wall to act on.
    let mut b = build([0.0, 0.0, 1.0], 0.0, stress);
    b.euler_bernoulli("beam");
    let e = error(&b);
    assert!(e.contains("section `beam`: an Euler-Bernoulli beam is not supported in a 2D model"), "{e}");
    let mut b = build([0.0, 0.0, 1.0], 0.0, stress);
    b.load("steam", Load::PipePressure { set: "beam".into(), internal: 1.0, external: 0.0 });
    let e = error(&b);
    assert!(e.contains("a pipe pressure is not supported in a 2D model"), "{e}");
}

/// A cantilever beam along X in a 2D model (node 1 clamped, node 2 its tip) on a held plate.
fn planar_beam_cantilever(profile: BeamProfile, length: f64, procedure: Procedure) -> Build {
    let mut b = Build::new(Elastic::Isotropic { young: 1.0, poisson: 0.3 }, 1e-9, procedure);
    b.node(1, [0.0; 3]);
    b.node(2, [length, 0.0, 0.0]);
    let beam = b.element(ElementKind::Beam2, vec![1, 2]);
    b.section("beam", vec![beam], SectionKind::Beam { profile, orientation: [0.0, 0.0, 1.0] });
    for (label, xyz) in [(901, [0.0, -1.0, 0.0]), (902, [1.0, -1.0, 0.0]), (903, [1.0, -2.0, 0.0]), (904, [0.0, -2.0, 0.0])] {
        b.node(label, xyz);
    }
    foundation(&mut b, [904, 903, 902, 901], [Some(0.0); 6]);
    b.hold_nodes("clamp", [Some(0.0); 6], vec![1]);
    b
}

/// A rotary inertia I33 at the tip of a planar cantilever beam acts on Rz alone. The tip's
/// translations carry no mass, so they condense out exactly: under a pure tip moment the tip turns
/// by M L / (E I), a rotational stiffness E I / L whatever the shear flexibility, and the one
/// mode of the disc is f = sqrt(E I / (L J)) / 2π. The beam's own mass (density 1e-9) is a
/// part in 1e6 of the disc's.
#[test]
fn rotary_inertia_at_a_planar_beam_node_acts_about_z() {
    let profile = BeamProfile::Circle { radius: 0.1 };
    let (length, disc) = (1.0, 1.0);
    let mut b = planar_beam_cantilever(profile, length, Procedure::Frequency { modes: 1, min_hz: None, max_hz: None });
    b.model.add_set("tip", Set { nodes: vec![2], ..Set::default() }).unwrap();
    b.model.add_mass("disc", PointMass { target: Target::Set("tip".into()), mass: None, inertia: Some([0.0, 0.0, disc, 0.0, 0.0, 0.0]) }).unwrap();
    let r = solve_frequency(&b.model, None, &FrequencyOptions { mass: MassKind::Consistent, shift: None }, &Limits::default()).unwrap();
    let exact = (profile.properties(0.3).iy / (length * disc)).sqrt() / (2.0 * std::f64::consts::PI);
    let f = r.modes[0].frequency;
    eprintln!("disc on a planar cantilever: {f:.9e} Hz vs {exact:.9e} ({:+.1e})", f / exact - 1.0);
    assert!(rel(f, exact) < 1e-5, "{f} vs {exact}");
    assert_eq!(r.stats.dofs, 3 * 2 + 2 * 4, "Ux, Uy and Rz at the beam's nodes, Ux and Uy at the plate's");
}

/// What leaves the plane of a 2D model is refused by name: an inertia about X or Y (or a product),
/// a torsional spring about X or Y, a spring along Z.
#[test]
fn discrete_elements_that_leave_the_plane_of_a_2d_model_are_refused() {
    let profile = BeamProfile::Circle { radius: 0.1 };
    let procedure = || Procedure::Frequency { modes: 1, min_hz: None, max_hz: None };
    let solve = |b: &Build| {
        solve_frequency(&b.model, None, &FrequencyOptions { mass: MassKind::Consistent, shift: None }, &Limits::default()).unwrap_err().to_string()
    };
    for (inertia, what) in [([1.0, 1.0, 1.0, 0.0, 0.0, 0.0], "I11, I22"), ([1.0, 1.0, 1.0, 0.0, 0.0, 0.5], "I23")] {
        let mut b = planar_beam_cantilever(profile, 1.0, procedure());
        b.model.add_set("tip", Set { nodes: vec![2], ..Set::default() }).unwrap();
        b.model.add_mass("disc", PointMass { target: Target::Set("tip".into()), mass: None, inertia: Some(inertia) }).unwrap();
        let e = solve(&b);
        assert!(e.contains("point mass `disc`") && e.contains("2D model") && e.contains("about Z only"), "{what}: {e}");
    }
    for dof in [3u8, 4, 5] {
        let mut b = planar_beam_cantilever(profile, 1.0, procedure());
        b.model.add_set("tip", Set { nodes: vec![2], ..Set::default() }).unwrap();
        b.model
            .add_spring("k", Spring::Spring { first: Target::Set("tip".into()), second: None, direction: Direction::Dof(dof), stiffness: 1.0 })
            .unwrap();
        let e = solve(&b);
        assert!(e.contains("`k`") && e.contains("2D model") && e.contains("leaves the plane"), "dof {dof}: {e}");
    }
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

    // A quad8 shell is not supported yet; a degenerate quad4 is an error naming the element.
    for (kind, nodes, want) in [(ElementKind::Quad8, 8, "not supported"), (ElementKind::Quad4, 4, "degenerate")] {
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
