//! Benchmark suite: shell problems compared with reference values published in papers, each
//! tagged with its source. These are not closed forms, so they are kept apart from the analytic
//! verification suite (`structural.rs`, `verification.rs`), which is the gate; a benchmark checks
//! convergence towards the published value within a stated margin.
//!
//! - **Scordelis–Lo roof** [benchmark: MacNeal & Harder 1985]: a cylindrical roof, R = 25,
//!   L = 50, t = 0.25, 80° arc, E = 4.32 · 10⁸, ν = 0, self-weight 90 per area, rigid
//!   diaphragms at the curved ends, straight edges free. Vertical deflection at the middle of a
//!   free edge 0.3024 (R. H. MacNeal, R. L. Harder, *A proposed standard set of problems to test
//!   finite element accuracy*, Finite Elements in Analysis and Design 1(1):3–20, 1985).
//! - **Pinched cylinder** [benchmark: Belytschko et al. 1985]: R = 300, L = 600, t = 3,
//!   E = 3 · 10⁶, ν = 0.3, rigid diaphragms at both ends, two opposite unit forces at the middle.
//!   Deflection under a load 1.8248 · 10⁻⁵ (T. Belytschko, H. Stolarski, W. K. Liu, N. Carpenter,
//!   J. S.-J. Ong, *Stress projection for membrane and shear locking in shell finite elements*,
//!   CMAME 51(1–3):221–258, 1985).
//!
//! Both are modelled on one symmetric part (a quarter of the roof, an eighth of the cylinder)
//! with quad4 (MITC4) shells on a structured mesh whose quadrilaterals are exactly flat.

use fieldcraft_model::{
    BoundaryCondition, Distribution, Elastic, ElementKind, Formulation, Load, Material, Model, Procedure, Section, SectionKind, Set, Step, Symmetry,
    Target,
};
use fieldcraft_solver::{Limits, solve_static};

/// A cylindrical shell patch: x from 0 to `length`, angle φ from 0 (top, z = R) to `arc`,
/// y = R sin φ, z = R cos φ; n × n quad4 cells. Returns the model and the node at (x, φ).
struct Cylinder {
    model: Model,
    n: usize,
}

impl Cylinder {
    fn new(n: usize, radius: f64, length: f64, arc: f64, thickness: f64, elastic: Elastic) -> Cylinder {
        let mut model = Model::default();
        let label = |i: usize, j: usize| (1 + i + (n + 1) * j) as u32;
        for j in 0..=n {
            for i in 0..=n {
                let (x, phi) = (length * i as f64 / n as f64, arc * j as f64 / n as f64);
                model.set_node(label(i, j), [x, radius * phi.sin(), radius * phi.cos()]).unwrap();
            }
        }
        let mut elements = Vec::new();
        for j in 0..n {
            for i in 0..n {
                let e = (1 + i + n * j) as u32;
                // Counter-clockwise seen from outside: the normal points away from the axis.
                model.set_element(e, ElementKind::Quad4, vec![label(i, j), label(i + 1, j), label(i + 1, j + 1), label(i, j + 1)]).unwrap();
                elements.push(e);
            }
        }
        model.add_set("shell", Set { elements, ..Set::default() }).unwrap();
        model.add_material("Mat", Material { elastic: Some(elastic), density: None, ..Material::default() }).unwrap();
        model
            .add_section(
                "shell",
                Section {
                    orientation: None,
                    set: "shell".into(),
                    material: "Mat".into(),
                    kind: SectionKind::Shell { thickness, offset: 0.0, points: 5 },
                    formulation: Formulation::default(),
                },
            )
            .unwrap();
        model.add_step(Step::new("Step", Procedure::Static {}), None).unwrap();
        Cylinder { model, n }
    }

    fn label(&self, i: usize, j: usize) -> u32 {
        (1 + i + (self.n + 1) * j) as u32
    }

    fn nodes(&self, f: impl Fn(usize, usize) -> bool) -> Vec<u32> {
        let mut out = Vec::new();
        for j in 0..=self.n {
            for i in 0..=self.n {
                if f(i, j) {
                    out.push(self.label(i, j));
                }
            }
        }
        out.sort_unstable();
        out
    }

    fn bc(&mut self, name: &str, nodes: Vec<u32>, bc: impl Fn(Target) -> BoundaryCondition) {
        self.model.add_set(name, Set { nodes, ..Set::default() }).unwrap();
        self.model.add_boundary_condition(name, bc(Target::Set(name.into()))).unwrap();
        self.model.set_boundary_condition_active("Step", name, true).unwrap();
    }

    /// A rigid diaphragm at x = length: no motion in its plane (y, z) and no rotation about x.
    fn diaphragm(&mut self) {
        let n = self.n;
        let values = [None, Some(0.0), Some(0.0), Some(0.0), None, None];
        self.bc("diaphragm", self.nodes(|i, _| i == n), |target| BoundaryCondition::Displacement { target, values });
    }

    fn symmetry(&mut self, name: &str, plane: Symmetry, f: impl Fn(usize, usize) -> bool) {
        let nodes = self.nodes(f);
        self.bc(name, nodes, |target| BoundaryCondition::Symmetry { target, plane });
    }

    fn load(&mut self, name: &str, load: Load) {
        self.model.add_load(name, load).unwrap();
        self.model.set_load_active("Step", name, true).unwrap();
    }
}

/// The Scordelis–Lo roof's free-edge mid-span deflection on an n × n quarter mesh.
fn scordelis_lo(n: usize) -> f64 {
    // A quarter: x from the mid-span (x = 0) to the diaphragm (x = 25), φ from the crown to the
    // free edge at 40°.
    let mut c = Cylinder::new(n, 25.0, 25.0, 40f64.to_radians(), 0.25, Elastic::Isotropic { young: 4.32e8, poisson: 0.0 });
    c.diaphragm();
    c.symmetry("midspan", Symmetry::X, |i, _| i == 0);
    c.symmetry("crown", Symmetry::Y, |_, j| j == 0);
    // Self-weight 90 per area: 360 per volume over t = 0.25.
    c.load("weight", Load::BodyForce { set: "shell".into(), vector: [0.0, 0.0, -360.0] });
    let r = solve_static(&c.model, None, &Limits::default()).unwrap();
    assert!((r.reaction_total[2] + r.applied_total[2]).abs() < 1e-9 * r.applied_total[2].abs());
    let k = r.nodes.binary_search(&c.label(0, n)).unwrap();
    -r.displacement[k][2]
}

/// Scordelis–Lo reference deflection (MacNeal & Harder 1985).
const SCORDELIS_LO: f64 = 0.3024;

/// [benchmark: MacNeal & Harder 1985] Scordelis–Lo roof on a 16 × 16 quarter mesh: within 2 %.
#[test]
fn benchmark_scordelis_lo_roof_coarse() {
    let ratio = scordelis_lo(16) / SCORDELIS_LO;
    assert!((ratio - 1.0).abs() < 0.02, "{ratio}");
}

/// [benchmark: MacNeal & Harder 1985] Scordelis–Lo roof, free-edge mid-span deflection 0.3024:
/// monotone convergence to within 1 % at 64 × 64.
#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn benchmark_scordelis_lo_roof() {
    let results: Vec<(usize, f64, f64)> = [8, 16, 32, 64]
        .iter()
        .map(|&n| {
            let w = scordelis_lo(n);
            (n, w, w / SCORDELIS_LO)
        })
        .collect();
    eprintln!("Scordelis–Lo roof (reference {SCORDELIS_LO}): (cells per side, deflection, ratio) {results:.4?}");
    let (_, _, finest) = results[results.len() - 1];
    assert!((finest - 1.0).abs() < 0.01, "{results:?}");
    assert!(results.windows(2).all(|w| w[1].2 > w[0].2), "converges from below: {results:?}");
}

/// [benchmark: Belytschko et al. 1985] Pinched cylinder with diaphragms, deflection under the
/// load 1.8248 · 10⁻⁵.
#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn benchmark_pinched_cylinder() {
    const REFERENCE: f64 = 1.8248e-5;
    let mut results = Vec::new();
    for n in [8, 16, 32] {
        // An eighth: x from the loaded section (x = 0) to the diaphragm (x = 300), φ from the
        // loaded top to the side (z = 0).
        let mut c = Cylinder::new(n, 300.0, 300.0, std::f64::consts::FRAC_PI_2, 3.0, Elastic::Isotropic { young: 3e6, poisson: 0.3 });
        c.diaphragm();
        c.symmetry("middle", Symmetry::X, |i, _| i == 0);
        c.symmetry("top", Symmetry::Y, |_, j| j == 0);
        c.symmetry("side", Symmetry::Z, |_, j| j == n);
        let top = c.label(0, 0);
        c.model.add_set("pinch", Set { nodes: vec![top], ..Set::default() }).unwrap();
        // A quarter of the unit force acts on the eighth.
        c.load("pinch", Load::Force { target: Target::Set("pinch".into()), vector: [0.0, 0.0, -0.25], distribution: Distribution::Total });
        let r = solve_static(&c.model, None, &Limits::default()).unwrap();
        let k = r.nodes.binary_search(&top).unwrap();
        let w = -r.displacement[k][2];
        results.push((n, w, w / REFERENCE));
    }
    eprintln!("pinched cylinder (reference {REFERENCE:e}): (cells per side, deflection, ratio) {results:.4?}");
    let (_, _, finest) = results[results.len() - 1];
    assert!((finest - 1.0).abs() < 0.02, "{results:?}");
    assert!(results.windows(2).all(|w| w[1].2 > w[0].2), "converges from below: {results:?}");
}
