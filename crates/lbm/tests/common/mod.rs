//! Shared helpers for the lattice Boltzmann tests: models, analytic shapes generated in code, and
//! running a flow to a time.
#![allow(dead_code)]

pub mod curved_cht;

use std::collections::BTreeMap;

use fieldcraft_lbm::{Flow, Limits, Obstacle, Setup};
use fieldcraft_mesh::surface::{BoundarySurface, SurfaceOptions};
use fieldcraft_model::{BoxFaces, Collision, Fluid, FluidDomain, Material, Model, Procedure, Step};

/// A domain over [min, max] with walls on every face.
pub fn domain(min: [f64; 3], max: [f64; 3], spacing: f64) -> FluidDomain {
    FluidDomain {
        material: "Fluid".into(),
        min,
        max,
        spacing,
        faces: BoxFaces::default(),
        obstacles: vec![],
        surfaces: vec![],
        initial_velocity: [0.0; 3],
        acceleration: [0.0; 3],
        tolerance: None,
        thermal: None,
        porous: vec![],
        fans: vec![],
    }
}

/// A fluid step `Run` over `end_time` with a given time step.
pub fn step(end_time: f64, time_step: f64, magic: Option<f64>) -> Procedure {
    collide(end_time, time_step, magic, Collision::Trt)
}

/// [`step`] with a collision operator.
pub fn collide(end_time: f64, time_step: f64, magic: Option<f64>, collision: Collision) -> Procedure {
    Procedure::Fluid {
        domain: "Flow".into(),
        end_time: Some(end_time),
        steps: None,
        time_step: Some(time_step),
        frames: 2,
        magic,
        sigma: None,
        collision: Some(collision),
        turbulence: None,
        steady: None,
    }
}

/// A channel 2 × 1 × 1 at spacing 0.05 with a parabolic inlet, a pressure outlet, a slip face, a
/// ball (obstacle `Ball`, face `skin`) and a body force, run 0.5 s: every kind of neighbour a
/// cell can have. `sigma` and `turbulence` as the fluid step takes them.
pub fn ball_channel(
    collision: Collision,
    sigma: Option<fieldcraft_model::HrrSigma>,
    turbulence: Option<fieldcraft_model::Turbulence>,
) -> (Model, Obstacle) {
    let mut d = domain([0.0; 3], [2.0, 1.0, 1.0], 0.05);
    d.faces.x_min = fieldcraft_model::FaceCondition::Velocity { velocity: [0.1, 0.02, 0.0], profile: fieldcraft_model::Profile::Parabolic };
    d.faces.x_max = fieldcraft_model::FaceCondition::Pressure { pressure: 0.0, outflow: fieldcraft_model::Outflow::Uniform };
    d.faces.z_max = fieldcraft_model::FaceCondition::Slip {};
    d.obstacles = vec!["Ball".into()];
    d.acceleration = [0.0, 0.0, -0.01];
    let ball = obstacle("Ball", sphere([0.7, 0.5, 0.45], 0.2, 3), &["skin"]);
    let p = Procedure::Fluid {
        domain: "Flow".into(),
        end_time: Some(0.5),
        steps: None,
        time_step: None,
        frames: 2,
        magic: None,
        sigma,
        collision: Some(collision),
        turbulence,
        steady: None,
    };
    (model(1.0, 0.01, d, p), ball)
}

/// Both collision operators.
pub const COLLISIONS: [Collision; 2] = [Collision::Trt, Collision::Regularized];

/// A model with fluid `Fluid` (density ρ, kinematic viscosity ν), domain `Flow` and step `Run`.
pub fn model(density: f64, nu: f64, domain: FluidDomain, procedure: Procedure) -> Model {
    let mut m = Model::default();
    m.add_material("Fluid", Material { density: Some(density), fluid: Some(Fluid::Newtonian { viscosity: nu * density }), ..Material::default() })
        .unwrap();
    m.add_fluid_domain("Flow", domain).unwrap();
    m.add_step(Step::new("Run", procedure), None).unwrap();
    m
}

/// A fluid with heat: density ρ, kinematic viscosity ν, thermal diffusivity α (conductivity
/// ρ α, unit specific heat), so k / (ρ c) = α.
pub fn hot_model(density: f64, nu: f64, alpha: f64, domain: fieldcraft_model::FluidDomain, procedure: Procedure) -> Model {
    let mut m = Model::default();
    m.add_material("Fluid", hot_fluid(density, nu, alpha)).unwrap();
    m.add_fluid_domain("Flow", domain).unwrap();
    m.add_step(Step::new("Run", procedure), None).unwrap();
    m
}

/// A fluid with heat: density ρ, kinematic viscosity ν, thermal diffusivity α (conductivity ρ α,
/// unit specific heat).
pub fn hot_fluid(density: f64, nu: f64, alpha: f64) -> Material {
    Material {
        density: Some(density),
        fluid: Some(Fluid::Newtonian { viscosity: nu * density }),
        thermal: Some(fieldcraft_model::Thermal {
            conductivity: Some(fieldcraft_model::Directional::Isotropic(density * alpha)),
            specific_heat: Some(1.0),
            ..Default::default()
        }),
        ..Material::default()
    }
}

/// The setup and the flow at its initial state.
pub fn prepare(m: &Model, obstacles: &[Obstacle]) -> (Setup, Flow) {
    let setup = Setup::new(m, None, obstacles, &Limits::default()).unwrap();
    let flow = Flow::new(&setup, obstacles).unwrap();
    (setup, flow)
}

/// Run every step of the setup.
pub fn run_all(setup: &Setup, flow: &mut Flow) {
    for _ in 0..setup.steps {
        flow.step().unwrap();
    }
}

/// The convergence order of errors at resolutions growing by `ratio` each.
pub fn orders(errors: &[f64], ratio: f64) -> Vec<f64> {
    errors.windows(2).map(|w| (w[0] / w[1]).ln() / ratio.ln()).collect()
}

/// A closed surface from triangles, each turned to face `outward(centroid)`.
pub fn surface(
    positions: Vec<[f64; 3]>,
    mut triangles: Vec<[u32; 3]>,
    tags: Vec<u32>,
    outward: impl Fn([f64; 3], u32) -> [f64; 3],
) -> BoundarySurface {
    for (t, tag) in triangles.iter_mut().zip(&tags) {
        let [a, b, c] = t.map(|v| positions[v as usize]);
        let (u, v) = ([b[0] - a[0], b[1] - a[1], b[2] - a[2]], [c[0] - a[0], c[1] - a[1], c[2] - a[2]]);
        let n = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
        let centroid = [0, 1, 2].map(|k| (a[k] + b[k] + c[k]) / 3.0);
        let o = outward(centroid, *tag);
        if n[0] * o[0] + n[1] * o[1] + n[2] * o[2] < 0.0 {
            t.swap(1, 2);
        }
    }
    BoundarySurface::new(&positions, &triangles, &tags, &SurfaceOptions::default()).unwrap()
}

/// A sphere of radius `r` about `centre`: an icosahedron subdivided `levels` times, one face.
pub fn sphere(centre: [f64; 3], r: f64, levels: u32) -> BoundarySurface {
    let t = (1.0 + 5f64.sqrt()) / 2.0;
    let mut unit: Vec<[f64; 3]> = vec![
        [-1.0, t, 0.0],
        [1.0, t, 0.0],
        [-1.0, -t, 0.0],
        [1.0, -t, 0.0],
        [0.0, -1.0, t],
        [0.0, 1.0, t],
        [0.0, -1.0, -t],
        [0.0, 1.0, -t],
        [t, 0.0, -1.0],
        [t, 0.0, 1.0],
        [-t, 0.0, -1.0],
        [-t, 0.0, 1.0],
    ];
    let normalise = |p: [f64; 3]| {
        let l = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
        p.map(|c| c / l)
    };
    unit = unit.into_iter().map(normalise).collect();
    let mut tris: Vec<[u32; 3]> = vec![
        [0, 11, 5],
        [0, 5, 1],
        [0, 1, 7],
        [0, 7, 10],
        [0, 10, 11],
        [1, 5, 9],
        [5, 11, 4],
        [11, 10, 2],
        [10, 7, 6],
        [7, 1, 8],
        [3, 9, 4],
        [3, 4, 2],
        [3, 2, 6],
        [3, 6, 8],
        [3, 8, 9],
        [4, 9, 5],
        [2, 4, 11],
        [6, 2, 10],
        [8, 6, 7],
        [9, 8, 1],
    ];
    for _ in 0..levels {
        let mut mid: BTreeMap<(u32, u32), u32> = BTreeMap::new();
        let mut next = Vec::with_capacity(tris.len() * 4);
        for [a, b, c] in tris {
            let mut m = |p: u32, q: u32| {
                *mid.entry((p.min(q), p.max(q))).or_insert_with(|| {
                    let (x, y) = (unit[p as usize], unit[q as usize]);
                    unit.push(normalise([(x[0] + y[0]) / 2.0, (x[1] + y[1]) / 2.0, (x[2] + y[2]) / 2.0]));
                    (unit.len() - 1) as u32
                })
            };
            let (ab, bc, ca) = (m(a, b), m(b, c), m(c, a));
            next.extend([[a, ab, ca], [b, bc, ab], [c, ca, bc], [ab, bc, ca]]);
        }
        tris = next;
    }
    let positions: Vec<[f64; 3]> = unit.iter().map(|p| [centre[0] + r * p[0], centre[1] + r * p[1], centre[2] + r * p[2]]).collect();
    let tags = vec![0; tris.len()];
    surface(positions, tris, tags, |p, _| [p[0] - centre[0], p[1] - centre[1], p[2] - centre[2]])
}

/// A tube along x from `x0` to `x1` about the axis through (`y`, `z`): inner radius `r`, outer
/// radius `outer`, `m` segments. Faces: 0 inner wall, 1 outer wall, 2 and 3 the ends.
pub fn tube(x0: f64, x1: f64, axis: [f64; 2], r: f64, outer: f64, m: usize) -> BoundarySurface {
    let mut positions = Vec::new();
    for x in [x0, x1] {
        for radius in [r, outer] {
            for k in 0..m {
                let phi = 2.0 * std::f64::consts::PI * k as f64 / m as f64;
                positions.push([x, axis[0] + radius * phi.cos(), axis[1] + radius * phi.sin()]);
            }
        }
    }
    // Ring index: end e (0, 1), radius kind k (0 inner, 1 outer), segment s.
    let id = |e: usize, k: usize, s: usize| ((e * 2 + k) * m + s % m) as u32;
    let mut tris = Vec::new();
    let mut tags = Vec::new();
    for s in 0..m {
        for k in 0..2 {
            tris.push([id(0, k, s), id(1, k, s), id(1, k, s + 1)]);
            tris.push([id(0, k, s), id(1, k, s + 1), id(0, k, s + 1)]);
            tags.extend([k as u32; 2]);
        }
        for e in 0..2 {
            tris.push([id(e, 0, s), id(e, 1, s), id(e, 1, s + 1)]);
            tris.push([id(e, 0, s), id(e, 1, s + 1), id(e, 0, s + 1)]);
            tags.extend([2 + e as u32; 2]);
        }
    }
    surface(positions, tris, tags, move |p, tag| {
        let radial = [0.0, p[1] - axis[0], p[2] - axis[1]];
        match tag {
            0 => radial.map(|c| -c),
            1 => radial,
            2 => [-1.0, 0.0, 0.0],
            _ => [1.0, 0.0, 0.0],
        }
    })
}

/// A circular cylinder along z from `z0` to `z1` about (`x`, `y`), radius `r`, `m` segments.
/// Faces: 0 the side, 1 the bottom, 2 the top.
pub fn cylinder(centre: [f64; 2], r: f64, z0: f64, z1: f64, m: usize) -> BoundarySurface {
    let mut positions = Vec::new();
    for z in [z0, z1] {
        for k in 0..m {
            let phi = 2.0 * std::f64::consts::PI * k as f64 / m as f64;
            positions.push([centre[0] + r * phi.cos(), centre[1] + r * phi.sin(), z]);
        }
    }
    positions.push([centre[0], centre[1], z0]);
    positions.push([centre[0], centre[1], z1]);
    let id = |e: usize, s: usize| (e * m + s % m) as u32;
    let (bottom, top) = ((2 * m) as u32, (2 * m + 1) as u32);
    let mut tris = Vec::new();
    let mut tags = Vec::new();
    for s in 0..m {
        tris.push([id(0, s), id(0, s + 1), id(1, s + 1)]);
        tris.push([id(0, s), id(1, s + 1), id(1, s)]);
        tags.extend([0, 0]);
        tris.push([bottom, id(0, s + 1), id(0, s)]);
        tags.push(1);
        tris.push([top, id(1, s), id(1, s + 1)]);
        tags.push(2);
    }
    surface(positions, tris, tags, move |p, tag| match tag {
        0 => [p[0] - centre[0], p[1] - centre[1], 0.0],
        1 => [0.0, 0.0, -1.0],
        _ => [0.0, 0.0, 1.0],
    })
}

pub fn obstacle(body: &str, surface: BoundarySurface, faces: &[&str]) -> Obstacle {
    Obstacle { body: body.into(), surface, faces: faces.iter().map(|s| s.to_string()).collect() }
}

/// A temperature given in degrees Celsius, in kelvin (the model's unit): the lattice runs in
/// degrees Celsius, so these are the lattice's own temperatures.
pub fn kelvin(celsius: f64) -> f64 {
    fieldcraft_model::units::celsius_to_kelvin(celsius)
}

/// An axis-aligned box from `lo` to `hi`. Faces: 0 x_min, 1 x_max, 2 y_min, 3 y_max, 4 z_min,
/// 5 z_max.
pub fn cuboid(lo: [f64; 3], hi: [f64; 3]) -> BoundarySurface {
    let corner = |k: usize| [0, 1, 2].map(|a| if k >> a & 1 == 1 { hi[a] } else { lo[a] });
    let positions: Vec<[f64; 3]> = (0..8).map(corner).collect();
    let mut tris = Vec::new();
    let mut tags = Vec::new();
    for axis in 0..3 {
        let (t1, t2) = ((axis + 1) % 3, (axis + 2) % 3);
        for side in 0..2u32 {
            let v = |a: u32, b: u32| ((side as usize) << axis | (a as usize) << t1 | (b as usize) << t2) as u32;
            tris.push([v(0, 0), v(1, 0), v(1, 1)]);
            tris.push([v(0, 0), v(1, 1), v(0, 1)]);
            tags.extend([2 * axis as u32 + side; 2]);
        }
    }
    surface(positions, tris, tags, |_, tag| {
        let mut n = [0.0; 3];
        n[tag as usize / 2] = if tag % 2 == 1 { 1.0 } else { -1.0 };
        n
    })
}

/// The face names of a [`cuboid`].
pub const CUBOID_FACES: [&str; 6] = ["x_min", "x_max", "y_min", "y_max", "z_min", "z_max"];

/// Add a conducting material `name`: density ρ, conductivity k (per axis), specific heat c.
pub fn add_solid_material(m: &mut Model, name: &str, density: f64, k: [f64; 3], c: f64) {
    let conductivity =
        if k[0] == k[1] && k[1] == k[2] { fieldcraft_model::Directional::Isotropic(k[0]) } else { fieldcraft_model::Directional::Orthotropic(k) };
    m.add_material(
        name,
        Material {
            density: Some(density),
            thermal: Some(fieldcraft_model::Thermal { conductivity: Some(conductivity), specific_heat: Some(c), ..Default::default() }),
            ..Material::default()
        },
    )
    .unwrap();
}

/// A model with fluid `Fluid` (`[density, kinematic viscosity, diffusivity]`, see [`hot_fluid`]),
/// conducting material `Solid` (density, conductivity per axis, specific heat), domain `Flow`
/// and step `Run`.
pub fn cht_model(fluid: [f64; 3], solid: (f64, [f64; 3], f64), domain: FluidDomain, procedure: Procedure) -> Model {
    let [rho, nu, alpha] = fluid;
    let mut m = Model::default();
    m.add_material("Fluid", hot_fluid(rho, nu, alpha)).unwrap();
    add_solid_material(&mut m, "Solid", solid.0, solid.1, solid.2);
    m.add_fluid_domain("Flow", domain).unwrap();
    m.add_step(Step::new("Run", procedure), None).unwrap();
    m
}

/// Obstacle `body` conducting in material `Solid` with `power`.
pub fn conducting(body: &str, power: f64) -> fieldcraft_model::ConductingSolid {
    fieldcraft_model::ConductingSolid { body: body.into(), material: "Solid".into(), power, limit: None }
}

// ---- The hybrid collision's shear wave (crates/lbm/tests/hrr.rs, and the GPU's) ----

/// A complex number, for the shear wave's 2 × 2 map.
#[derive(Clone, Copy, Debug)]
pub struct Cx(f64, f64);

impl Cx {
    fn add(self, o: Cx) -> Cx {
        Cx(self.0 + o.0, self.1 + o.1)
    }
    fn sub(self, o: Cx) -> Cx {
        Cx(self.0 - o.0, self.1 - o.1)
    }
    fn mul(self, o: Cx) -> Cx {
        Cx(self.0 * o.0 - self.1 * o.1, self.0 * o.1 + self.1 * o.0)
    }
    fn scale(self, s: f64) -> Cx {
        Cx(self.0 * s, self.1 * s)
    }
    fn abs(self) -> f64 {
        self.0.hypot(self.1)
    }
    /// The principal square root.
    fn sqrt(self) -> Cx {
        let r = self.abs();
        let re = (0.5 * (r + self.0)).max(0.0).sqrt();
        let im = (0.5 * (r - self.0)).max(0.0).sqrt();
        Cx(re, if self.1 < 0.0 { -im } else { im })
    }
}

/// The eigenvalue of larger modulus of the complex matrix [[a, b], [c, d]].
fn dominant_eigenvalue([a, b, c, d]: [Cx; 4]) -> Cx {
    let half_trace = a.add(d).scale(0.5);
    let det = a.mul(d).sub(b.mul(c));
    let root = half_trace.mul(half_trace).sub(det).sqrt();
    let (l1, l2) = (half_trace.add(root), half_trace.sub(root));
    if l1.abs() >= l2.abs() { l1 } else { l2 }
}

/// The exact decay rate per step, −ln |λ|, of the shear wave `u_x ∝ sin(k y)` under the hybrid
/// collision with relaxation time `tau` (τ̄ = 1/ω) and weight `sigma`, `k` in radians per cell.
///
/// The wave is a single Fourier mode e^{iky} of a periodic lattice with no density
/// perturbation, so at linear order in its amplitude its state before a collision is two complex
/// amplitudes: the momentum `J` (x component) and the non-equilibrium flux `Π = Π_xy`. Products of
/// the velocity with the flux or with itself are second order, so the equilibrium's `ρ u u`, the
/// third-order Hermite terms (`u P̃`) and Π_xx, Π_yy drop out. The collision rebuilds
///
/// `f*_q = w_q [ρ + 3 c_q·j + 9 (1 − ω) c_qx c_qy P̃]`
///
/// (the D3Q19 second-order Hermite term for xy: `1/(2 c_s⁴) · 2 · c_x c_y P_xy = 9 c_x c_y P_xy`),
/// with `P̃ = σ Π + (1 − σ)(−2 c_s² τ̄ S^FD)` and `S^FD_xy = ½ ∂_y u_x` from centred differences:
/// `∂_y e^{iky} = (e^{ik} − e^{−ik})/2 = i sin k`, so `S^FD_xy = ½ i sin k J` and
///
/// `P̃ = σ Π − (1 − σ)(i/(3ω)) sin k · J`.
///
/// Streaming is `f_q(y, t+1) = f*_q(y − c_qy)`, a factor `e^{−ik c_qy}` per direction. Taking
/// the moments `J′ = Σ c_qx f_q` and `Π′ = Σ c_qx c_qy f_q` needs the D3Q19 directions with
/// `c_x = ±1`: (±1, 0, 0) of weight 1/18 and (±1, 0, ±1) of weight 1/36, together `Σ w = 2/9`
/// with `c_y = 0`; (±1, ±1, 0) of weight 1/36, `Σ w = 1/18` for each of `c_y = +1` and `−1`. So
///
/// - `Σ w c_x² e^{−ik c_y} = 2/9 + (1/18)(2 cos k) = (2 + cos k)/9`,
/// - `Σ w c_x² c_y e^{−ik c_y} = (1/18)(e^{−ik} − e^{ik}) = −(i/9) sin k`,
/// - `Σ w c_x² c_y² e^{−ik c_y} = (1/18)(2 cos k) = (1/9) cos k`,
///
/// and (the `ρ` term carries no x momentum or xy flux)
///
/// - `J′ = 3 J (2 + cos k)/9 + 9 (1 − ω) P̃ (−(i/9) sin k) = (2 + cos k)/3 · J − i (1 − ω) sin k · P̃`,
/// - `Π′ = 3 J (−(i/9) sin k) + 9 (1 − ω) P̃ (1/9) cos k = −(i/3) sin k · J + (1 − ω) cos k · P̃`.
///
/// Substituting P̃ gives the 2 × 2 map of (J, Π); after the fast mode has died out the wave's
/// amplitude is multiplied by its dominant eigenvalue λ every step. With σ = 1 this is the
/// regularised collision's map.
pub fn exact_rate(k: f64, tau: f64, sigma: f64) -> f64 {
    let omega = 1.0 / tau;
    let (s, c) = (k.sin(), k.cos());
    // P̃ = σ Π + β J
    let beta = Cx(0.0, -(1.0 - sigma) * s / (3.0 * omega));
    let gamma = Cx(0.0, -(1.0 - omega) * s); // J′ += γ P̃
    let delta = Cx((1.0 - omega) * c, 0.0); // Π′ += δ P̃
    let a = Cx((2.0 + c) / 3.0, 0.0).add(gamma.mul(beta));
    let b = gamma.scale(sigma);
    let cc = Cx(0.0, -s / 3.0).add(delta.mul(beta));
    let d = delta.scale(sigma);
    -dominant_eigenvalue([a, b, cc, d]).abs().ln()
}
