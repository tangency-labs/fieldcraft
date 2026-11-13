//! A fluid step written as an OpenFOAM case: the dictionaries of a `blockMesh` + `snappyHexMesh`
//! mesh of the fluid domain around its obstacles, the initial and boundary fields, and the
//! controls of an incompressible transient run (`pimpleFoam`), as a map of case-relative paths
//! to text. The formats come from OpenFOAM's public user guide (case structure, dictionary
//! syntax, boundary conditions); nothing of OpenFOAM's source is read or copied.
//!
//! - The domain box is the background block (cells of `spacing × cell_factor`); its six faces
//!   are the patches `x_min` … `z_max`.
//! - Each obstacle body is a multi-solid ASCII STL in `constant/triSurface`, one solid per named
//!   face, so `snappyHexMesh` makes a patch per face (`<body>_<face>`, sanitised) and the
//!   obstacle's surface is refined back to the lattice spacing.
//! - Face conditions: velocity → fixed velocity; pressure → fixed kinematic pressure (p/ρ) with
//!   an inflow/outflow velocity; slip → `symmetryPlane` (box faces are planar); wall → no-slip;
//!   moving wall → fixed velocity. Obstacle faces are no-slip unless the domain names them.
//! - The step's length and frames give `endTime` and the write interval; the time step adapts
//!   to a Courant limit (the lattice's time step is the LBM's own). A Smagorinsky model becomes
//!   OpenFOAM's (its `Ck` from the constant, `Ce` 1.048), with Van Driest damping when the step
//!   has it.
//! - With a [`Gas`] the case is compressible (`rhoCentralFoam`, laminar, ideal gas): pressure
//!   faces above the ambient hold their total pressure and temperature (inlets), the others are
//!   non-reflecting outlets at the ambient; walls are adiabatic.
//! - Fields are written as ASCII, and `Allrun` writes the cell centres (`0/C`), so the results
//!   can be read back ([`read_field`]).
//!
//! Everything is SI, as the model stores it and OpenFOAM reads it: no unit conversion happens
//! here. What the case cannot express yet is refused by name ([`crate::Error::Unsupported`]).

use std::collections::BTreeMap;
use std::fmt::Write as _;

use fieldcraft_model::{BoxFace, FaceCondition, Fluid, FluidDomain, Model, Procedure, Profile, Step, SurfaceCondition, Turbulence};

use crate::{Error, Result};

const FORMAT: &str = "OpenFOAM";

/// The most triangles the obstacle surfaces may have together (the STL text grows with them).
pub const MAX_TRIANGLES: usize = 5_000_000;

/// One obstacle body's surface, in metres: its triangles, the face each triangle belongs to
/// (an index into `faces`) and the faces' persistent names.
#[derive(Debug, Clone, Copy)]
pub struct ObstacleSurface<'a> {
    pub body: &'a str,
    pub positions: &'a [[f64; 3]],
    pub triangles: &'a [[u32; 3]],
    pub tri_face: &'a [u32],
    pub faces: &'a [String],
}

/// A box refined to `spacing` (metres; the lattice spacing when left out, rounded to the
/// background cell over a power of two).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RefineBox {
    pub min: [f64; 3],
    pub max: [f64; 3],
    pub spacing: Option<f64>,
}

/// An ideal gas for a compressible case (SI).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Gas {
    /// The absolute pressure the model's gauge pressures are above (Pa).
    pub ambient: f64,
    /// The total temperature at inlets and the initial temperature (K).
    pub temperature: f64,
    /// kg/mol.
    pub molar_mass: f64,
    /// At constant pressure (J/(kg·K)).
    pub specific_heat: f64,
    pub prandtl: f64,
}

/// The universal gas constant (J/(mol·K)).
pub const GAS_CONSTANT: f64 = 8.314_462_618;

impl Gas {
    /// Dry air at 1 atm and 20 °C.
    pub const AIR: Gas = Gas { ambient: 101_325.0, temperature: 293.15, molar_mass: 0.028_964_7, specific_heat: 1004.5, prandtl: 0.71 };

    /// The ratio of specific heats.
    pub fn gamma(&self) -> f64 {
        self.specific_heat / (self.specific_heat - GAS_CONSTANT / self.molar_mass)
    }

    fn check(&self) -> Result<()> {
        let positive = [self.ambient, self.temperature, self.molar_mass, self.specific_heat, self.prandtl].iter().all(|v| v.is_finite() && *v > 0.0);
        let g = self.gamma();
        if !positive || !(g.is_finite() && g > 1.0 && g < 2.0) {
            return Err(unsupported(
                "a gas whose ambient pressure, temperature, molar mass, specific heat and Prandtl number are not positive, or whose ratio of specific heats is not between 1 and 2",
            ));
        }
        Ok(())
    }
}

/// How the case is meshed and run.
#[derive(Debug, Clone, PartialEq)]
pub struct CaseOptions {
    /// The background cells are the domain's spacing times this (1, 2, 4 or 8); obstacle
    /// surfaces and refinement boxes are refined back to the spacing.
    pub cell_factor: u32,
    pub refine_boxes: Vec<RefineBox>,
    /// A point in the fluid (metres) that tells `snappyHexMesh` which side to keep; found from
    /// the open faces when left out.
    pub location: Option<[f64; 3]>,
    /// MPI ranks `Allrun` runs the solver on (1: serial).
    pub processes: u32,
    /// The Courant number the time step adapts to.
    pub max_courant: f64,
    /// Compressible flow of this gas (`rhoCentralFoam`) instead of incompressible (`pimpleFoam`).
    pub gas: Option<Gas>,
}

impl Default for CaseOptions {
    fn default() -> Self {
        CaseOptions { cell_factor: 2, refine_boxes: Vec::new(), location: None, processes: 1, max_courant: 0.8, gas: None }
    }
}

/// A written case.
#[derive(Debug, Clone, PartialEq)]
pub struct Case {
    /// Case-relative path (with `/`) → text, in path order.
    pub files: BTreeMap<String, String>,
    /// OpenFOAM patch → what it is in the model (`x_min`, or `Body1/F2:top`).
    pub patches: BTreeMap<String, String>,
    /// The background block's cells along x, y and z.
    pub background_cells: [u64; 3],
    /// The point in the fluid `snappyHexMesh` keeps.
    pub location: [f64; 3],
    /// The solver application.
    pub solver: &'static str,
    /// What the case leaves out or does differently from the lattice Boltzmann step.
    pub notes: Vec<String>,
}

fn unsupported(what: impl Into<String>) -> Error {
    Error::Unsupported { format: FORMAT, what: what.into() }
}

/// The fluid step `step` (the model's only fluid step when `None`) as an OpenFOAM case.
pub fn write_case(model: &Model, step: Option<&str>, obstacles: &[ObstacleSurface<'_>], options: &CaseOptions) -> Result<Case> {
    let step = fluid_step(model, step)?;
    let Procedure::Fluid { domain: domain_name, end_time, steps, time_step, frames, turbulence, steady, .. } = &step.procedure else {
        return Err(unsupported(format!("step `{}` is not a fluid step", step.name)));
    };
    let domain =
        model.fluid_domains().get(domain_name).ok_or_else(|| unsupported(format!("step `{}` names no fluid domain `{domain_name}`", step.name)))?;
    let mut notes = Vec::new();
    refuse_unsupported(domain, &mut notes)?;
    let material = model
        .materials()
        .get(&domain.material)
        .ok_or_else(|| unsupported(format!("fluid domain `{domain_name}` names no material `{}`", domain.material)))?;
    let density = material.density.filter(|d| d.is_finite() && *d > 0.0).ok_or_else(|| unsupported("a fluid material without a density"))?;
    let viscosity = match material.fluid {
        Some(Fluid::Newtonian { viscosity }) if viscosity.is_finite() && viscosity >= 0.0 => viscosity,
        Some(Fluid::Newtonian { .. }) => return Err(unsupported("a fluid viscosity that is not a finite, non-negative number")),
        Some(_) => return Err(unsupported("shear-dependent viscosity (only Newtonian fluids so far)")),
        None => return Err(unsupported(format!("material `{}` is not a fluid", domain.material))),
    };
    let duration = match (end_time, steps, time_step) {
        (Some(t), _, _) => *t,
        (None, Some(n), Some(dt)) => *n as f64 * dt,
        _ => return Err(unsupported("a step length in lattice steps without a time step (give end_time)")),
    };
    if !(duration.is_finite() && duration > 0.0) {
        return Err(unsupported("a step length that is not a positive time"));
    }
    if time_step.is_some() {
        notes.push("the lattice time step is not used: OpenFOAM's adapts to the Courant limit".into());
    }
    if steady.is_some() {
        notes.push("the steady-state stop rule is not mapped: the run goes to the end time".into());
    }
    if !matches!(options.cell_factor, 1 | 2 | 4 | 8) {
        return Err(unsupported("a cell factor other than 1, 2, 4 or 8"));
    }
    if !(options.max_courant.is_finite() && options.max_courant > 0.0 && options.max_courant <= 5.0) {
        return Err(unsupported("a Courant limit outside (0, 5]"));
    }
    let processes = options.processes.clamp(1, 4096);

    let mut files = BTreeMap::new();
    let mut patches = BTreeMap::new();
    for face in BoxFace::ALL {
        patches.insert(face.name().to_string(), face.name().to_string());
    }

    // Obstacles: one STL per body, a solid per named face.
    let total: usize = obstacles.iter().map(|o| o.triangles.len()).sum();
    if total > MAX_TRIANGLES {
        return Err(Error::Limit(format!("{total} obstacle triangles (at most {MAX_TRIANGLES})")));
    }
    let mut surfaces = Vec::new(); // (surface name, file, [(region, patch, condition)])
    let mut used = BTreeMap::<String, ()>::new();
    for o in obstacles {
        if !domain.obstacles.iter().any(|b| b == o.body) {
            continue;
        }
        let name = unique(&sanitize(o.body), &mut used);
        let (stl, regions) = stl(o)?;
        let mut rs = Vec::new();
        for (region, face) in regions {
            let patch = format!("{name}_{region}");
            let condition = surface_condition(domain, o.body, &face);
            patches.insert(patch.clone(), format!("{}/{face}", o.body));
            rs.push((region, patch, condition));
        }
        files.insert(format!("constant/triSurface/{name}.stl"), stl);
        surfaces.push((name, rs));
    }
    for body in &domain.obstacles {
        if !obstacles.iter().any(|o| o.body == body) {
            return Err(unsupported(format!("obstacle `{body}` has no surface")));
        }
    }

    let spacing = domain.spacing;
    let background = spacing * f64::from(options.cell_factor);
    let mut cells = [0u64; 3];
    for (axis, c) in cells.iter_mut().enumerate() {
        let extent = domain.max.get(axis).copied().unwrap_or(0.0) - domain.min.get(axis).copied().unwrap_or(0.0);
        let n = (extent / background).round();
        if !(n.is_finite() && (1.0..=1.0e5).contains(&n)) {
            return Err(unsupported("a domain box whose extent is not 1 to 100000 background cells"));
        }
        *c = n as u64;
    }
    let level = options.cell_factor.trailing_zeros();
    let mut boxes = Vec::with_capacity(options.refine_boxes.len());
    for b in &options.refine_boxes {
        let box_level = match b.spacing {
            None => level,
            Some(h) if h.is_finite() && h > 0.0 => {
                let l = (background / h).log2().round();
                if !(0.0..=10.0).contains(&l) {
                    return Err(unsupported("a refinement box spacing that is not the background cell over 1 to 1024"));
                }
                l as u32
            }
            Some(_) => return Err(unsupported("a refinement box spacing that is not positive")),
        };
        if b.min.iter().chain(&b.max).any(|v| !v.is_finite()) || (0..3).any(|i| b.min[i] >= b.max[i]) {
            return Err(unsupported("a refinement box whose min is not below its max"));
        }
        boxes.push((*b, box_level));
    }
    let location = match options.location {
        Some(p) => p,
        None => find_location(domain, background, obstacles)?,
    };

    files.insert("system/blockMeshDict".into(), block_mesh(domain, cells)?);
    files.insert("system/surfaceFeatureExtractDict".into(), feature_extract(&surfaces));
    files.insert("system/snappyHexMeshDict".into(), snappy(&surfaces, level, &boxes, location));
    let frames = (*frames).max(1);
    let solver = match options.gas {
        None => {
            for (field, text) in fields(domain, density, &surfaces) {
                files.insert(format!("0/{field}"), text);
            }
            files.insert(
                "constant/transportProperties".into(),
                dictionary("transportProperties", &format!("transportModel Newtonian;\nnu {};\n", num(viscosity / density))),
            );
            files.insert("constant/turbulenceProperties".into(), turbulence_properties(turbulence.as_ref(), &mut notes)?);
            files.insert("system/fvSchemes".into(), dictionary("fvSchemes", FV_SCHEMES));
            files.insert("system/fvSolution".into(), dictionary("fvSolution", FV_SOLUTION));
            "pimpleFoam"
        }
        Some(gas) => {
            gas.check()?;
            notes.push(format!("compressible: ideal gas, γ {:.4}; the material's density is not used", gas.gamma()));
            if turbulence.is_some() {
                notes.push("the subgrid model is not mapped in a compressible case: the flow is laminar".into());
            }
            for (field, text) in gas_fields(domain, &gas, &surfaces) {
                files.insert(format!("0/{field}"), text);
            }
            files.insert("constant/thermophysicalProperties".into(), thermophysical(&gas, viscosity));
            files.insert("constant/turbulenceProperties".into(), dictionary("turbulenceProperties", "simulationType laminar;\n"));
            files.insert("system/fvSchemes".into(), dictionary("fvSchemes", GAS_SCHEMES));
            files.insert("system/fvSolution".into(), dictionary("fvSolution", GAS_SOLUTION));
            "rhoCentralFoam"
        }
    };
    files.insert("system/controlDict".into(), control(solver, duration, frames, options.max_courant));
    files.insert("system/decomposeParDict".into(), dictionary("decomposeParDict", &format!("numberOfSubdomains {processes};\nmethod scotch;\n")));
    files.insert("Allrun".into(), allrun(&step.name, solver, processes));
    files.insert("case.foam".into(), String::new());

    Ok(Case { files, patches, background_cells: cells, location, solver, notes })
}

fn fluid_step<'m>(model: &'m Model, name: Option<&str>) -> Result<&'m Step> {
    let mut fluid = model.steps().iter().filter(|s| matches!(s.procedure, Procedure::Fluid { .. }));
    match name {
        Some(n) => model.steps().iter().find(|s| s.name == n).ok_or_else(|| unsupported(format!("no step `{n}`"))),
        None => match (fluid.next(), fluid.next()) {
            (Some(s), None) => Ok(s),
            (None, _) => Err(unsupported("a model without a fluid step")),
            (Some(_), Some(_)) => Err(unsupported("several fluid steps without naming one")),
        },
    }
}

fn refuse_unsupported(d: &FluidDomain, notes: &mut Vec<String>) -> Result<()> {
    let refused = [
        (d.thermal.is_some(), "heat in the fluid"),
        (!d.porous.is_empty(), "porous zones"),
        (!d.fans.is_empty(), "fans"),
        (!d.baffles.is_empty(), "thin walls (baffles)"),
        (d.free_surface.is_some(), "a free surface"),
        (d.acceleration.iter().any(|a| *a != 0.0), "a body acceleration"),
    ];
    if let Some((_, what)) = refused.iter().find(|(on, _)| *on) {
        return Err(unsupported(format!("{what} (not mapped yet)")));
    }
    for face in BoxFace::ALL {
        match d.faces.get(face) {
            FaceCondition::Periodic {} => return Err(unsupported(format!("a periodic face ({}; not mapped yet)", face.name()))),
            FaceCondition::Velocity { profile: Profile::Parabolic, .. } => {
                return Err(unsupported(format!("a parabolic inlet profile ({}; not mapped yet)", face.name())));
            }
            FaceCondition::Pressure { outflow, .. } if !outflow.is_uniform() => {
                notes.push(format!("{}: the convective outflow is written as an inflow/outflow pressure face", face.name()));
            }
            _ => {}
        }
    }
    Ok(())
}

fn surface_condition(d: &FluidDomain, body: &str, face: &str) -> SurfaceCondition {
    d.surfaces
        .iter()
        .find(|s| {
            s.face.body == body && s.face.entity.name().is_some_and(|n| n == face || face.strip_prefix(n).is_some_and(|rest| rest.starts_with('#')))
        })
        .map_or(SurfaceCondition::Wall {}, |s| s.condition)
}

/// A name OpenFOAM takes as a word and a patch: letters, digits and `_`, starting with a letter
/// (`+` and `-` of face names such as `F1:+x` become `p` and `m`: `F1_px`).
fn sanitize(name: &str) -> String {
    let mut out = String::with_capacity(name.len());
    for c in name.chars() {
        match c {
            '+' => out.push('p'),
            '-' => out.push('m'),
            c if c.is_ascii_alphanumeric() || c == '_' => out.push(c),
            _ => out.push('_'),
        }
    }
    if !out.starts_with(|c: char| c.is_ascii_alphabetic()) {
        out.insert_str(0, "s_");
    }
    out
}

fn unique(base: &str, used: &mut BTreeMap<String, ()>) -> String {
    let mut name = base.to_string();
    let mut k = 2;
    while used.contains_key(&name) {
        name = format!("{base}_{k}");
        k += 1;
    }
    used.insert(name.clone(), ());
    name
}

/// The body as ASCII STL, one solid per face that has triangles (in face order), and the
/// solids' (region, face name) pairs.
fn stl(o: &ObstacleSurface<'_>) -> Result<(String, Vec<(String, String)>)> {
    if o.tri_face.len() != o.triangles.len() {
        return Err(unsupported(format!("obstacle `{}`: a face index per triangle", o.body)));
    }
    let mut by_face = BTreeMap::<u32, Vec<usize>>::new();
    for (t, f) in o.tri_face.iter().enumerate() {
        by_face.entry(*f).or_default().push(t);
    }
    let mut used = BTreeMap::new();
    let mut text = String::new();
    let mut regions = Vec::new();
    for (f, tris) in by_face {
        let face = o.faces.get(f as usize).cloned().unwrap_or_else(|| format!("face{f}"));
        let region = unique(&sanitize(&face), &mut used);
        let _ = writeln!(text, "solid {region}");
        for t in tris {
            let tri = o.triangles.get(t).ok_or_else(|| unsupported("a triangle index out of range"))?;
            let p = |i: u32| {
                o.positions.get(i as usize).copied().ok_or_else(|| unsupported(format!("obstacle `{}`: a vertex index out of range", o.body)))
            };
            let (a, b, c) = (p(tri[0])?, p(tri[1])?, p(tri[2])?);
            if a.iter().chain(&b).chain(&c).any(|x| !x.is_finite()) {
                return Err(unsupported(format!("obstacle `{}`: a vertex that is not finite", o.body)));
            }
            let n = normal(a, b, c);
            let _ = writeln!(text, "  facet normal {} {} {}\n    outer loop", num(n[0]), num(n[1]), num(n[2]));
            for v in [a, b, c] {
                let _ = writeln!(text, "      vertex {} {} {}", num(v[0]), num(v[1]), num(v[2]));
            }
            let _ = writeln!(text, "    endloop\n  endfacet");
        }
        let _ = writeln!(text, "endsolid {region}");
        regions.push((region, face));
    }
    Ok((text, regions))
}

fn normal(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> [f64; 3] {
    let u = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let v = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
    let n = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
    let l = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
    if l > 0.0 { [n[0] / l, n[1] / l, n[2] / l] } else { [0.0; 3] }
}

/// A point in the fluid: just inside the middle of an open box face (pressure faces first),
/// else the box's centre, that is in no obstacle (ray parity).
fn find_location(d: &FluidDomain, background: f64, obstacles: &[ObstacleSurface<'_>]) -> Result<[f64; 3]> {
    let inset = 1.37 * background;
    let centre = [0, 1, 2].map(|i| d.min[i] + 0.5137 * (d.max[i] - d.min[i]));
    let mut candidates = Vec::new();
    let open = |c: &FaceCondition| matches!(c, FaceCondition::Pressure { .. } | FaceCondition::Velocity { .. });
    for pressure_first in [true, false] {
        for face in BoxFace::ALL {
            let c = d.faces.get(face);
            if open(c) && matches!(c, FaceCondition::Pressure { .. }) == pressure_first {
                let mut p = centre;
                let a = face.axis();
                p[a] = if face.is_max() { d.max[a] - inset } else { d.min[a] + inset };
                candidates.push(p);
            }
        }
    }
    candidates.push(centre);
    candidates
        .into_iter()
        .find(|p| !obstacles.iter().any(|o| inside(*p, o)))
        .ok_or_else(|| unsupported("finding a point in the fluid (give `location`)"))
}

/// Whether `p` is inside the closed surface: crossings of a ray in a fixed oblique direction.
fn inside(p: [f64; 3], o: &ObstacleSurface<'_>) -> bool {
    let dir = [0.8017837257372732, 0.5345224838248488, 0.2672612419124244];
    let mut crossings = 0u32;
    for tri in o.triangles {
        let (Some(a), Some(b), Some(c)) = (o.positions.get(tri[0] as usize), o.positions.get(tri[1] as usize), o.positions.get(tri[2] as usize))
        else {
            continue;
        };
        if ray_hits(p, dir, *a, *b, *c) {
            crossings += 1;
        }
    }
    crossings % 2 == 1
}

fn ray_hits(o: [f64; 3], d: [f64; 3], a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> bool {
    let sub = |x: [f64; 3], y: [f64; 3]| [x[0] - y[0], x[1] - y[1], x[2] - y[2]];
    let cross = |x: [f64; 3], y: [f64; 3]| [x[1] * y[2] - x[2] * y[1], x[2] * y[0] - x[0] * y[2], x[0] * y[1] - x[1] * y[0]];
    let dot = |x: [f64; 3], y: [f64; 3]| x[0] * y[0] + x[1] * y[1] + x[2] * y[2];
    let (e1, e2) = (sub(b, a), sub(c, a));
    let h = cross(d, e2);
    let det = dot(e1, h);
    if det.abs() < 1e-30 {
        return false;
    }
    let s = sub(o, a);
    let u = dot(s, h) / det;
    if !(0.0..=1.0).contains(&u) {
        return false;
    }
    let q = cross(s, e1);
    let v = dot(d, q) / det;
    if v < 0.0 || u + v > 1.0 {
        return false;
    }
    dot(e2, q) / det > 0.0
}

/// Shortest round-trip text of a number (the same on every platform).
fn num(x: f64) -> String {
    if x == 0.0 { "0".into() } else { format!("{x:e}") }
}

fn vector(v: [f64; 3]) -> String {
    format!("({} {} {})", num(v[0]), num(v[1]), num(v[2]))
}

fn dictionary(object: &str, body: &str) -> String {
    format!("FoamFile\n{{\n    version 2.0;\n    format ascii;\n    class dictionary;\n    object {object};\n}}\n// Written by FieldCraft.\n\n{body}")
}

fn field_file(class: &str, object: &str, body: &str) -> String {
    format!("FoamFile\n{{\n    version 2.0;\n    format ascii;\n    class {class};\n    object {object};\n}}\n// Written by FieldCraft.\n\n{body}")
}

/// Hex vertex order (OpenFOAM's block convention) and each box face's vertices, outward.
const BOX_FACES: [(BoxFace, [usize; 4]); 6] = [
    (BoxFace::XMin, [0, 4, 7, 3]),
    (BoxFace::XMax, [2, 6, 5, 1]),
    (BoxFace::YMin, [1, 5, 4, 0]),
    (BoxFace::YMax, [3, 7, 6, 2]),
    (BoxFace::ZMin, [0, 3, 2, 1]),
    (BoxFace::ZMax, [4, 5, 6, 7]),
];

fn patch_type(c: &FaceCondition) -> &'static str {
    match c {
        FaceCondition::Slip {} => "symmetryPlane",
        FaceCondition::Wall {} | FaceCondition::MovingWall { .. } => "wall",
        _ => "patch",
    }
}

fn block_mesh(d: &FluidDomain, n: [u64; 3]) -> Result<String> {
    let (lo, hi) = (d.min, d.max);
    let corners = [
        [lo[0], lo[1], lo[2]],
        [hi[0], lo[1], lo[2]],
        [hi[0], hi[1], lo[2]],
        [lo[0], hi[1], lo[2]],
        [lo[0], lo[1], hi[2]],
        [hi[0], lo[1], hi[2]],
        [hi[0], hi[1], hi[2]],
        [lo[0], hi[1], hi[2]],
    ];
    let mut s = String::from("scale 1;\n\nvertices\n(\n");
    for c in corners {
        let _ = writeln!(s, "    {}", vector(c));
    }
    let _ = writeln!(s, ");\n\nblocks\n(\n    hex (0 1 2 3 4 5 6 7) ({} {} {}) simpleGrading (1 1 1)\n);\n\nboundary\n(", n[0], n[1], n[2]);
    for (face, v) in BOX_FACES {
        let _ = writeln!(
            s,
            "    {}\n    {{\n        type {};\n        faces (({} {} {} {}));\n    }}",
            face.name(),
            patch_type(d.faces.get(face)),
            v[0],
            v[1],
            v[2],
            v[3]
        );
    }
    s.push_str(");\n");
    Ok(dictionary("blockMeshDict", &s))
}

type Surfaces = [(String, Vec<(String, String, SurfaceCondition)>)];

fn feature_extract(surfaces: &Surfaces) -> String {
    let mut s = String::new();
    for (name, _) in surfaces {
        let _ = writeln!(s, "{name}.stl\n{{\n    extractionMethod extractFromSurface;\n    includedAngle 150;\n    writeObj no;\n}}\n");
    }
    dictionary("surfaceFeatureExtractDict", &s)
}

fn snappy(surfaces: &Surfaces, level: u32, boxes: &[(RefineBox, u32)], location: [f64; 3]) -> String {
    let mut geometry = String::new();
    let mut features = String::new();
    let mut refinement = String::new();
    for (name, _) in surfaces {
        let _ = writeln!(geometry, "    {name}.stl {{ type triSurfaceMesh; name {name}; }}");
        let _ = writeln!(features, "        {{ file \"{name}.eMesh\"; level {level}; }}");
        let _ = writeln!(refinement, "        {name} {{ level ({level} {level}); patchInfo {{ type wall; }} }}");
    }
    let mut regions = String::new();
    for (i, (b, box_level)) in boxes.iter().enumerate() {
        let _ = writeln!(geometry, "    refine{i} {{ type searchableBox; min {}; max {}; }}", vector(b.min), vector(b.max));
        let _ = writeln!(regions, "        refine{i} {{ mode inside; levels ((1e15 {box_level})); }}");
    }
    let body = format!(
        "castellatedMesh true;\nsnap true;\naddLayers false;\n\n\
geometry\n{{\n{geometry}}}\n\n\
castellatedMeshControls\n{{\n    maxLocalCells 2000000;\n    maxGlobalCells 20000000;\n    minRefinementCells 0;\n    maxLoadUnbalance 0.1;\n    nCellsBetweenLevels 3;\n\
    features\n    (\n{features}    );\n    refinementSurfaces\n    {{\n{refinement}    }}\n    resolveFeatureAngle 30;\n\
    refinementRegions\n    {{\n{regions}    }}\n    locationInMesh {};\n    allowFreeStandingZoneFaces true;\n}}\n\n\
snapControls\n{{\n    nSmoothPatch 3;\n    tolerance 2;\n    nSolveIter 50;\n    nRelaxIter 5;\n    nFeatureSnapIter 10;\n    implicitFeatureSnap false;\n    explicitFeatureSnap true;\n    multiRegionFeatureSnap false;\n}}\n\n\
addLayersControls\n{{\n    relativeSizes true;\n    layers {{}}\n    expansionRatio 1.2;\n    finalLayerThickness 0.4;\n    minThickness 0.1;\n    nGrow 0;\n    featureAngle 60;\n    nRelaxIter 3;\n    nSmoothSurfaceNormals 1;\n    nSmoothNormals 3;\n    nSmoothThickness 10;\n    maxFaceThicknessRatio 0.5;\n    maxThicknessToMedialRatio 0.3;\n    minMedialAxisAngle 90;\n    nBufferCellsNoExtrude 0;\n    nLayerIter 50;\n}}\n\n\
meshQualityControls\n{{\n    #includeEtc \"caseDicts/meshQualityDict\"\n    nSmoothScale 4;\n    errorReduction 0.75;\n    relaxed {{ maxNonOrtho 75; }}\n}}\n\nmergeTolerance 1e-6;\n",
        vector(location)
    );
    dictionary("snappyHexMeshDict", &body)
}

/// `0/U`, `0/p` and `0/nut`: each box face and each obstacle face (by patch name).
fn fields(d: &FluidDomain, density: f64, surfaces: &Surfaces) -> [(&'static str, String); 3] {
    let mut u = String::new();
    let mut p = String::new();
    let mut nut = String::new();
    for face in BoxFace::ALL {
        let (fu, fp, fn_) = match *d.faces.get(face) {
            FaceCondition::Velocity { velocity, .. } => {
                (format!("type fixedValue; value uniform {};", vector(velocity)), "type zeroGradient;".to_string(), "type zeroGradient;")
            }
            FaceCondition::Pressure { pressure, .. } => (
                "type pressureInletOutletVelocity; value uniform (0 0 0);".to_string(),
                format!("type fixedValue; value uniform {};", num(pressure / density)),
                "type zeroGradient;",
            ),
            FaceCondition::Slip {} => ("type symmetryPlane;".into(), "type symmetryPlane;".into(), "type symmetryPlane;"),
            FaceCondition::Wall {} => ("type noSlip;".into(), "type zeroGradient;".into(), "type zeroGradient;"),
            FaceCondition::MovingWall { velocity } => {
                (format!("type fixedValue; value uniform {};", vector(velocity)), "type zeroGradient;".into(), "type zeroGradient;")
            }
            // Refused before the fields are written.
            FaceCondition::Periodic {} => ("type cyclic;".into(), "type cyclic;".into(), "type cyclic;"),
        };
        let n = face.name();
        let _ = writeln!(u, "    {n} {{ {fu} }}");
        let _ = writeln!(p, "    {n} {{ {fp} }}");
        let _ = writeln!(nut, "    {n} {{ {fn_} }}");
    }
    for (_, regions) in surfaces {
        for (_, patch, condition) in regions {
            let fu = match condition {
                SurfaceCondition::Wall {} => "type noSlip;".to_string(),
                SurfaceCondition::MovingWall { velocity } => format!("type fixedValue; value uniform {};", vector(*velocity)),
            };
            let _ = writeln!(u, "    {patch} {{ {fu} }}");
            let _ = writeln!(p, "    {patch} {{ type zeroGradient; }}");
            let _ = writeln!(nut, "    {patch} {{ type zeroGradient; }}");
        }
    }
    let file = |class: &str, object: &str, dim: &str, internal: String, b: &str| {
        field_file(class, object, &format!("dimensions {dim};\n\ninternalField uniform {internal};\n\nboundaryField\n{{\n{b}}}\n"))
    };
    [
        ("U", file("volVectorField", "U", "[0 1 -1 0 0 0 0]", vector(d.initial_velocity), &u)),
        ("p", file("volScalarField", "p", "[0 2 -2 0 0 0 0]", "0".into(), &p)),
        ("nut", file("volScalarField", "nut", "[0 2 -1 0 0 0 0]", "0".into(), &nut)),
    ]
}

/// `0/U`, `0/p` (absolute) and `0/T` of a compressible case.
fn gas_fields(d: &FluidDomain, gas: &Gas, surfaces: &Surfaces) -> [(&'static str, String); 3] {
    let g = num(gas.gamma());
    let t0 = num(gas.temperature);
    let ambient = num(gas.ambient);
    let far = (0..3).map(|i| d.max[i] - d.min[i]).fold(0.0, f64::max);
    let mut u = String::new();
    let mut p = String::new();
    let mut t = String::new();
    for face in BoxFace::ALL {
        let (fu, fp, ft) = match *d.faces.get(face) {
            FaceCondition::Velocity { velocity, .. } => (
                format!("type fixedValue; value uniform {};", vector(velocity)),
                "type zeroGradient;".to_string(),
                format!("type fixedValue; value uniform {t0};"),
            ),
            FaceCondition::Pressure { pressure, .. } if pressure > 0.0 => {
                let p0 = num(gas.ambient + pressure);
                (
                    "type pressureInletOutletVelocity; value uniform (0 0 0);".to_string(),
                    format!("type totalPressure; rho none; psi thermo:psi; gamma {g}; p0 uniform {p0}; value uniform {p0};"),
                    format!("type totalTemperature; gamma {g}; T0 uniform {t0}; value uniform {t0};"),
                )
            }
            FaceCondition::Pressure { pressure, .. } => {
                let pa = num(gas.ambient + pressure);
                (
                    "type inletOutlet; inletValue uniform (0 0 0); value uniform (0 0 0);".to_string(),
                    format!(
                        "type waveTransmissive; field p; phi phi; rho rho; psi thermo:psi; gamma {g}; fieldInf {pa}; lInf {}; value uniform {pa};",
                        num(far)
                    ),
                    format!("type inletOutlet; inletValue uniform {t0}; value uniform {t0};"),
                )
            }
            FaceCondition::Slip {} | FaceCondition::Periodic {} => {
                ("type symmetryPlane;".into(), "type symmetryPlane;".into(), "type symmetryPlane;".into())
            }
            FaceCondition::Wall {} => ("type noSlip;".into(), "type zeroGradient;".into(), "type zeroGradient;".into()),
            FaceCondition::MovingWall { velocity } => {
                (format!("type fixedValue; value uniform {};", vector(velocity)), "type zeroGradient;".into(), "type zeroGradient;".into())
            }
        };
        let n = face.name();
        let _ = writeln!(u, "    {n} {{ {fu} }}");
        let _ = writeln!(p, "    {n} {{ {fp} }}");
        let _ = writeln!(t, "    {n} {{ {ft} }}");
    }
    for (_, regions) in surfaces {
        for (_, patch, condition) in regions {
            let fu = match condition {
                SurfaceCondition::Wall {} => "type noSlip;".to_string(),
                SurfaceCondition::MovingWall { velocity } => format!("type fixedValue; value uniform {};", vector(*velocity)),
            };
            let _ = writeln!(u, "    {patch} {{ {fu} }}");
            let _ = writeln!(p, "    {patch} {{ type zeroGradient; }}");
            let _ = writeln!(t, "    {patch} {{ type zeroGradient; }}");
        }
    }
    let file = |class: &str, object: &str, dim: &str, internal: String, b: &str| {
        field_file(class, object, &format!("dimensions {dim};\n\ninternalField uniform {internal};\n\nboundaryField\n{{\n{b}}}\n"))
    };
    [
        ("U", file("volVectorField", "U", "[0 1 -1 0 0 0 0]", vector(d.initial_velocity), &u)),
        ("p", file("volScalarField", "p", "[1 -1 -2 0 0 0 0]", ambient, &p)),
        ("T", file("volScalarField", "T", "[0 0 0 1 0 0 0]", t0, &t)),
    ]
}

fn thermophysical(gas: &Gas, viscosity: f64) -> String {
    dictionary(
        "thermophysicalProperties",
        &format!(
            "thermoType\n{{\n    type hePsiThermo;\n    mixture pureMixture;\n    transport const;\n    thermo hConst;\n    equationOfState perfectGas;\n    specie specie;\n    energy sensibleInternalEnergy;\n}}\n\n\
mixture\n{{\n    specie {{ molWeight {}; }}\n    thermodynamics {{ Cp {}; Hf 0; }}\n    transport {{ mu {}; Pr {}; }}\n}}\n",
            num(gas.molar_mass * 1000.0),
            num(gas.specific_heat),
            num(viscosity),
            num(gas.prandtl)
        ),
    )
}

const GAS_SCHEMES: &str = "fluxScheme Kurganov;\n\
ddtSchemes { default Euler; }\n\
gradSchemes { default Gauss linear; }\n\
divSchemes { default none; div(tauMC) Gauss linear; }\n\
laplacianSchemes { default Gauss linear corrected; }\n\
interpolationSchemes\n{\n    default linear;\n    reconstruct(rho) vanLeer;\n    reconstruct(U) vanLeerV;\n    reconstruct(T) vanLeer;\n}\n\
snGradSchemes { default corrected; }\n";

const GAS_SOLUTION: &str = "solvers\n{\n    \"(rho|rhoU|rhoE)\" { solver diagonal; }\n\
    U { solver smoothSolver; smoother GaussSeidel; nSweeps 2; tolerance 1e-9; relTol 0.01; }\n\
    e { $U; tolerance 1e-10; relTol 0; }\n}\n";

/// OpenFOAM's Smagorinsky model takes `Ck` and `Ce` with `Cs² = Ck √(Ck/Ce)`; the model's
/// constant is `Cs`, so `Ck = (Cs² √Ce)^(2/3)` at OpenFOAM's usual `Ce` 1.048.
const CE: f64 = 1.048;

fn turbulence_properties(t: Option<&Turbulence>, notes: &mut Vec<String>) -> Result<String> {
    let body = match t {
        None => "simulationType laminar;\n".to_string(),
        Some(Turbulence::Smagorinsky { constant, van_driest, wall_model, .. }) => {
            if !(constant.is_finite() && *constant > 0.0) {
                return Err(unsupported("a Smagorinsky constant that is not positive"));
            }
            if wall_model.is_some() {
                notes.push("the wall model is not mapped: walls are resolved no-slip".into());
            }
            let ck = (constant * constant * CE.sqrt()).powf(2.0 / 3.0);
            let delta = match van_driest {
                Some(a) => format!(
                    "delta vanDriest;\n    vanDriestCoeffs\n    {{\n        delta cubeRootVol;\n        cubeRootVolCoeffs {{ deltaCoeff 1; }}\n        Aplus {};\n        Cdelta 0.158;\n    }}",
                    num(*a)
                ),
                None => "delta cubeRootVol;\n    cubeRootVolCoeffs { deltaCoeff 1; }".into(),
            };
            format!(
                "simulationType LES;\n\nLES\n{{\n    LESModel Smagorinsky;\n    turbulence on;\n    printCoeffs on;\n    {delta}\n    SmagorinskyCoeffs {{ Ck {}; Ce {}; }}\n}}\n",
                num(ck),
                num(CE)
            )
        }
        Some(Turbulence::Vreman { .. }) => return Err(unsupported("Vreman's subgrid model (not mapped yet)")),
    };
    Ok(dictionary("turbulenceProperties", &body))
}

fn control(solver: &str, duration: f64, frames: u32, courant: f64) -> String {
    let interval = duration / f64::from(frames);
    let first = (interval / 1000.0).min(1e-6);
    // A running mean of the incompressible LES from a third of the run; a compressible start-up
    // has no such state.
    let average = if solver == "pimpleFoam" {
        format!(
            "    fieldAverage\n    {{\n        type fieldAverage;\n        libs (fieldFunctionObjects);\n        timeStart {};\n        writeControl writeTime;\n\
        fields ( U {{ mean on; prime2Mean on; base time; }} p {{ mean on; prime2Mean off; base time; }} );\n    }}\n",
            num(duration / 3.0)
        )
    } else {
        String::new()
    };
    dictionary(
        "controlDict",
        &format!(
            "application {solver};\nstartFrom latestTime;\nstartTime 0;\nstopAt endTime;\nendTime {};\ndeltaT {};\n\
writeControl adjustable;\nwriteInterval {};\npurgeWrite 0;\nwriteFormat ascii;\nwriteCompression off;\nwritePrecision 8;\ntimeFormat general;\ntimePrecision 8;\n\
runTimeModifiable true;\nadjustTimeStep yes;\nmaxCo {};\n\n\
functions\n{{\n{average}}}\n",
            num(duration),
            num(first),
            num(interval),
            num(courant),
        ),
    )
}

const FV_SCHEMES: &str = "ddtSchemes { default backward; }\n\
gradSchemes { default Gauss linear; grad(U) cellLimited Gauss linear 1; }\n\
divSchemes\n{\n    default none;\n    div(phi,U) Gauss linearUpwind grad(U);\n    div((nuEff*dev2(T(grad(U))))) Gauss linear;\n}\n\
laplacianSchemes { default Gauss linear limited corrected 0.33; }\n\
interpolationSchemes { default linear; }\n\
snGradSchemes { default limited corrected 0.33; }\n\
wallDist { method meshWave; }\n";

const FV_SOLUTION: &str = "solvers\n{\n    p { solver GAMG; smoother DIC; tolerance 1e-6; relTol 0.01; }\n    pFinal { $p; relTol 0; }\n\
    \"(U|UFinal)\" { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-6; relTol 0; }\n}\n\n\
PIMPLE\n{\n    momentumPredictor yes;\n    nOuterCorrectors 1;\n    nCorrectors 2;\n    nNonOrthogonalCorrectors 1;\n}\n";

fn allrun(step: &str, solver: &str, processes: u32) -> String {
    let step: String = step.chars().filter(|c| !c.is_control()).collect();
    let mut s = format!(
        "#!/bin/sh\n# FieldCraft fluid step `{step}` as an OpenFOAM case: mesh, solve, reconstruct.\n\
# Run in this directory with OpenFOAM's environment; NP overrides the ranks ({processes}).\n\
set -e\ncd \"$(dirname \"$0\")\"\nNP=${{NP:-{processes}}}\n\
run() {{ log=\"log.$1\"; echo \"$*  > $log\"; \"$@\" > \"$log\" 2>&1; }}\n\
run blockMesh\nrun surfaceFeatureExtract\nrun snappyHexMesh -overwrite\nrun checkMesh\n\
run postProcess -func writeCellCentres -time 0\n"
    );
    let _ = write!(
        s,
        "if [ \"$NP\" -gt 1 ]; then\n    rm -rf processor*\n    run decomposePar -force\n\
    echo \"mpirun -np $NP {solver} -parallel  > log.{solver}\"\n    mpirun -np \"$NP\" {solver} -parallel > log.{solver} 2>&1\n\
    run reconstructPar\n    rm -rf processor*\nelse\n    run {solver}\nfi\n",
    );
    s
}

/// The most values a field file may hold.
pub const MAX_FIELD_VALUES: usize = 50_000_000;

/// The internal field of an OpenFOAM field file (`internalField`), as written in ASCII.
#[derive(Debug, Clone, PartialEq)]
pub enum FieldValues {
    UniformScalar(f64),
    UniformVector([f64; 3]),
    Scalars(Vec<f64>),
    Vectors(Vec<[f64; 3]>),
}

impl FieldValues {
    /// The scalar of cell `i`.
    pub fn scalar(&self, i: usize) -> Option<f64> {
        match self {
            FieldValues::UniformScalar(v) => Some(*v),
            FieldValues::Scalars(v) => v.get(i).copied(),
            _ => None,
        }
    }

    /// The vector of cell `i`.
    pub fn vector(&self, i: usize) -> Option<[f64; 3]> {
        match self {
            FieldValues::UniformVector(v) => Some(*v),
            FieldValues::Vectors(v) => v.get(i).copied(),
            _ => None,
        }
    }

    /// The number of values (`None` when uniform).
    pub fn count(&self) -> Option<usize> {
        match self {
            FieldValues::Scalars(v) => Some(v.len()),
            FieldValues::Vectors(v) => Some(v.len()),
            _ => None,
        }
    }
}

fn read_error(message: impl Into<String>) -> Error {
    Error::Read { format: FORMAT, message: message.into() }
}

/// The internal field of an ASCII field file: `uniform` with a scalar or a vector, or a
/// `nonuniform List<scalar>` / `List<vector>` (binary files are refused by name).
pub fn read_field(text: &str) -> Result<FieldValues> {
    let start = text.find("internalField").ok_or_else(|| read_error("no internalField"))?;
    let header = text.get(..start).unwrap_or("").replace(';', " ; ");
    let words: Vec<&str> = header.split_whitespace().collect();
    if words.windows(2).any(|w| w == ["format", "binary"]) {
        return Err(read_error("a binary field file (write the case with `writeFormat ascii`)"));
    }
    let rest = text.get(start + "internalField".len()..).unwrap_or("");
    let end = rest.find("boundaryField").unwrap_or(rest.len());
    let spaced = rest.get(..end).unwrap_or("").replace('(', " ( ").replace(')', " ) ").replace(';', " ; ");
    let mut t = spaced.split_whitespace();
    let mut next = |what: &str| t.next().ok_or_else(|| read_error(format!("the internal field ends before {what}")));
    let number = |s: &str| {
        s.parse::<f64>()
            .ok()
            .filter(|v| v.is_finite())
            .ok_or_else(|| read_error(format!("`{}` is not a finite number", s.chars().take(32).collect::<String>())))
    };
    match next("its kind")? {
        "uniform" => match next("its value")? {
            "(" => {
                let v = [number(next("a vector")?)?, number(next("a vector")?)?, number(next("a vector")?)?];
                Ok(FieldValues::UniformVector(v))
            }
            s => Ok(FieldValues::UniformScalar(number(s)?)),
        },
        "nonuniform" => {
            let kind = next("its list type")?;
            let n: usize = next("its count")?.parse().map_err(|_| read_error("a list count that is not a whole number"))?;
            if n > MAX_FIELD_VALUES {
                return Err(Error::Limit(format!("{n} values in a field (at most {MAX_FIELD_VALUES})")));
            }
            if next("its list")? != "(" {
                return Err(read_error("a list that does not start with `(`"));
            }
            let out = match kind {
                "List<scalar>" => {
                    let mut v = Vec::with_capacity(n);
                    for _ in 0..n {
                        v.push(number(next("its values end")?)?);
                    }
                    FieldValues::Scalars(v)
                }
                "List<vector>" => {
                    let mut v = Vec::with_capacity(n);
                    for _ in 0..n {
                        if next("a vector")? != "(" {
                            return Err(read_error("a vector that does not start with `(`"));
                        }
                        let x = [number(next("a vector")?)?, number(next("a vector")?)?, number(next("a vector")?)?];
                        if next("a vector")? != ")" {
                            return Err(read_error("a vector of more than three components"));
                        }
                        v.push(x);
                    }
                    FieldValues::Vectors(v)
                }
                other => return Err(read_error(format!("a field of `{}` (scalars and vectors only)", other.chars().take(32).collect::<String>()))),
            };
            if next("its list's end")? != ")" {
                return Err(read_error(format!("a list longer than its count {n}")));
            }
            Ok(out)
        }
        other => Err(read_error(format!("an internal field that is `{}`, not uniform or nonuniform", other.chars().take(32).collect::<String>()))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A 10 × 4 × 4 mm channel (SI) with a 2 mm cube in it, air, a velocity inlet at x_min and a
    /// pressure outlet at x_max.
    fn channel(extra: serde_json::Value) -> Model {
        let mut domain = serde_json::json!({
            "material": "Air", "min": [0.0, 0.0, 0.0], "max": [0.01, 0.004, 0.004], "spacing": 0.0005,
            "faces": {"x_min": {"type": "velocity", "velocity": [2.0, 0.0, 0.0]},
                      "x_max": {"type": "pressure", "pressure": 0.0},
                      "y_min": {"type": "slip"}},
            "obstacles": ["Cube"]
        });
        if let (Some(d), Some(e)) = (domain.as_object_mut(), extra.as_object()) {
            for (k, v) in e {
                d.insert(k.clone(), v.clone());
            }
        }
        let model = serde_json::json!({
            "materials": {"Air": {"density": 1.2, "fluid": {"type": "newtonian", "viscosity": 1.8e-5}}},
            "fluid_domains": {"Flow": domain},
            "steps": [{"name": "Jet", "procedure": {"type": "fluid", "domain": "Flow", "end_time": 0.01, "frames": 5,
                        "turbulence": {"type": "smagorinsky", "constant": 0.16}}, "output": []}]
        });
        serde_json::from_value(model).unwrap()
    }

    /// A cube from 4 to 6 mm in x and 1 to 3 mm in y and z: 12 triangles, a face per side.
    fn cube() -> (Vec<[f64; 3]>, Vec<[u32; 3]>, Vec<u32>, Vec<String>) {
        let (x0, x1, y0, y1) = (0.004, 0.006, 0.001, 0.003);
        let p = vec![[x0, y0, y0], [x1, y0, y0], [x1, y1, y0], [x0, y1, y0], [x0, y0, y1], [x1, y0, y1], [x1, y1, y1], [x0, y1, y1]];
        let quads = [[0, 3, 2, 1], [4, 5, 6, 7], [0, 1, 5, 4], [3, 7, 6, 2], [0, 4, 7, 3], [1, 2, 6, 5]];
        let mut t = Vec::new();
        let mut f = Vec::new();
        for (i, q) in quads.iter().enumerate() {
            t.push([q[0], q[1], q[2]]);
            t.push([q[0], q[2], q[3]]);
            f.extend([i as u32, i as u32]);
        }
        let names = ["F1:-z", "F1:+z", "F1:-y", "F1:+y", "F1:-x", "F1:+x"].map(String::from).to_vec();
        (p, t, f, names)
    }

    fn write(model: &Model, options: &CaseOptions) -> Result<Case> {
        let (p, t, f, n) = cube();
        let o = ObstacleSurface { body: "Cube", positions: &p, triangles: &t, tri_face: &f, faces: &n };
        write_case(model, None, &[o], options)
    }

    #[test]
    fn writes_a_complete_case_with_a_patch_per_face() {
        let case = write(&channel(serde_json::json!({})), &CaseOptions::default()).unwrap();
        let names: Vec<&str> = case.files.keys().map(String::as_str).collect();
        for f in [
            "0/U",
            "0/p",
            "0/nut",
            "Allrun",
            "case.foam",
            "constant/transportProperties",
            "constant/triSurface/Cube.stl",
            "constant/turbulenceProperties",
            "system/blockMeshDict",
            "system/controlDict",
            "system/decomposeParDict",
            "system/fvSchemes",
            "system/fvSolution",
            "system/snappyHexMeshDict",
            "system/surfaceFeatureExtractDict",
        ] {
            assert!(names.contains(&f), "{f} missing: {names:?}");
        }
        assert_eq!(case.background_cells, [10, 4, 4], "a 1 mm background over 10 × 4 × 4 mm");
        let block = &case.files["system/blockMeshDict"];
        assert!(block.contains("hex (0 1 2 3 4 5 6 7) (10 4 4)"), "{block}");
        assert!(block.contains("y_min\n    {\n        type symmetryPlane;"), "{block}");
        assert!(block.contains("y_max\n    {\n        type wall;"), "unset faces are walls: {block}");
        // One STL solid and one patch per face, with names OpenFOAM accepts.
        let stl = &case.files["constant/triSurface/Cube.stl"];
        assert_eq!(stl.matches("endsolid").count(), 6);
        assert_eq!(stl.matches("endfacet").count(), 12);
        assert!(stl.contains("solid F1_mx\n") && stl.contains("solid F1_px\n"), "{stl}");
        assert_eq!(case.patches.get("Cube_F1_mz").map(String::as_str), Some("Cube/F1:-z"));
        assert_eq!(case.patches.len(), 12);
        // Kinematic viscosity and pressure, the inlet velocity, the obstacle walls.
        assert!(case.files["constant/transportProperties"].contains("nu 1.5e-5;"));
        assert!(case.files["0/U"].contains("x_min { type fixedValue; value uniform (2e0 0 0); }"), "{}", case.files["0/U"]);
        assert!(case.files["0/U"].contains("Cube_F1_mz { type noSlip; }"));
        assert!(case.files["0/p"].contains("x_max { type fixedValue; value uniform 0; }"));
        // Smagorinsky 0.16 → Ck ≈ 0.0882.
        let les = &case.files["constant/turbulenceProperties"];
        assert!(les.contains("LESModel Smagorinsky") && les.contains("Ck 8.82"), "{les}");
        // The step's length and frames.
        let control = &case.files["system/controlDict"];
        assert!(control.contains("endTime 1e-2;") && control.contains("writeInterval 2e-3;"), "{control}");
        // The kept point is in the fluid by the outlet, not in the cube.
        assert!(case.location[0] > 0.006 && case.location[0] < 0.01, "{:?}", case.location);
        assert!(case.files["system/snappyHexMeshDict"].contains("level (1 1)"));
    }

    #[test]
    fn the_same_model_gives_the_same_case() {
        let m = channel(serde_json::json!({}));
        assert_eq!(write(&m, &CaseOptions::default()).unwrap(), write(&m, &CaseOptions::default()).unwrap());
    }

    #[test]
    fn named_moving_faces_and_refinement_boxes() {
        let m = channel(serde_json::json!({"surfaces": [{"face": {"body": "Cube", "face": "F1:+z"},
            "condition": {"type": "moving_wall", "velocity": [0.5, 0.0, 0.0]}}]}));
        let options = CaseOptions {
            cell_factor: 4,
            refine_boxes: vec![RefineBox { min: [0.0, 0.0, 0.0], max: [0.01, 0.002, 0.002], spacing: None }],
            processes: 8,
            ..CaseOptions::default()
        };
        let case = write(&m, &options).unwrap();
        assert!(case.files["0/U"].contains("Cube_F1_pz { type fixedValue; value uniform (5e-1 0 0); }"), "{}", case.files["0/U"]);
        let snappy = &case.files["system/snappyHexMeshDict"];
        assert!(snappy.contains("refine0 { type searchableBox;") && snappy.contains("levels ((1e15 2))"), "{snappy}");
        assert!(case.files["system/decomposeParDict"].contains("numberOfSubdomains 8;"));
        assert!(case.files["Allrun"].contains("NP=${NP:-8}"));
    }

    #[test]
    fn refuses_what_it_cannot_map() {
        let cases = [
            (serde_json::json!({"faces": {"x_min": {"type": "periodic"}, "x_max": {"type": "periodic"}}}), "periodic"),
            (serde_json::json!({"faces": {"x_min": {"type": "velocity", "velocity": [1.0, 0.0, 0.0], "profile": "parabolic"}}}), "parabolic"),
            (serde_json::json!({"acceleration": [0.0, 0.0, -9.81]}), "acceleration"),
        ];
        for (extra, what) in cases {
            let e = write(&channel(extra), &CaseOptions::default()).unwrap_err().to_string();
            assert!(e.contains("not supported by OpenFOAM") && e.contains(what), "{e}");
        }
        let e = write(&channel(serde_json::json!({})), &CaseOptions { cell_factor: 3, ..CaseOptions::default() }).unwrap_err();
        assert!(e.to_string().contains("cell factor"), "{e}");
        // An obstacle with no surface.
        let e = write_case(&channel(serde_json::json!({})), None, &[], &CaseOptions::default()).unwrap_err();
        assert!(e.to_string().contains("obstacle `Cube` has no surface"), "{e}");
    }

    #[test]
    fn hostile_surfaces_never_panic() {
        let m = channel(serde_json::json!({}));
        let (p, t, f, n) = cube();
        let bad_index = [[0u32, 1, 99]];
        let short_faces = [0u32; 3];
        let nan = [[f64::NAN; 3]; 8];
        for (pos, tris, faces) in [(&p[..], &bad_index[..], &f[..1]), (&p[..], &t[..], &short_faces[..]), (&nan[..], &t[..], &f[..])] {
            let o = ObstacleSurface { body: "Cube", positions: pos, triangles: tris, tri_face: faces, faces: &n };
            assert!(write_case(&m, None, &[o], &CaseOptions { location: Some([0.009, 0.002, 0.002]), ..CaseOptions::default() }).is_err());
        }
    }

    #[test]
    fn a_compressible_case_holds_total_pressure_inlets_and_wave_outlets() {
        let m = channel(serde_json::json!({"faces": {"x_min": {"type": "pressure", "pressure": 2.99e7},
            "x_max": {"type": "pressure", "pressure": 0.0}, "y_min": {"type": "slip"}}}));
        let gas = Gas { ambient: 1e5, ..Gas::AIR };
        let options = CaseOptions {
            gas: Some(gas),
            refine_boxes: vec![RefineBox { min: [0.003, 0.0, 0.0], max: [0.01, 0.001, 0.001], spacing: Some(0.000_125) }],
            ..CaseOptions::default()
        };
        let case = write(&m, &options).unwrap();
        assert_eq!(case.solver, "rhoCentralFoam");
        assert!((gas.gamma() - 1.4).abs() < 2e-3, "{}", gas.gamma());
        let p = &case.files["0/p"];
        assert!(p.contains("x_min { type totalPressure;") && p.contains("p0 uniform 3e7;"), "{p}");
        assert!(p.contains("x_max { type waveTransmissive;") && p.contains("fieldInf 1e5;"), "{p}");
        assert!(p.contains("internalField uniform 1e5;") && p.contains("[1 -1 -2 0 0 0 0]"), "{p}");
        let t = &case.files["0/T"];
        assert!(t.contains("x_min { type totalTemperature;") && t.contains("Cube_F1_mz { type zeroGradient; }"), "{t}");
        let thermo = &case.files["constant/thermophysicalProperties"];
        assert!(thermo.contains("molWeight 2.89647e1;") && thermo.contains("mu 1.8e-5;"), "{thermo}");
        assert!(!case.files.contains_key("0/nut") && !case.files.contains_key("constant/transportProperties"));
        assert!(case.files["system/controlDict"].contains("application rhoCentralFoam;"));
        assert!(case.files["Allrun"].contains("mpirun -np \"$NP\" rhoCentralFoam -parallel") && case.files["Allrun"].contains("writeCellCentres"));
        // 1 mm background over 0.125 mm: level 3.
        assert!(case.files["system/snappyHexMeshDict"].contains("levels ((1e15 3))"));
        assert!(case.notes.iter().any(|n| n.contains("laminar")), "{:?}", case.notes);
        let bad = CaseOptions { gas: Some(Gas { molar_mass: 0.0, ..Gas::AIR }), ..CaseOptions::default() };
        assert!(write(&m, &bad).is_err());
    }

    #[test]
    fn reads_ascii_internal_fields() {
        let head = "FoamFile\n{\n    format ascii;\n    class volVectorField;\n}\n// * * //\ndimensions [0 1 -1 0 0 0 0];\n\n";
        let vectors = format!(
            "{head}internalField nonuniform List<vector> \n3\n(\n(1 2 3)\n(4 5 6)\n(-1e-3 0 7.5)\n)\n;\n\nboundaryField\n{{\n    a {{ type fixedValue; value uniform (9 9 9); }}\n}}\n"
        );
        assert_eq!(read_field(&vectors).unwrap(), FieldValues::Vectors(vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [-1e-3, 0.0, 7.5]]));
        let scalars = format!("{head}internalField   nonuniform List<scalar> 2(101325 1.5e5);\nboundaryField {{}}\n");
        assert_eq!(read_field(&scalars).unwrap(), FieldValues::Scalars(vec![101_325.0, 1.5e5]));
        assert_eq!(read_field(&format!("{head}internalField uniform (0 0 1);\n")).unwrap(), FieldValues::UniformVector([0.0, 0.0, 1.0]));
        assert_eq!(read_field(&format!("{head}internalField uniform 293.15;\n")).unwrap().scalar(7), Some(293.15));
        for bad in [
            "internalField nonuniform List<scalar> 3(1 2);",
            "internalField nonuniform List<scalar> 1(1 2);",
            "internalField nonuniform List<vector> 1((1 2));",
            "internalField nonuniform List<tensor> 1((1 2 3 4 5 6 7 8 9));",
            "internalField nonuniform List<scalar> 1(nan);",
            "internalField nonuniform List<scalar> 99999999999(1);",
            "no field here",
        ] {
            assert!(read_field(bad).is_err(), "{bad}");
        }
        assert!(read_field("FoamFile { format binary; }\ninternalField uniform 1;").unwrap_err().to_string().contains("binary"));
    }
}
