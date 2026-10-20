//! From a fluid step of the model to lattice units: the grid, the time step, the relaxation
//! times, and the checks that keep the scheme in its stable and accurate range.
//!
//! # Units (Krüger et al. 2017, ch. 7)
//! Lattice units take the spacing `Δx`, the time step `Δt` and the fluid's density `ρ₀` as 1. A
//! velocity converts by `Δx/Δt`, a kinematic viscosity by `Δx²/Δt`, an acceleration by `Δx/Δt²`,
//! a pressure by `ρ₀ Δx²/Δt²` and a force by `ρ₀ Δx⁴/Δt²`. The lattice viscosity fixes the
//! symmetric relaxation time, `ν = c_s² (τ⁺ − ½)` with `c_s² = 1/3`.
//!
//! # Time step
//! Unless the step gives one, `Δt` is the largest that keeps the characteristic lattice velocity
//! at [`TARGET_LATTICE_VELOCITY`] (Mach 0.087: the weakly compressible error is O(Ma²)) and
//! `τ⁺` at most [`AUTO_TAU_MAX`] (beyond it the scheme only gets less accurate in time). The
//! characteristic velocity is the largest of the inlet and wall velocities, the initial velocity,
//! `√(2 Δp / ρ)` across pressure faces and, for a driving acceleration `a`, the smaller of the
//! viscous plane-channel peak `a H² / (8 ν)` (H the box's smallest extent) and the free-fall speed
//! `√(2 a L)` (L its largest).
//!
//! # Limits
//! A run is refused, naming the numbers and what to change, when the lattice Mach number of the
//! characteristic velocity exceeds [`MAX_MACH`] (the scheme is weakly compressible and loses
//! stability toward Ma ≈ 0.5), when `τ⁺ < `[`MIN_TAU`] (the grid cannot resolve this viscosity at
//! this velocity; linear stability analysis, Lallemand & Luo 2000, shows the stable velocity range
//! shrinks to nothing as τ → ½), or when the lattice acceleration exceeds [`MAX_LATTICE_ACCELERATION`].
//! While running, a lattice speed above [`MAX_LATTICE_SPEED`] or a non-finite or non-positive
//! density stops the run with [`crate::Error::Unstable`].

use fieldcraft_mesh::surface::BoundarySurface;
use fieldcraft_model::fluid_zones::DEFAULT_FAN_RELAXATION;
use fieldcraft_model::units::kelvin_above;
use fieldcraft_model::{
    BoxFace, Collision, Directional, FaceCondition, FluidDomain, Model, Procedure, Resistance, SurfaceCondition, ThermalCondition,
};

use crate::voxel::Grid;
use crate::{Error, Limits, Result};

/// The characteristic lattice velocity the automatic time step aims at.
pub const TARGET_LATTICE_VELOCITY: f64 = 0.05;
/// The largest symmetric relaxation time the automatic time step allows.
pub const AUTO_TAU_MAX: f64 = 1.0;
/// The largest lattice Mach number `u Δt / (Δx c_s)` of the characteristic velocity.
pub const MAX_MACH: f64 = 0.3;
/// The smallest symmetric relaxation time τ⁺.
pub const MIN_TAU: f64 = 0.505;
/// The largest lattice acceleration `a Δt² / Δx`.
pub const MAX_LATTICE_ACCELERATION: f64 = 0.01;
/// A lattice speed above this during a run means it went unstable (Mach 0.87).
pub const MAX_LATTICE_SPEED: f64 = 0.5;
/// Below this τ⁺ the regularized collision is warned about: on the flow-past sample's 2.5 mm grid
/// it ran at τ⁺ 0.508 and went unstable at 0.5075, the reflecting inlet and outlet failing
/// (`docs/lbm-collision.md`).
pub const REGULARIZED_TAU_LIMIT: f64 = 0.508;
/// The default two-relaxation-time parameter Λ: bounce-back walls exactly halfway between cells
/// (Ginzburg, Verhaeghe & d'Humières 2008).
pub const DEFAULT_MAGIC: f64 = 3.0 / 16.0;
/// A lattice Mach number above this (with its round-off) is warned about: the compressibility
/// error, of order Ma², becomes visible.
const MACH_WARNING: f64 = 0.1 + 1e-12; // unit-factors: a Mach number, with its round-off

/// The surface of an obstacle body and the persistent name of each of its faces (by the
/// triangles' face tags); prepared by the caller from the geometry.
#[derive(Debug, Clone)]
pub struct Obstacle {
    pub body: String,
    pub surface: BoundarySurface,
    /// Face names by tag; tags without a name are called `face0`, `face1`, ….
    pub faces: Vec<String>,
}

/// An obstacle of the domain as the run sees it.
#[derive(Debug, Clone, PartialEq)]
pub struct ObstacleFaces {
    /// The obstacle's index in the caller's obstacle list.
    pub index: usize,
    /// Every face tag's name.
    pub names: Vec<String>,
    /// Every face tag's wall velocity in lattice units (zero for walls at rest).
    pub velocities: Vec<[f64; 3]>,
    /// Every face tag's thermal condition (fluxes in lattice units); adiabatic without heat.
    pub thermal: Vec<ThermalCondition>,
}

/// Heat in the fluid, in lattice units: a D3Q7 scalar lattice for the temperature
/// (`flow::ThermalState`).
///
/// **The lattice's temperatures are the model's kelvin less [`ThermalSetup::zero`], the fluid
/// domain's `temperature_zero`**: 0 °C, or for a model from before SI storage the zero of the
/// system it was typed in (0 K for m-kg-s-N-Pa, 0 °F for in-lbf-s-psi). The scheme's walls are
/// not invariant under a shift of the temperature scale (a wall's bounce-back carries `T u`),
/// and before SI storage the lattice ran on the numbers as typed: this keeps every model's
/// results as they were, and how a model is shown never changes them. (A change of the scale's
/// step, °F against K, only scales the lattice's temperatures, which the scheme is invariant
/// under to round-off.) Single precision (the frames, the GPU backend) keeps the differences that
/// drive the heat.
#[derive(Debug, Clone, PartialEq)]
pub struct ThermalSetup {
    /// The scalar lattice's relaxation time: its diffusivity is `(τ_T − ½) / 4`.
    pub tau: f64,
    /// The fluid's thermal diffusivity k / (ρ c), model units.
    pub diffusivity: f64,
    pub initial_temperature: f64,
    /// The temperature, in kelvin, the lattice's 0 stands for (see above).
    pub zero: f64,
    /// `−β g` in lattice units: the buoyant acceleration per temperature above the initial one.
    pub buoyancy: [f64; 3],
    /// The box faces' conditions (fluxes in lattice units: temperature per step through a face).
    pub faces: [ThermalCondition; 6],
    /// Power per lattice heat flow (a temperature per step through one cell face): ρ c Δx³ / Δt.
    pub heat_scale: f64,
}

/// A fluid step in lattice units.
#[derive(Debug, Clone, PartialEq)]
pub struct Setup {
    pub step: String,
    pub domain: FluidDomain,
    pub grid: Grid,
    /// ρ₀ and the kinematic viscosity ν = μ / ρ₀, in model units.
    pub density: f64,
    pub viscosity: f64,
    pub time_step: f64,
    pub steps: u64,
    pub frames: u32,
    /// τ⁺ and τ⁻ (Λ = (τ⁺ − ½)(τ⁻ − ½)).
    pub tau_plus: f64,
    pub tau_minus: f64,
    pub magic: f64,
    /// The collision operator.
    pub collision: Collision,
    /// The characteristic velocity in model and lattice units, and its lattice Mach number.
    pub velocity: f64,
    pub lattice_velocity: f64,
    pub mach: f64,
    /// The body acceleration and initial velocity in lattice units.
    pub acceleration: [f64; 3],
    pub initial_velocity: [f64; 3],
    /// The domain's obstacles, in its order.
    pub obstacles: Vec<ObstacleFaces>,
    /// Heat in the fluid, when the domain has thermal settings.
    pub thermal: Option<ThermalSetup>,
    /// The stop rule, when the step runs to steady state.
    pub steady: Option<crate::steady::SteadySetup>,
    /// The domain's porous zones, then its fans, in lattice units.
    pub zones: Vec<ZoneSetup>,
    pub warnings: Vec<String>,
}

/// What a zone of the domain does, in lattice units.
#[derive(Debug, Clone, PartialEq)]
pub enum ZoneKind {
    /// Porous resistance: the acceleration `−α u − β |u| u` (α = ν/K, β = C₂/2).
    Porous { alpha: f64, beta: f64 },
    /// A fan pushing along its axis with `sign`; its curve's `[Q, Δp]` points in lattice units
    /// (cell volumes per step, `ρ₀ Δx²/Δt²`) and its relaxation per step.
    Fan { sign: f64, curve: Vec<[f64; 2]>, relaxation: f64 },
}

/// A fan or porous zone: the cells whose centres its box contains, `lo..=hi` along each axis.
#[derive(Debug, Clone, PartialEq)]
pub struct ZoneSetup {
    pub name: String,
    pub kind: ZoneKind,
    pub lo: [usize; 3],
    pub hi: [usize; 3],
    /// The axis the flow goes through it along (its thickness and its reported flow).
    pub axis: usize,
}

impl ZoneSetup {
    /// The cell layers along its axis.
    pub fn layers(&self) -> usize {
        let a = self.axis.min(2);
        self.hi[a] + 1 - self.lo[a]
    }

    pub fn is_fan(&self) -> bool {
        matches!(self.kind, ZoneKind::Fan { .. })
    }
}

/// The cells whose centres lie in `[min, max]`, per axis `lo..=hi`, or `None` when there are none.
fn cell_range(grid: &Grid, min: [f64; 3], max: [f64; 3]) -> Option<([usize; 3], [usize; 3])> {
    const TOL: f64 = 1e-9; // unit-factors: a fraction of a cell
    let (mut lo, mut hi) = ([0usize; 3], [0usize; 3]);
    for a in 0..3 {
        let h = grid.spacing;
        let first = ((min[a] - grid.origin[a]) / h - 0.5 - TOL).ceil().max(0.0);
        let last = ((max[a] - grid.origin[a]) / h - 0.5 + TOL).floor().min(grid.cells[a] as f64 - 1.0);
        if !(first <= last) {
            return None;
        }
        (lo[a], hi[a]) = (first as usize, last as usize);
    }
    Some((lo, hi))
}

/// The domain's porous zones and fans in lattice units.
fn zones(domain: &FluidDomain, grid: &Grid, viscosity: f64, dt: f64, density: f64) -> Result<Vec<ZoneSetup>> {
    let h = grid.spacing;
    let pressure_scale = density * h * h / (dt * dt);
    let volume_flow_scale = h * h * h / dt;
    let place = |kind: &str, name: &str, min, max| {
        cell_range(grid, min, max).ok_or_else(|| {
            Error::Invalid(format!("{kind} `{name}` contains no cell centre: make it at least one spacing ({h}) thick and put it inside the box"))
        })
    };
    let mut out = Vec::with_capacity(domain.porous.len() + domain.fans.len());
    for z in &domain.porous {
        let (lo, hi) = place("porous zone", &z.name, z.min, z.max)?;
        let axis = z.axis.index();
        let layers = (hi[axis] + 1 - lo[axis]) as f64;
        let (alpha, beta) = match z.resistance {
            Resistance::DarcyForchheimer { permeability, inertial } => (viscosity * dt / permeability, 0.5 * inertial * h),
            Resistance::LossCoefficient { coefficient } => (0.0, 0.5 * coefficient / layers),
        };
        if !(alpha.is_finite() && beta.is_finite()) {
            return Err(Error::Invalid(format!("porous zone `{}`: its resistance is too large for this lattice", z.name)));
        }
        out.push(ZoneSetup { name: z.name.clone(), kind: ZoneKind::Porous { alpha, beta }, lo, hi, axis });
    }
    for f in &domain.fans {
        let (lo, hi) = place("fan", &f.name, f.min, f.max)?;
        let (axis, sign) = f.direction.axis_sign();
        let curve: Vec<[f64; 2]> = f.curve.iter().map(|p| [p[0] / volume_flow_scale, p[1] / pressure_scale]).collect();
        if curve.iter().flatten().any(|v| !v.is_finite()) {
            return Err(Error::Invalid(format!("fan `{}`: its curve does not fit this lattice", f.name)));
        }
        let relaxation = f.relaxation.unwrap_or(DEFAULT_FAN_RELAXATION);
        out.push(ZoneSetup { name: f.name.clone(), kind: ZoneKind::Fan { sign, curve, relaxation }, lo, hi, axis });
    }
    Ok(out)
}

fn norm(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

/// The fluid step named `step` (or the only fluid step when `step` is `None`) and its domain.
pub fn fluid_step<'m>(model: &'m Model, step: Option<&str>) -> Result<(&'m fieldcraft_model::Step, &'m FluidDomain)> {
    let s = model.select_step(step, |p| p.fluid_domain().is_some(), "fluid")?;
    let name = s.procedure.fluid_domain().unwrap_or_default();
    let domain = model.fluid_domain(name).ok_or_else(|| Error::NotFound(format!("fluid domain `{name}`")))?;
    Ok((s, domain))
}

/// A thermal condition in lattice units: a flux (power per area) becomes the temperature it adds
/// to the cell behind the face per step.
fn lattice_thermal(c: &ThermalCondition, flux_scale: f64, zero: f64) -> Result<ThermalCondition> {
    Ok(match *c {
        ThermalCondition::Temperature { value } => ThermalCondition::Temperature { value: kelvin_above(value, zero) },
        ThermalCondition::Flux { value } => {
            let lattice = value * flux_scale;
            if !lattice.is_finite() {
                return Err(Error::Invalid(format!("a heat flux of {value} is too large for this lattice")));
            }
            ThermalCondition::Flux { value: lattice }
        }
        other => other,
    })
}

impl Setup {
    /// Convert a fluid step, checking it against the stability limits and `limits`.
    pub fn new(model: &Model, step: Option<&str>, obstacles: &[Obstacle], limits: &Limits) -> Result<Setup> {
        // The model checked the domain and its material when they were defined.
        let (s, domain) = fluid_step(model, step)?;
        let domain = domain.clone();
        let Procedure::Fluid { end_time, steps, time_step, frames, magic, collision, steady, .. } = &s.procedure else {
            return Err(Error::Invalid(format!("step `{}` is not a fluid step", s.name)));
        };
        let material = model.material(&domain.material);
        let (Some(density), Some(fluid)) = (material.and_then(|m| m.density), material.and_then(|m| m.fluid)) else {
            return Err(Error::Invalid(format!("material `{}` needs a density and fluid behaviour", domain.material)));
        };
        let viscosity = fluid.viscosity() / density;
        if !(viscosity > 0.0 && viscosity.is_finite()) {
            return Err(Error::Invalid(format!("the kinematic viscosity μ/ρ of material `{}` is {viscosity}", domain.material)));
        }

        // The grid.
        let counts = domain.cells()?;
        let mut cells = [0usize; 3];
        let mut total: u64 = 1;
        for a in 0..3 {
            cells[a] = usize::try_from(counts[a]).map_err(|_| Error::Limit("the grid is too large".into()))?;
            total = total.saturating_mul(counts[a]);
        }
        if total > limits.max_cells {
            return Err(Error::Limit(format!(
                "the domain has {total} cells ({} × {} × {}), above the cap of {}; coarsen the spacing or raise `max_cells`",
                cells[0], cells[1], cells[2], limits.max_cells
            )));
        }
        for face in BoxFace::ALL {
            let open = matches!(domain.faces.get(face), FaceCondition::Velocity { .. } | FaceCondition::Pressure { .. });
            if open && cells[face.axis()] < 2 {
                return Err(Error::Invalid(format!("box face {} is an inlet or outlet, so the box needs at least 2 cells across it", face.name())));
            }
        }
        let grid = Grid { origin: domain.min, spacing: domain.spacing, cells };
        let h = domain.spacing;

        // The characteristic velocity.
        let mut u_char: f64 = norm(domain.initial_velocity);
        let mut pressures = Vec::new();
        for face in BoxFace::ALL {
            match *domain.faces.get(face) {
                FaceCondition::MovingWall { velocity } | FaceCondition::Velocity { velocity, .. } => u_char = u_char.max(norm(velocity)),
                FaceCondition::Pressure { pressure, .. } => pressures.push(pressure),
                FaceCondition::Wall {} | FaceCondition::Slip {} | FaceCondition::Periodic {} => {}
            }
        }
        for sb in &domain.surfaces {
            if let SurfaceCondition::MovingWall { velocity } = sb.condition {
                u_char = u_char.max(norm(velocity));
            }
        }
        for fan in &domain.fans {
            u_char = u_char.max((2.0 * fan.max_pressure_rise() / density).sqrt());
        }
        let lo = pressures.iter().copied().fold(f64::INFINITY, f64::min);
        let hi = pressures.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        if hi > lo {
            u_char = u_char.max((2.0 * (hi - lo) / density).sqrt());
        }
        let a = norm(domain.acceleration);
        if a > 0.0 {
            let extents: Vec<f64> = (0..3).map(|k| domain.max[k] - domain.min[k]).collect();
            let small = extents.iter().copied().fold(f64::INFINITY, f64::min);
            let large = extents.iter().copied().fold(0.0, f64::max);
            u_char = u_char.max((a * small * small / (8.0 * viscosity)).min((2.0 * a * large).sqrt()));
        }
        if !u_char.is_finite() {
            return Err(Error::Invalid("the flow's velocity scale is not a finite number".into()));
        }

        // The time step and the number of steps.
        let diffusive = (AUTO_TAU_MAX - 0.5) * h * h / (3.0 * viscosity);
        let auto = if u_char > 0.0 { (TARGET_LATTICE_VELOCITY * h / u_char).min(diffusive) } else { diffusive };
        let (dt, n) = match (*time_step, *end_time, *steps) {
            (Some(dt), Some(t), _) => {
                let n = (t / dt).round();
                if !(n >= 1.0) {
                    return Err(Error::Invalid(format!("the end time {t} is shorter than the time step {dt}")));
                }
                (dt, n)
            }
            (Some(dt), None, Some(n)) => (dt, n as f64),
            (None, Some(t), _) => {
                let n = (t / auto).ceil().max(1.0);
                (t / n, n)
            }
            (None, None, Some(n)) => (auto, n as f64),
            (_, None, None) => return Err(Error::Invalid(format!("step `{}` needs `end_time` or `steps`", s.name))),
        };
        if !(dt > 0.0 && dt.is_finite()) {
            return Err(Error::Invalid(format!("the time step {dt} is not a positive number")));
        }
        if !(n <= limits.max_steps as f64) {
            return Err(Error::Limit(format!(
                "the step needs {n} time steps of {dt}, above the cap of {}; shorten it or raise `max_steps`",
                limits.max_steps
            )));
        }
        let n = n as u64;
        let updates = (total as f64) * (n as f64);
        if !(updates <= limits.max_cell_updates as f64) {
            return Err(Error::Limit(format!(
                "{total} cells for {n} time steps is {updates:e} cell updates, above the cap of {:e}; coarsen, shorten or raise `max_cell_updates`",
                limits.max_cell_updates as f64
            )));
        }

        // Lattice parameters and the stability limits.
        let nu = viscosity * dt / (h * h);
        let tau_plus = 3.0 * nu + 0.5;
        let lattice_velocity = u_char * dt / h;
        let mach = lattice_velocity / crate::lattice::CS2.sqrt();
        if mach > MAX_MACH {
            return Err(Error::Invalid(format!(
                "the lattice Mach number is {mach:.3} (velocity {u_char} with time step {dt} and spacing {h}), above {MAX_MACH}: lower the time step"
            )));
        }
        if !(tau_plus >= MIN_TAU) {
            return Err(Error::Invalid(format!(
                "the relaxation time τ = {tau_plus:.5} is below {MIN_TAU}: the grid cannot resolve viscosity {viscosity} at this time step (cell Reynolds number {:.1}); refine the spacing or raise the time step",
                if nu > 0.0 { lattice_velocity / nu } else { f64::INFINITY }
            )));
        }
        let magic_given = magic.is_some();
        let magic = magic.unwrap_or(DEFAULT_MAGIC);
        let tau_minus = 0.5 + magic / (tau_plus - 0.5);
        let acceleration = domain.acceleration.map(|c| c * dt * dt / h);
        if norm(acceleration) > MAX_LATTICE_ACCELERATION {
            return Err(Error::Invalid(format!(
                "the lattice acceleration is {:.3e}, above {MAX_LATTICE_ACCELERATION}: lower the time step",
                norm(acceleration)
            )));
        }

        // Heat: the scalar lattice's relaxation time from the diffusivity, the conditions in
        // lattice units (a flux is a temperature per step through a cell face).
        let mut warnings = Vec::new();
        let thermal = match &domain.thermal {
            None => None,
            Some(t) => {
                let props = material.and_then(|m| m.thermal.as_ref());
                let (Some(Directional::Isotropic(k)), Some(c)) = (props.and_then(|p| p.conductivity), props.and_then(|p| p.specific_heat)) else {
                    return Err(Error::Invalid(format!(
                        "material `{}` needs an isotropic conductivity and a specific heat for heat in the fluid",
                        domain.material
                    )));
                };
                let rho_c = density * c;
                let diffusivity = k / rho_c;
                if !(diffusivity > 0.0 && diffusivity.is_finite()) {
                    return Err(Error::Invalid(format!("the thermal diffusivity k / (ρ c) of material `{}` is {diffusivity}", domain.material)));
                }
                let alpha = diffusivity * dt / (h * h);
                let tau = 4.0 * alpha + 0.5;
                if !(tau >= MIN_TAU) {
                    return Err(Error::Invalid(format!(
                        "the thermal relaxation time τ_T = {tau:.5} is below {MIN_TAU}: the grid cannot resolve diffusivity {diffusivity} at this time step (cell Péclet number {:.1}); refine the spacing or raise the time step",
                        if alpha > 0.0 { lattice_velocity / alpha } else { f64::INFINITY }
                    )));
                }
                if tau > 2.0 {
                    warnings.push(format!("the thermal relaxation time τ_T = {tau:.3} is large: heat diffuses fast on this grid and time step; a shorter time step is more accurate"));
                }
                let flux_scale = dt / (rho_c * h);
                // The temperature the lattice counts from (the model's, never how it is shown).
                let zero = t.temperature_zero;
                let mut faces = [ThermalCondition::Adiabatic {}; 6];
                for face in BoxFace::ALL {
                    let c = t.faces.get(face);
                    if matches!(c, ThermalCondition::Flux { .. })
                        && matches!(
                            domain.faces.get(face),
                            FaceCondition::Velocity { .. } | FaceCondition::Pressure { .. } | FaceCondition::Periodic {}
                        )
                    {
                        return Err(Error::Invalid(format!(
                            "box face {}: a heat flux goes on a wall, not on an inlet, outlet or periodic face",
                            face.name()
                        )));
                    }
                    if let Some(slot) = faces.get_mut(face.index()) {
                        *slot = lattice_thermal(c, flux_scale, zero)?;
                    }
                }
                let beta = t.expansion.unwrap_or(0.0);
                Some(ThermalSetup {
                    tau,
                    diffusivity,
                    initial_temperature: kelvin_above(t.initial_temperature, zero),
                    zero,
                    buoyancy: t.gravity.map(|g| -beta * g * dt * dt / h),
                    faces,
                    heat_scale: rho_c * h * h * h / dt,
                })
            }
        };

        // Obstacles, their face names and their faces' wall velocities.
        let mut domain_obstacles = Vec::with_capacity(domain.obstacles.len());
        for body in &domain.obstacles {
            let (index, o) = obstacles
                .iter()
                .enumerate()
                .find(|(_, o)| &o.body == body)
                .ok_or_else(|| Error::NotFound(format!("obstacle body `{body}` has no surface (open the design or import its surface)")))?;
            let tags = o.surface.tri_face().iter().copied().max().map_or(1, |m| m as usize + 1).max(o.faces.len());
            let names: Vec<String> = (0..tags).map(|t| o.faces.get(t).cloned().unwrap_or_else(|| format!("face{t}"))).collect();
            let mut velocities = vec![[0.0; 3]; tags];
            let mut thermal_faces = vec![ThermalCondition::Adiabatic {}; tags];
            if let (Some(t), Some(ts)) = (&domain.thermal, &thermal) {
                // power per area → temperature per step through a cell face, as for the box faces
                let flux_scale = h * h / ts.heat_scale;
                for sb in t.surfaces.iter().filter(|sb| &sb.face.body == body) {
                    let mut hit = false;
                    for (tag, name) in names.iter().enumerate() {
                        if sb.face.selects(&fieldcraft_model::GeometryRef::face(body.as_str(), name.as_str())?) {
                            hit = true;
                            if let Some(slot) = thermal_faces.get_mut(tag) {
                                *slot = lattice_thermal(&sb.condition, flux_scale, ts.zero)?;
                            }
                        }
                    }
                    if !hit {
                        return Err(Error::NotFound(format!("the {} is not a face of the obstacle's surface", sb.face)));
                    }
                }
            }
            for sb in domain.surfaces.iter().filter(|sb| &sb.face.body == body) {
                let mut hit = false;
                for (tag, name) in names.iter().enumerate() {
                    if sb.face.selects(&fieldcraft_model::GeometryRef::face(body.as_str(), name.as_str())?) {
                        hit = true;
                        if let (Some(slot), SurfaceCondition::MovingWall { velocity }) = (velocities.get_mut(tag), sb.condition) {
                            *slot = velocity.map(|c| c * dt / h);
                        }
                    }
                }
                if !hit {
                    return Err(Error::NotFound(format!("the {} is not a face of the obstacle's surface", sb.face)));
                }
            }
            domain_obstacles.push(ObstacleFaces { index, names, velocities, thermal: thermal_faces });
        }
        let zones = zones(&domain, &grid, viscosity, dt, density)?;
        if mach > MACH_WARNING {
            warnings.push(format!(
                "lattice Mach number {mach:.3}: the compressibility error, of order Ma² ≈ {:.1}%, may be visible",
                100.0 * mach * mach
            ));
        }
        let collision = collision.unwrap_or_default();
        if tau_plus < 0.51 && collision == Collision::Trt {
            warnings.push(format!(
                "τ = {tau_plus:.4} is close to ½: under-resolved flows may go unstable; a finer spacing or the regularized collision (`collision: \"regularized\"`) is safer"
            ));
        }
        if tau_plus < REGULARIZED_TAU_LIMIT && collision == Collision::Regularized {
            warnings.push(format!(
                "τ = {tau_plus:.4} is below {REGULARIZED_TAU_LIMIT}, where even the regularized collision went unstable on a coarse bluff-body grid: the inlet and outlet still reflect pressure waves there; a finer spacing is safer"
            ));
        }
        if collision == Collision::Regularized && magic_given {
            warnings.push("the two-relaxation-time parameter Λ (`magic`) has no effect with the regularized collision".into());
        }

        Ok(Setup {
            step: s.name.clone(),
            grid,
            density,
            viscosity,
            time_step: dt,
            steps: n,
            frames: *frames,
            tau_plus,
            tau_minus,
            magic,
            collision,
            velocity: u_char,
            lattice_velocity,
            mach,
            acceleration,
            initial_velocity: domain.initial_velocity.map(|c| c * dt / h),
            obstacles: domain_obstacles,
            thermal,
            steady: steady.as_ref().map(|r| crate::steady::SteadySetup::new(r, n, dt)).transpose()?,
            zones,
            warnings,
            domain,
        })
    }

    /// Model units per lattice unit of volume flow (cell volumes per step).
    pub fn volume_flow_scale(&self) -> f64 {
        self.grid.spacing.powi(3) / self.time_step
    }

    /// Model units per lattice unit of velocity, pressure and force.
    pub fn velocity_scale(&self) -> f64 {
        self.grid.spacing / self.time_step
    }

    pub fn pressure_scale(&self) -> f64 {
        self.density * self.grid.spacing * self.grid.spacing / (self.time_step * self.time_step)
    }

    pub fn force_scale(&self) -> f64 {
        self.pressure_scale() * self.grid.spacing * self.grid.spacing
    }

    /// A velocity in lattice units.
    pub fn to_lattice(&self, v: [f64; 3]) -> [f64; 3] {
        let k = 1.0 / self.velocity_scale();
        v.map(|c| c * k)
    }

    /// The lattice density of a gauge pressure: `1 + p / (c_s² ρ₀ Δx²/Δt²)`.
    pub fn lattice_density(&self, pressure: f64) -> f64 {
        1.0 + pressure / (crate::lattice::CS2 * self.pressure_scale())
    }

    /// The gauge pressure of a lattice density (the inverse of [`Setup::lattice_density`]).
    pub fn gauge_pressure(&self, rho: f64) -> f64 {
        (rho - 1.0) * crate::lattice::CS2 * self.pressure_scale()
    }
}
