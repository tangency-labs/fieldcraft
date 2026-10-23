//! Fluid domains: the box a fluid step solves the flow in.
//!
//! A [`FluidDomain`] is an axis-aligned box filled with a fluid material and divided into cubic
//! cells of one spacing. Each of the six box faces carries a [`FaceCondition`] (a wall, a moving
//! or slip wall, a velocity inlet, a pressure outlet, or a periodic pair). Solid bodies of the
//! geometry placed in the box ([`FluidDomain::obstacles`], by body name) take their cells out of
//! the fluid; their faces are no-slip walls unless a [`SurfaceBoundary`] names one as a moving
//! wall. The forces the fluid exerts on every face are a fluid step's result.
//!
//! Like the rest of the model this is solver-neutral and unit-consistent: lengths, velocities,
//! pressures and accelerations are in the model's unit system, and the solver chooses its own
//! scaled units from them.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

use crate::fluid_zones::{Fan, PorousZone};
use crate::geometry::{Entity, GeometryRef};
use crate::{Error, Material, Model, Result, context, finite, insert_new, not_in_use, positive, remove_existing, replace_existing, require};

/// Most obstacle bodies in one domain.
pub const MAX_OBSTACLES: usize = 1000;

/// Most surface conditions in one domain.
pub const MAX_SURFACE_CONDITIONS: usize = 10_000;

/// Most cells along one side of a domain: the cell count is the product of the three, so this
/// bounds every count derived from a domain before the solver applies its own caps.
pub const MAX_CELLS_PER_SIDE: u64 = 1 << 20;

/// How far the box extent may be from a whole number of spacings, relative to the spacing.
const WHOLE_CELLS_TOLERANCE: f64 = 1e-6; // unit-factors: relative to the spacing

/// A face of the domain box.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum BoxFace {
    XMin,
    XMax,
    YMin,
    YMax,
    ZMin,
    ZMax,
}

impl BoxFace {
    pub const ALL: [BoxFace; 6] = [BoxFace::XMin, BoxFace::XMax, BoxFace::YMin, BoxFace::YMax, BoxFace::ZMin, BoxFace::ZMax];

    /// The axis normal to the face: 0, 1 or 2.
    pub fn axis(self) -> usize {
        match self {
            BoxFace::XMin | BoxFace::XMax => 0,
            BoxFace::YMin | BoxFace::YMax => 1,
            BoxFace::ZMin | BoxFace::ZMax => 2,
        }
    }

    /// Whether the face is at the upper end of its axis.
    pub fn is_max(self) -> bool {
        matches!(self, BoxFace::XMax | BoxFace::YMax | BoxFace::ZMax)
    }

    /// The face across the box.
    pub fn opposite(self) -> BoxFace {
        match self {
            BoxFace::XMin => BoxFace::XMax,
            BoxFace::XMax => BoxFace::XMin,
            BoxFace::YMin => BoxFace::YMax,
            BoxFace::YMax => BoxFace::YMin,
            BoxFace::ZMin => BoxFace::ZMax,
            BoxFace::ZMax => BoxFace::ZMin,
        }
    }

    /// `x_min`, `x_max`, …
    pub fn name(self) -> &'static str {
        match self {
            BoxFace::XMin => "x_min",
            BoxFace::XMax => "x_max",
            BoxFace::YMin => "y_min",
            BoxFace::YMax => "y_max",
            BoxFace::ZMin => "z_min",
            BoxFace::ZMax => "z_max",
        }
    }

    /// The face's index in [`BoxFace::ALL`].
    pub fn index(self) -> usize {
        2 * self.axis() + usize::from(self.is_max())
    }

    /// The face of `axis` (0, 1, 2; larger values mean z) at its upper end when `max`.
    pub fn of(axis: usize, max: bool) -> BoxFace {
        match (axis, max) {
            (0, false) => BoxFace::XMin,
            (0, true) => BoxFace::XMax,
            (1, false) => BoxFace::YMin,
            (1, true) => BoxFace::YMax,
            (_, false) => BoxFace::ZMin,
            (_, true) => BoxFace::ZMax,
        }
    }
}

/// How a velocity inlet's velocity varies over its face.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Profile {
    /// The same velocity everywhere on the face.
    #[default]
    Uniform,
    /// The velocity times `4 s (1 − s)` across each direction of the face that ends at walls
    /// (not periodic), `s` from 0 to 1 across the box: fully developed laminar flow between
    /// plates (one such direction) or in a rectangular duct approximated by the product (two).
    /// The given velocity is the peak.
    Parabolic,
}

/// The condition on a face of the domain box.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case", deny_unknown_fields)]
pub enum FaceCondition {
    /// A no-slip wall at rest.
    Wall {},
    /// A no-slip wall sliding in its own plane (the velocity has no component normal to the
    /// face).
    MovingWall { velocity: [f64; 3] },
    /// A frictionless wall: no flow through it, no shear.
    Slip {},
    /// A velocity inlet (or outlet, if the velocity points out of the box).
    Velocity {
        velocity: [f64; 3],
        #[serde(default, skip_serializing_if = "is_uniform")]
        profile: Profile,
    },
    /// A pressure outlet (or inlet) at this gauge pressure.
    Pressure {
        pressure: f64,
        /// How the flow leaves through the face (default: uniform).
        #[serde(default, skip_serializing_if = "Outflow::is_uniform")]
        outflow: Outflow,
    },
    /// What leaves through this face enters through the opposite one, which must be periodic too.
    Periodic {},
}

/// How the flow leaves through a pressure face.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Outflow {
    /// The face at the given pressure, the flow crossing it normally, its velocity set by that
    /// pressure. Where a wake that has not recovered reaches the face, the speed jumps in the last
    /// cell layer.
    #[default]
    Uniform,
    /// The face at the given pressure, the flow leaving it as it arrives: its velocity, normal and
    /// tangential, is carried out from inside by a convective outflow condition. Where the flow
    /// comes back in, as the uniform outlet. For outlets that a wake or vortices reach.
    Convective,
}

impl Outflow {
    pub fn is_uniform(&self) -> bool {
        *self == Outflow::Uniform
    }
}

fn is_uniform(p: &Profile) -> bool {
    *p == Profile::Uniform
}

impl Default for FaceCondition {
    fn default() -> Self {
        FaceCondition::Wall {}
    }
}

impl FaceCondition {
    /// `wall`, `moving_wall`, `slip`, `velocity`, `pressure` or `periodic`.
    pub fn name(&self) -> &'static str {
        match self {
            FaceCondition::Wall {} => "wall",
            FaceCondition::MovingWall { .. } => "moving_wall",
            FaceCondition::Slip {} => "slip",
            FaceCondition::Velocity { .. } => "velocity",
            FaceCondition::Pressure { .. } => "pressure",
            FaceCondition::Periodic {} => "periodic",
        }
    }

    fn validate(&self, face: BoxFace) -> Result<()> {
        match *self {
            FaceCondition::Wall {} | FaceCondition::Slip {} | FaceCondition::Periodic {} => Ok(()),
            FaceCondition::MovingWall { velocity } => {
                vector("wall velocity", velocity)?;
                if velocity.get(face.axis()).is_some_and(|v| *v != 0.0) {
                    return Err(Error::Invalid(format!(
                        "a moving wall slides in its own plane: the {} component of its velocity must be 0",
                        ["x", "y", "z"].get(face.axis()).copied().unwrap_or("normal")
                    )));
                }
                Ok(())
            }
            FaceCondition::Velocity { velocity, .. } => vector("inlet velocity", velocity).map(|_| ()),
            FaceCondition::Pressure { pressure, .. } => finite("pressure", pressure).map(|_| ()),
        }
    }
}

/// The conditions on the six faces of the box; every face is a wall unless given.
#[derive(Debug, Clone, Copy, Default, PartialEq, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct BoxFaces {
    pub x_min: FaceCondition,
    pub x_max: FaceCondition,
    pub y_min: FaceCondition,
    pub y_max: FaceCondition,
    pub z_min: FaceCondition,
    pub z_max: FaceCondition,
}

impl BoxFaces {
    pub fn get(&self, face: BoxFace) -> &FaceCondition {
        match face {
            BoxFace::XMin => &self.x_min,
            BoxFace::XMax => &self.x_max,
            BoxFace::YMin => &self.y_min,
            BoxFace::YMax => &self.y_max,
            BoxFace::ZMin => &self.z_min,
            BoxFace::ZMax => &self.z_max,
        }
    }

    pub fn get_mut(&mut self, face: BoxFace) -> &mut FaceCondition {
        match face {
            BoxFace::XMin => &mut self.x_min,
            BoxFace::XMax => &mut self.x_max,
            BoxFace::YMin => &mut self.y_min,
            BoxFace::YMax => &mut self.y_max,
            BoxFace::ZMin => &mut self.z_min,
            BoxFace::ZMax => &mut self.z_max,
        }
    }

    /// Whether the faces of `axis` (0, 1, 2) are a periodic pair.
    pub fn periodic(&self, axis: usize) -> bool {
        BoxFace::ALL.iter().any(|f| f.axis() == axis && matches!(self.get(*f), FaceCondition::Periodic {}))
    }
}

/// The condition on a face of an obstacle body.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case", deny_unknown_fields)]
pub enum SurfaceCondition {
    /// A no-slip wall at rest (every obstacle face is one unless told otherwise).
    Wall {},
    /// A no-slip wall moving with this velocity (tangential motion of a rotating or sliding
    /// surface; the body itself does not move).
    MovingWall { velocity: [f64; 3] },
}

/// A condition on one named face of an obstacle body.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SurfaceBoundary {
    /// `{"body": "Body1", "face": "F2:top"}`; the face name also selects its split pieces
    /// (`F2:top#0`, …), as in sets.
    pub face: GeometryRef,
    pub condition: SurfaceCondition,
}

/// What a boundary does to the fluid's temperature.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case", deny_unknown_fields)]
pub enum ThermalCondition {
    /// No heat crosses it (every boundary, unless told otherwise); at an inlet the fluid enters
    /// at the initial temperature, at an outlet it leaves at the temperature it has.
    Adiabatic {},
    /// Held at this temperature (at an inlet, the fluid enters at it).
    Temperature { value: f64 },
    /// A heat flux into the fluid (power per area; negative draws heat out).
    Flux { value: f64 },
}

impl ThermalCondition {
    pub fn name(&self) -> &'static str {
        match self {
            ThermalCondition::Adiabatic {} => "adiabatic",
            ThermalCondition::Temperature { .. } => "temperature",
            ThermalCondition::Flux { .. } => "flux",
        }
    }

    fn validate(&self) -> Result<()> {
        match self {
            ThermalCondition::Adiabatic {} => Ok(()),
            ThermalCondition::Temperature { value } | ThermalCondition::Flux { value } => finite(self.name(), *value).map(|_| ()),
        }
    }
}

impl Default for ThermalCondition {
    fn default() -> Self {
        ThermalCondition::Adiabatic {}
    }
}

/// The thermal condition of each box face.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct BoxThermalFaces {
    pub x_min: ThermalCondition,
    pub x_max: ThermalCondition,
    pub y_min: ThermalCondition,
    pub y_max: ThermalCondition,
    pub z_min: ThermalCondition,
    pub z_max: ThermalCondition,
}

impl BoxThermalFaces {
    pub fn get(&self, face: BoxFace) -> &ThermalCondition {
        match face {
            BoxFace::XMin => &self.x_min,
            BoxFace::XMax => &self.x_max,
            BoxFace::YMin => &self.y_min,
            BoxFace::YMax => &self.y_max,
            BoxFace::ZMin => &self.z_min,
            BoxFace::ZMax => &self.z_max,
        }
    }

    pub fn get_mut(&mut self, face: BoxFace) -> &mut ThermalCondition {
        match face {
            BoxFace::XMin => &mut self.x_min,
            BoxFace::XMax => &mut self.x_max,
            BoxFace::YMin => &mut self.y_min,
            BoxFace::YMax => &mut self.y_max,
            BoxFace::ZMin => &mut self.z_min,
            BoxFace::ZMax => &mut self.z_max,
        }
    }
}

/// A thermal condition on one named face of an obstacle body.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SurfaceThermalBoundary {
    pub face: GeometryRef,
    pub condition: ThermalCondition,
}

/// An obstacle that conducts heat (conjugate heat transfer): its cells carry a temperature that
/// diffuses with its material's conductivity and heat capacity, coupled to the fluid's
/// temperature across every face it shares with the fluid (the temperature and the heat flux are
/// continuous there), and it may generate heat (a component's power).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ConductingSolid {
    /// An obstacle body of the domain.
    pub body: String,
    /// A material with a density, a conductivity (isotropic, or orthotropic along the domain's
    /// axes: a board's in-plane and through-thickness values) and a specific heat.
    pub material: String,
    /// The heat it generates (power), spread evenly over its volume; negative draws heat out.
    #[serde(default, skip_serializing_if = "crate::contact::is_zero")]
    pub power: f64,
    /// The highest temperature the component may reach (kelvin; a junction or case limit): the
    /// solve reports its peak against it. None for no limit.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub limit: Option<f64>,
}

/// A thin layer of thermal interface material between two touching conducting solids (a paste,
/// a pad, a gap filler): not resolved in the lattice, it adds its resistance to every face
/// between cells of the two bodies. Given either as a resistance per area or as the layer's
/// thickness and conductivity (`R = t / k`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ThermalContact {
    /// Two conducting solids of the domain.
    pub bodies: [String; 2],
    /// The resistance per area (K·m²/W).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub resistance: Option<f64>,
    /// The layer's thickness, with its `conductivity`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub thickness: Option<f64>,
    /// The layer's conductivity, with its `thickness`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub conductivity: Option<f64>,
}

impl ThermalContact {
    /// The resistance per area, from whichever form the contact is given in (zero when it is
    /// given in neither: [`ThermalContact::validate`] refuses that).
    pub fn resistance(&self) -> f64 {
        match (self.resistance, self.thickness, self.conductivity) {
            (Some(r), ..) => r,
            (None, Some(t), Some(k)) => t / k,
            _ => 0.0,
        }
    }

    fn validate(&self) -> Result<()> {
        match (self.resistance, self.thickness, self.conductivity) {
            (Some(r), None, None) => {
                positive("resistance", r)?;
            }
            (None, Some(t), Some(k)) => {
                positive("thickness", t)?;
                positive("conductivity", k)?;
                positive("resistance (thickness over conductivity)", t / k)?;
            }
            _ => {
                return Err(Error::Invalid(
                    "give either a `resistance` per area or a layer's `thickness` and `conductivity`, not both or neither".into(),
                ));
            }
        }
        if self.bodies[0] == self.bodies[1] {
            return Err(Error::Invalid(format!("a contact joins two bodies; `{}` is named twice", self.bodies[0])));
        }
        Ok(())
    }

    /// Whether it joins `a` and `b`, in either order.
    pub fn joins(&self, a: &str, b: &str) -> bool {
        let [x, y] = &self.bodies;
        (x == a && y == b) || (x == b && y == a)
    }

    fn name(&self) -> String {
        format!("{}–{}", self.bodies[0], self.bodies[1])
    }
}

/// Heat in the fluid: a temperature field advected with the flow and diffused with the
/// material's diffusivity k / (ρ c), its conditions on the boundaries, and buoyancy.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct FluidThermal {
    /// The temperature everywhere at the start (kelvin; 0 °C when left out), and the reference
    /// of the buoyancy.
    pub initial_temperature: f64,
    /// Gravity for the buoyancy (Boussinesq: a force per mass of −β (T − T₀) g on the fluid);
    /// zero for none. Separate from the domain's `acceleration`, which is a plain body force.
    pub gravity: [f64; 3],
    /// The fluid's volumetric thermal expansion coefficient β (per temperature); needed with
    /// a gravity.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub expansion: Option<f64>,
    pub faces: BoxThermalFaces,
    /// Conditions on named faces of the obstacles (adiabatic unless listed).
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub surfaces: Vec<SurfaceThermalBoundary>,
    /// Obstacles that conduct heat; the others are walls with the conditions in `surfaces`.
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub solids: Vec<ConductingSolid>,
    /// Thermal interface layers between touching conducting solids; touching solids not listed
    /// are in perfect contact.
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub contacts: Vec<ThermalContact>,
    /// The temperature (kelvin) the lattice-Boltzmann scheme's temperature lattice counts from:
    /// its walls are not invariant under a shift of the temperature scale, so this is part of
    /// the model, never of how it is shown. 0 °C, the scale the scheme was verified on; a model
    /// from before SI storage keeps the zero of the system it was typed in (0 K for
    /// m-kg-s-N-Pa, 0 °F for in-lbf-s-psi), so it solves as it did.
    #[serde(skip_serializing_if = "is_celsius_zero")]
    pub temperature_zero: f64,
}

fn is_celsius_zero(t: &f64) -> bool {
    *t == crate::units::CELSIUS_ZERO
}

impl Default for FluidThermal {
    fn default() -> Self {
        FluidThermal {
            initial_temperature: crate::units::CELSIUS_ZERO,
            gravity: [0.0; 3],
            expansion: None,
            faces: BoxThermalFaces::default(),
            surfaces: Vec::new(),
            solids: Vec::new(),
            contacts: Vec::new(),
            temperature_zero: crate::units::CELSIUS_ZERO,
        }
    }
}

impl FluidThermal {
    fn validate(&self, obstacles: &[String]) -> Result<()> {
        finite("initial temperature", self.initial_temperature)?;
        finite("temperature zero", self.temperature_zero)?;
        vector("gravity", self.gravity)?;
        if let Some(b) = self.expansion {
            finite("expansion coefficient", b)?;
            if b < 0.0 {
                return Err(Error::Invalid(format!("the expansion coefficient {b} is negative")));
            }
        }
        if self.gravity != [0.0; 3] && self.expansion.is_none() {
            return Err(Error::Invalid("buoyancy needs the fluid's expansion coefficient (`expansion`) with its gravity".into()));
        }
        for face in BoxFace::ALL {
            self.faces.get(face).validate().map_err(|e| context("box face", face.name(), e))?;
        }
        if self.surfaces.len() > MAX_SURFACE_CONDITIONS {
            return Err(Error::Limit(format!("a fluid domain holds at most {MAX_SURFACE_CONDITIONS} thermal surface conditions")));
        }
        for (i, s) in self.surfaces.iter().enumerate() {
            s.condition.validate()?;
            if !matches!(s.face.entity, Entity::Face(_)) {
                return Err(Error::Invalid(format!("thermal conditions go on faces, not on the {}", s.face)));
            }
            if !obstacles.contains(&s.face.body) {
                return Err(Error::Invalid(format!(
                    "the {} is not on an obstacle of the domain (list body `{}` in `obstacles`)",
                    s.face, s.face.body
                )));
            }
            if self.surfaces.get(..i).is_some_and(|before| before.iter().any(|b| b.face == s.face)) {
                return Err(Error::Invalid(format!("the {} has two thermal conditions", s.face)));
            }
            if self.solids.iter().any(|c| c.body == s.face.body) {
                return Err(Error::Invalid(format!(
                    "the {} is on a conducting solid: its faces exchange heat with the fluid, so they take no thermal condition",
                    s.face
                )));
            }
        }
        for (i, c) in self.solids.iter().enumerate() {
            finite("power", c.power).map_err(|e| context("conducting solid", &c.body, e))?;
            if let Some(limit) = c.limit {
                finite("temperature limit", limit).map_err(|e| context("conducting solid", &c.body, e))?;
                if limit <= 0.0 {
                    return Err(context(
                        "conducting solid",
                        &c.body,
                        Error::Invalid(format!("the temperature limit {limit} K is not above absolute zero")),
                    ));
                }
            }
            if !obstacles.contains(&c.body) {
                return Err(Error::Invalid(format!("conducting solid `{}` is not an obstacle of the domain (list it in `obstacles`)", c.body)));
            }
            if self.solids.get(..i).is_some_and(|before| before.iter().any(|b| b.body == c.body)) {
                return Err(Error::Invalid(format!("conducting solid `{}` is listed twice", c.body)));
            }
        }
        if self.contacts.len() > MAX_SURFACE_CONDITIONS {
            return Err(Error::Limit(format!("a fluid domain holds at most {MAX_SURFACE_CONDITIONS} thermal contacts")));
        }
        for (i, c) in self.contacts.iter().enumerate() {
            c.validate().map_err(|e| context("thermal contact", &c.name(), e))?;
            for body in &c.bodies {
                if !self.solids.iter().any(|s| s.body == *body) {
                    return Err(context(
                        "thermal contact",
                        &c.name(),
                        Error::Invalid(format!("`{body}` is not a conducting solid of the domain (list it in `solids`)")),
                    ));
                }
            }
            if self.contacts.get(..i).is_some_and(|before| before.iter().any(|b| b.joins(&c.bodies[0], &c.bodies[1]))) {
                return Err(context("thermal contact", &c.name(), Error::Invalid("the two bodies have a contact already".into())));
            }
        }
        Ok(())
    }
}

/// What a conducting solid needs of its material: a density, a conductivity and a specific heat
/// (the material keeps them positive).
fn material_fits_solid(m: &Material, name: &str) -> Result<()> {
    let thermal = m.thermal.as_ref();
    if !(thermal.is_some_and(|t| t.conductivity.is_some() && t.specific_heat.is_some()) && m.density.is_some()) {
        return Err(Error::Invalid(format!("material `{name}` needs a density, a positive conductivity and a specific heat to conduct heat")));
    }
    Ok(())
}

/// What a fluid domain needs of its material: a density and fluid behaviour, and with heat a
/// positive isotropic conductivity and specific heat (the diffusivity k / (ρ c)).
fn material_fits_fluid_domain(m: &Material, name: &str, d: &FluidDomain) -> Result<()> {
    if m.fluid.is_none() {
        return Err(Error::Invalid(format!("material `{name}` has no fluid behaviour (a viscosity)")));
    }
    if m.density.is_none() {
        return Err(Error::Invalid(format!("material `{name}` has no density")));
    }
    if d.thermal.is_some() {
        let thermal = m.thermal.as_ref();
        match thermal.and_then(|t| t.conductivity) {
            Some(crate::Directional::Isotropic(k)) if k > 0.0 => {}
            Some(crate::Directional::Isotropic(_)) | None => {
                return Err(Error::Invalid(format!("material `{name}` needs a positive isotropic conductivity for heat in the fluid")));
            }
            Some(crate::Directional::Orthotropic(_)) => {
                return Err(Error::Invalid(format!("material `{name}`: a fluid's conductivity is isotropic")));
            }
        }
        if !thermal.and_then(|t| t.specific_heat).is_some_and(|c| c > 0.0) {
            return Err(Error::Invalid(format!("material `{name}` needs a positive specific heat for heat in the fluid")));
        }
    }
    Ok(())
}

/// Where a boundary condition goes: a face of the box, or a named face of an obstacle.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(untagged)]
pub enum FluidBoundary {
    Box(BoxFace),
    Surface(GeometryRef),
}

/// A fluid domain.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FluidDomain {
    /// A material with a density and fluid behaviour ([`crate::Fluid`]).
    pub material: String,
    /// The box's lower and upper corners.
    pub min: [f64; 3],
    pub max: [f64; 3],
    /// The side of the cubic cells; each extent of the box must be a whole number of spacings.
    pub spacing: f64,
    #[serde(default)]
    pub faces: BoxFaces,
    /// Solid bodies of the geometry in the box (by name): their cells are not fluid.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub obstacles: Vec<String>,
    /// Conditions on named faces of the obstacles (each face once).
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub surfaces: Vec<SurfaceBoundary>,
    /// The fluid's velocity at the start.
    #[serde(default, skip_serializing_if = "is_zero")]
    pub initial_velocity: [f64; 3],
    /// A body acceleration on the fluid (gravity, or a uniform driving force per mass).
    #[serde(default, skip_serializing_if = "is_zero")]
    pub acceleration: [f64; 3],
    /// The chordal tolerance obstacles are triangulated with (default: a tenth of the spacing).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub tolerance: Option<f64>,
    /// Heat in the fluid (a temperature field, buoyancy); `None` for an isothermal flow.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub thermal: Option<FluidThermal>,
    /// Vents, grilles, filters and perforated plates as porous resistance zones.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub porous: Vec<PorousZone>,
    /// Fans: zones that raise the pressure by their pressure–flow curve.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub fans: Vec<Fan>,
}

fn is_zero(v: &[f64; 3]) -> bool {
    v.iter().all(|c| *c == 0.0)
}

fn vector(what: &str, v: [f64; 3]) -> Result<[f64; 3]> {
    for c in v {
        finite(what, c)?;
    }
    Ok(v)
}

impl FluidDomain {
    /// The number of cells along each axis, checked to be whole and within
    /// [`MAX_CELLS_PER_SIDE`].
    pub fn cells(&self) -> Result<[u64; 3]> {
        let h = positive("spacing", self.spacing)?;
        let mut n = [0u64; 3];
        for (axis, slot) in n.iter_mut().enumerate() {
            let (lo, hi) = (self.min.get(axis).copied().unwrap_or(0.0), self.max.get(axis).copied().unwrap_or(0.0));
            let name = ["x", "y", "z"].get(axis).copied().unwrap_or("?");
            if !(lo.is_finite() && hi.is_finite() && hi > lo) {
                return Err(Error::Invalid(format!("the box's {name} range {lo}..{hi} is empty or not finite")));
            }
            let cells = (hi - lo) / h;
            let whole = cells.round();
            if !(whole >= 1.0 && whole <= MAX_CELLS_PER_SIDE as f64) {
                return Err(Error::Invalid(format!("the box spans {cells} spacings along {name}; it must span 1 to {MAX_CELLS_PER_SIDE} cells")));
            }
            if (cells - whole).abs() > WHOLE_CELLS_TOLERANCE * whole.max(1.0) {
                return Err(Error::Invalid(format!("the box's {name} extent {} is not a whole number of spacings {h} ({cells} cells)", hi - lo)));
            }
            *slot = whole as u64;
        }
        Ok(n)
    }

    /// The obstacles that conduct heat (none without thermal settings).
    pub fn solids(&self) -> impl Iterator<Item = &ConductingSolid> {
        self.thermal.iter().flat_map(|t| t.solids.iter())
    }

    /// Check everything that does not need the rest of the model.
    pub fn validate(&self) -> Result<()> {
        for (what, v) in [("min", self.min), ("max", self.max), ("initial velocity", self.initial_velocity), ("acceleration", self.acceleration)] {
            vector(what, v)?;
        }
        self.cells()?;
        for face in BoxFace::ALL {
            let c = self.faces.get(face);
            c.validate(face).map_err(|e| context("box face", face.name(), e))?;
            let other = self.faces.get(face.opposite());
            // Both faces of a pair are visited, so checking one direction finds every broken pair.
            if matches!(c, FaceCondition::Periodic {}) && !matches!(other, FaceCondition::Periodic {}) {
                return Err(Error::Invalid(format!(
                    "box face {} is periodic but {} is {}: periodic faces come in opposite pairs",
                    face.name(),
                    face.opposite().name(),
                    other.name()
                )));
            }
        }
        if self.obstacles.len() > MAX_OBSTACLES {
            return Err(Error::Limit(format!("a fluid domain holds at most {MAX_OBSTACLES} obstacles")));
        }
        for (i, b) in self.obstacles.iter().enumerate() {
            GeometryRef::body(b.as_str())?;
            if self.obstacles.get(..i).is_some_and(|before| before.contains(b)) {
                return Err(Error::Invalid(format!("obstacle `{b}` is listed twice")));
            }
        }
        if let Some(t) = &self.thermal {
            t.validate(&self.obstacles).map_err(|e| context("thermal", "settings", e))?;
        }
        if self.surfaces.len() > MAX_SURFACE_CONDITIONS {
            return Err(Error::Limit(format!("a fluid domain holds at most {MAX_SURFACE_CONDITIONS} surface conditions")));
        }
        for (i, s) in self.surfaces.iter().enumerate() {
            if !matches!(s.face.entity, Entity::Face(_)) {
                return Err(Error::Invalid(format!("surface conditions go on faces, not on the {}", s.face)));
            }
            if !self.obstacles.contains(&s.face.body) {
                return Err(Error::Invalid(format!(
                    "the {} is not on an obstacle of the domain (list body `{}` in `obstacles`)",
                    s.face, s.face.body
                )));
            }
            if let SurfaceCondition::MovingWall { velocity } = s.condition {
                vector("wall velocity", velocity)?;
            }
            if self.surfaces.get(..i).is_some_and(|before| before.iter().any(|o| o.face == s.face)) {
                return Err(Error::Invalid(format!("the {} has two conditions", s.face)));
            }
        }
        if let Some(t) = self.tolerance {
            positive("tolerance", t)?;
        }
        crate::fluid_zones::validate_zones(&self.porous, &self.fans)?;
        Ok(())
    }

    /// Set the condition on one boundary: a box face, or a named face of an obstacle, which
    /// Set one boundary's thermal condition (an adiabatic obstacle face removes its entry,
    /// adiabatic being the default). Needs the domain's thermal settings; the caller validates
    /// the result ([`Model::set_fluid_thermal_boundary`] does).
    pub fn set_thermal_boundary(&mut self, at: FluidBoundary, condition: ThermalCondition) -> Result<()> {
        let Some(t) = &mut self.thermal else {
            return Err(Error::Invalid("the domain has no thermal settings (`thermal`), so no boundary has a thermal condition".into()));
        };
        match at {
            FluidBoundary::Box(face) => *t.faces.get_mut(face) = condition,
            FluidBoundary::Surface(face) => {
                t.surfaces.retain(|s| s.face != face);
                if condition != (ThermalCondition::Adiabatic {}) {
                    t.surfaces.push(SurfaceThermalBoundary { face, condition });
                }
            }
        }
        Ok(())
    }

    /// takes a wall or a moving wall (a wall there removes the entry, walls being the default).
    /// The caller validates the result ([`Model::set_fluid_boundary`] does).
    pub fn set_boundary(&mut self, at: FluidBoundary, condition: FaceCondition) -> Result<()> {
        match at {
            FluidBoundary::Box(face) => *self.faces.get_mut(face) = condition,
            FluidBoundary::Surface(face) => {
                let c = match condition {
                    FaceCondition::Wall {} => SurfaceCondition::Wall {},
                    FaceCondition::MovingWall { velocity } => SurfaceCondition::MovingWall { velocity },
                    other => return Err(Error::Invalid(format!("the {face} takes a wall or a moving wall, not {}", other.name()))),
                };
                self.surfaces.retain(|s| s.face != face);
                if c != (SurfaceCondition::Wall {}) {
                    self.surfaces.push(SurfaceBoundary { face, condition: c });
                }
            }
        }
        Ok(())
    }
}

impl Model {
    pub fn fluid_domains(&self) -> &BTreeMap<String, FluidDomain> {
        &self.fluid_domains
    }

    pub fn fluid_domain(&self, name: &str) -> Option<&FluidDomain> {
        self.fluid_domains.get(name)
    }

    fn check_fluid_domain(&self, d: &FluidDomain) -> Result<()> {
        d.validate()?;
        require(&self.materials, "material", &d.material)?;
        if let Some(m) = self.materials.get(&d.material) {
            material_fits_fluid_domain(m, &d.material, d)?;
        }
        for c in d.solids() {
            require(&self.materials, "material", &c.material)?;
            if let Some(m) = self.materials.get(&c.material) {
                material_fits_solid(m, &c.material).map_err(|e| context("conducting solid", &c.body, e))?;
            }
        }
        Ok(())
    }

    /// The fluid domains `material` must keep fitting when it is replaced.
    pub(crate) fn check_material_fits_fluid_domains(&self, name: &str, material: &Material) -> Result<()> {
        for (domain, d) in &self.fluid_domains {
            if d.material == name {
                material_fits_fluid_domain(material, name, d).map_err(|e| context("fluid domain", domain, e))?;
            }
            for c in d.solids().filter(|c| c.material == name) {
                material_fits_solid(material, name).map_err(|e| context("fluid domain", domain, context("conducting solid", &c.body, e)))?;
            }
        }
        Ok(())
    }

    /// Set the thermal condition on one boundary of a domain.
    pub fn set_fluid_thermal_boundary(&mut self, name: &str, at: FluidBoundary, condition: ThermalCondition) -> Result<()> {
        let mut d = self.fluid_domains.get(name).cloned().ok_or_else(|| Error::NotFound(format!("fluid domain `{name}`")))?;
        d.set_thermal_boundary(at, condition).map_err(|e| context("fluid domain", name, e))?;
        self.replace_fluid_domain(name, d)
    }

    /// Set (or, with `None`, clear) the temperature limit of one of a domain's conducting solids.
    pub fn set_component_limit(&mut self, name: &str, body: &str, limit: Option<f64>) -> Result<()> {
        let mut d = self.fluid_domains.get(name).cloned().ok_or_else(|| Error::NotFound(format!("fluid domain `{name}`")))?;
        let solid = d.thermal.as_mut().and_then(|t| t.solids.iter_mut().find(|c| c.body == body)).ok_or_else(|| {
            context(
                "fluid domain",
                name,
                Error::Invalid(format!(
                    "`{body}` is not a conducting solid of the domain, so it has no temperature to limit (list it in the thermal settings' `solids`)"
                )),
            )
        })?;
        solid.limit = limit;
        self.replace_fluid_domain(name, d)
    }

    pub fn add_fluid_domain(&mut self, name: &str, domain: FluidDomain) -> Result<()> {
        self.check_fluid_domain(&domain).map_err(|e| context("fluid domain", name, e))?;
        insert_new(&mut self.fluid_domains, "fluid domain", name, domain)
    }

    pub fn replace_fluid_domain(&mut self, name: &str, domain: FluidDomain) -> Result<()> {
        self.check_fluid_domain(&domain).map_err(|e| context("fluid domain", name, e))?;
        replace_existing(&mut self.fluid_domains, "fluid domain", name, domain)
    }

    /// The steps that solve a fluid domain, as `step \`name\`` strings.
    pub fn fluid_domain_users(&self, name: &str) -> Vec<String> {
        self.steps.iter().filter(|s| s.procedure.fluid_domain() == Some(name)).map(|s| format!("step `{}`", s.name)).collect()
    }

    pub fn remove_fluid_domain(&mut self, name: &str) -> Result<()> {
        not_in_use("fluid domain", name, self.fluid_domain_users(name))?;
        remove_existing(&mut self.fluid_domains, "fluid domain", name).map(|_| ())
    }

    /// Set the condition on one boundary of a domain.
    pub fn set_fluid_boundary(&mut self, name: &str, at: FluidBoundary, condition: FaceCondition) -> Result<()> {
        let mut d = self.fluid_domains.get(name).cloned().ok_or_else(|| Error::NotFound(format!("fluid domain `{name}`")))?;
        d.set_boundary(at, condition).map_err(|e| context("fluid domain", name, e))?;
        self.replace_fluid_domain(name, d)
    }

    /// The fluid domains that use a material.
    pub(crate) fn fluid_domains_of_material(&self, material: &str) -> impl Iterator<Item = &str> {
        self.fluid_domains.iter().filter(move |(_, d)| d.material == material || d.solids().any(|c| c.material == material)).map(|(n, _)| n.as_str())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Directional, Fluid, Material, Procedure, Step};

    fn water() -> Material {
        Material { density: Some(1.0e-9), fluid: Some(Fluid::Newtonian { viscosity: 1.0e-9 }), ..Material::default() }
    }

    pub(crate) fn channel() -> FluidDomain {
        FluidDomain {
            material: "Water".into(),
            min: [0.0, 0.0, 0.0],
            max: [10.0, 2.0, 1.0],
            spacing: 0.5,
            faces: BoxFaces {
                x_min: FaceCondition::Velocity { velocity: [1.0, 0.0, 0.0], profile: Profile::Parabolic },
                x_max: FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Uniform },
                z_min: FaceCondition::Periodic {},
                z_max: FaceCondition::Periodic {},
                ..BoxFaces::default()
            },
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

    #[test]
    fn thermal_settings_validate_and_need_the_materials_heat_properties() {
        let mut m = Model::default();
        m.add_material("Water", water()).unwrap();
        let thermal = |f: &dyn Fn(&mut FluidThermal)| {
            let mut d = channel();
            let mut t = FluidThermal { initial_temperature: 20.0, ..FluidThermal::default() };
            f(&mut t);
            d.thermal = Some(t);
            d
        };
        // heat in the fluid needs a conductivity and a specific heat on the material
        assert!(m.add_fluid_domain("Flow", thermal(&|_| {})).is_err());
        let mut hot = water();
        hot.thermal = Some(crate::Thermal { conductivity: Some(Directional::Isotropic(0.6)), specific_heat: Some(4180.0), ..Default::default() });
        m.replace_material("Water", hot.clone()).unwrap();
        m.add_fluid_domain("Flow", thermal(&|_| {})).unwrap();
        m.add_fluid_domain(
            "Hot",
            thermal(&|t| {
                t.faces.y_min = ThermalCondition::Temperature { value: 50.0 };
                t.gravity = [0.0, -9810.0, 0.0];
                t.expansion = Some(2e-4);
            }),
        )
        .unwrap();
        for d in [
            thermal(&|t| t.initial_temperature = f64::NAN),
            thermal(&|t| t.gravity = [0.0, -1.0, 0.0]),
            thermal(&|t| t.expansion = Some(-1.0)),
            thermal(&|t| t.faces.x_min = ThermalCondition::Flux { value: f64::INFINITY }),
            thermal(&|t| {
                t.surfaces = vec![SurfaceThermalBoundary {
                    face: GeometryRef::face("Nope", "F1").unwrap(),
                    condition: ThermalCondition::Temperature { value: 1.0 },
                }]
            }),
        ] {
            assert!(m.add_fluid_domain("Bad", d).is_err());
        }
        let mut aniso = hot;
        aniso.thermal =
            Some(crate::Thermal { conductivity: Some(Directional::Orthotropic([1.0; 3])), specific_heat: Some(1.0), ..Default::default() });
        m.add_material("Oil", aniso.clone()).unwrap();
        let mut d = thermal(&|_| {});
        d.material = "Oil".into();
        assert!(m.add_fluid_domain("Aniso", d).is_err(), "a fluid's conductivity is isotropic");
        assert!(m.replace_material("Water", water()).is_err(), "the domains with heat keep needing the heat properties");
        assert!(m.replace_material("Water", aniso).is_err());
        // setting conditions
        m.set_fluid_thermal_boundary("Flow", FluidBoundary::Box(BoxFace::XMin), ThermalCondition::Temperature { value: 80.0 }).unwrap();
        assert_eq!(m.fluid_domain("Flow").unwrap().thermal.as_ref().unwrap().faces.x_min, ThermalCondition::Temperature { value: 80.0 });
        assert!(m.set_fluid_thermal_boundary("Flow", FluidBoundary::Box(BoxFace::XMin), ThermalCondition::Flux { value: f64::NAN }).is_err());
        let mut plain = channel();
        plain.thermal = None;
        m.add_fluid_domain("Plain", plain).unwrap();
        assert!(
            m.set_fluid_thermal_boundary("Plain", FluidBoundary::Box(BoxFace::XMin), ThermalCondition::Adiabatic {}).is_err(),
            "no thermal settings"
        );
    }

    #[test]
    fn conducting_solids_need_an_obstacle_and_a_material_that_conducts() {
        let mut m = Model::default();
        let mut hot = water();
        hot.thermal = Some(crate::Thermal { conductivity: Some(Directional::Isotropic(0.6)), specific_heat: Some(4180.0), ..Default::default() });
        m.add_material("Water", hot).unwrap();
        let alu = Material {
            density: Some(2700.0),
            thermal: Some(crate::Thermal { conductivity: Some(Directional::Isotropic(237.0)), specific_heat: Some(900.0), ..Default::default() }),
            ..Material::default()
        };
        m.add_material("Alu", alu.clone()).unwrap();
        // a board: in-plane and through-thickness conductivities
        let mut fr4 = alu.clone();
        fr4.thermal = Some(crate::Thermal {
            conductivity: Some(Directional::Orthotropic([20.0, 20.0, 0.3])),
            specific_heat: Some(1100.0),
            ..Default::default()
        });
        m.add_material("Board", fr4).unwrap();
        m.add_material("Steel", Material { density: Some(7800.0), ..Material::default() }).unwrap();
        let with = |f: &dyn Fn(&mut FluidDomain)| {
            let mut d = channel();
            d.obstacles = vec!["Sink".into(), "Pcb".into()];
            d.thermal = Some(FluidThermal {
                initial_temperature: 293.15,
                solids: vec![
                    ConductingSolid { body: "Sink".into(), material: "Alu".into(), power: 0.0, limit: None },
                    ConductingSolid { body: "Pcb".into(), material: "Board".into(), power: 5.0, limit: Some(368.15) },
                ],
                ..FluidThermal::default()
            });
            f(&mut d);
            d
        };
        fn solids(d: &mut FluidDomain) -> &mut Vec<ConductingSolid> {
            &mut d.thermal.as_mut().unwrap().solids
        }
        m.add_fluid_domain("Flow", with(&|_| {})).unwrap();
        for (why, d) in [
            ("not an obstacle", with(&|d| d.obstacles.truncate(1))),
            ("listed twice", with(&|d| solids(d)[1].body = "Sink".into())),
            ("power not a number", with(&|d| solids(d)[0].power = f64::NAN)),
            ("limit not a number", with(&|d| solids(d)[0].limit = Some(f64::INFINITY))),
            ("limit at absolute zero", with(&|d| solids(d)[0].limit = Some(0.0))),
            ("no such material", with(&|d| solids(d)[0].material = "Nope".into())),
            ("no conductivity", with(&|d| solids(d)[0].material = "Steel".into())),
            (
                "a thermal condition on a conducting face",
                with(&|d| {
                    d.thermal.as_mut().unwrap().surfaces = vec![SurfaceThermalBoundary {
                        face: GeometryRef::face("Sink", "F1").unwrap(),
                        condition: ThermalCondition::Temperature { value: 1.0 },
                    }]
                }),
            ),
        ] {
            assert!(m.add_fluid_domain("Bad", d).is_err(), "{why}");
        }
        // the solid's material is in use, keeps fitting, and follows a rename
        assert!(m.material_users("Alu").contains(&"fluid domain `Flow`".to_string()));
        assert!(m.remove_material("Alu").is_err());
        assert!(m.replace_material("Alu", Material { density: Some(1.0), ..Material::default() }).is_err());
        m.rename(crate::rename::Kind::Material, "Alu", "Aluminium").unwrap();
        assert_eq!(m.fluid_domain("Flow").unwrap().solids().next().unwrap().material, "Aluminium");
        // it round-trips, and a domain without solids writes none
        let json = serde_json::to_string(m.fluid_domain("Flow").unwrap()).unwrap();
        assert_eq!(&serde_json::from_str::<FluidDomain>(&json).unwrap(), m.fluid_domain("Flow").unwrap());
        assert!(!serde_json::to_string(&FluidThermal::default()).unwrap().contains("solids"));
        // a solid without a limit writes none, and one written before limits reads as none
        assert!(json.contains(r#""limit":368.15"#) && json.matches("limit").count() == 1, "{json}");
        let old: ConductingSolid = serde_json::from_str(r#"{"body": "Sink", "material": "Aluminium"}"#).unwrap();
        assert_eq!(old.limit, None);
        // a limit set and cleared on a conducting solid; refused on a body that is not one
        m.set_component_limit("Flow", "Sink", Some(358.15)).unwrap();
        assert_eq!(m.fluid_domain("Flow").unwrap().solids().next().unwrap().limit, Some(358.15));
        m.set_component_limit("Flow", "Sink", None).unwrap();
        assert_eq!(m.fluid_domain("Flow").unwrap().solids().next().unwrap().limit, None);
        let e = m.set_component_limit("Flow", "Nope", Some(358.15)).unwrap_err().to_string();
        assert!(e.contains("`Nope` is not a conducting solid"), "{e}");
        assert!(m.set_component_limit("Flow", "Sink", Some(-1.0)).is_err());
        assert!(m.set_component_limit("Nope", "Sink", Some(358.15)).is_err());
    }

    #[test]
    fn thermal_contacts_join_two_conducting_solids_with_a_positive_resistance() {
        let mut m = Model::default();
        let mut hot = water();
        hot.thermal = Some(crate::Thermal { conductivity: Some(Directional::Isotropic(0.6)), specific_heat: Some(4180.0), ..Default::default() });
        m.add_material("Water", hot).unwrap();
        let alu = Material {
            density: Some(2700.0),
            thermal: Some(crate::Thermal { conductivity: Some(Directional::Isotropic(237.0)), specific_heat: Some(900.0), ..Default::default() }),
            ..Material::default()
        };
        m.add_material("Alu", alu).unwrap();
        let paste = |r: Option<f64>, t: Option<f64>, k: Option<f64>| ThermalContact {
            bodies: ["Die".into(), "Lid".into()],
            resistance: r,
            thickness: t,
            conductivity: k,
        };
        let with = |f: &dyn Fn(&mut FluidThermal)| {
            let mut d = channel();
            d.obstacles = vec!["Die".into(), "Lid".into(), "Case".into()];
            let mut t = FluidThermal {
                solids: vec![
                    ConductingSolid { body: "Die".into(), material: "Alu".into(), power: 1.0, limit: None },
                    ConductingSolid { body: "Lid".into(), material: "Alu".into(), power: 0.0, limit: None },
                ],
                contacts: vec![paste(Some(1e-5), None, None)],
                ..FluidThermal::default()
            };
            f(&mut t);
            d.thermal = Some(t);
            d
        };
        m.add_fluid_domain("Flow", with(&|_| {})).unwrap();
        m.add_fluid_domain("Pad", with(&|t| t.contacts = vec![paste(None, Some(1e-4), Some(5.0))])).unwrap();
        let pad = m.fluid_domain("Pad").unwrap().thermal.as_ref().unwrap().contacts[0].clone();
        assert!((pad.resistance() - 2e-5).abs() < 1e-20);
        for (why, d, says) in [
            ("neither form", with(&|t| t.contacts = vec![paste(None, None, None)]), "either"),
            ("both forms", with(&|t| t.contacts = vec![paste(Some(1.0), Some(1.0), Some(1.0))]), "either"),
            ("a thickness alone", with(&|t| t.contacts = vec![paste(None, Some(1.0), None)]), "either"),
            ("zero", with(&|t| t.contacts = vec![paste(Some(0.0), None, None)]), "positive"),
            ("negative", with(&|t| t.contacts = vec![paste(Some(-1.0), None, None)]), "positive"),
            ("not a number", with(&|t| t.contacts = vec![paste(Some(f64::NAN), None, None)]), "positive"),
            ("infinite", with(&|t| t.contacts = vec![paste(None, Some(1.0), Some(f64::INFINITY))]), "positive"),
            ("overflowing", with(&|t| t.contacts = vec![paste(None, Some(1e300), Some(1e-300))]), "positive"),
            ("itself", with(&|t| t.contacts[0].bodies[1] = "Die".into()), "named twice"),
            ("not conducting", with(&|t| t.contacts[0].bodies[1] = "Case".into()), "`Case` is not a conducting solid"),
            (
                "twice, either way round",
                with(&|t| t.contacts.push(ThermalContact { bodies: ["Lid".into(), "Die".into()], ..paste(Some(1.0), None, None) })),
                "already",
            ),
        ] {
            let e = m.add_fluid_domain("Bad", d).err().map(|e| e.to_string()).unwrap_or_default();
            assert!(e.contains("thermal contact `") && e.contains(says), "{why}: {e}");
        }
        // it round-trips in either form, and a domain without contacts writes none
        for name in ["Flow", "Pad"] {
            let json = serde_json::to_string(m.fluid_domain(name).unwrap()).unwrap();
            assert_eq!(&serde_json::from_str::<FluidDomain>(&json).unwrap(), m.fluid_domain(name).unwrap());
        }
        assert!(!serde_json::to_string(&FluidThermal::default()).unwrap().contains("contacts"));
    }

    #[test]
    fn domains_validate_and_refer_to_a_fluid_material() {
        let mut m = Model::default();
        assert!(m.add_fluid_domain("Flow", channel()).is_err(), "no material yet");
        m.add_material("Water", water()).unwrap();
        m.add_fluid_domain("Flow", channel()).unwrap();
        assert_eq!(m.fluid_domain("Flow").unwrap().cells().unwrap(), [20, 4, 2]);

        let bad = |f: &dyn Fn(&mut FluidDomain)| {
            let mut d = channel();
            f(&mut d);
            d
        };
        for d in [
            bad(&|d| d.spacing = 0.0),
            bad(&|d| d.spacing = f64::NAN),
            bad(&|d| d.spacing = 0.3),
            bad(&|d| d.max[0] = -1.0),
            bad(&|d| d.min[1] = f64::INFINITY),
            bad(&|d| d.spacing = 1e-300),
            bad(&|d| d.faces.z_max = FaceCondition::Wall {}),
            bad(&|d| d.faces.y_min = FaceCondition::MovingWall { velocity: [1.0, 1.0, 0.0] }),
            bad(&|d| d.faces.x_min = FaceCondition::Velocity { velocity: [f64::NAN, 0.0, 0.0], profile: Profile::Uniform }),
            bad(&|d| d.faces.x_max = FaceCondition::Pressure { pressure: f64::INFINITY, outflow: Outflow::Convective }),
            bad(&|d| d.obstacles = vec!["A".into(), "A".into()]),
            bad(&|d| d.obstacles = vec![String::new()]),
            bad(&|d| d.surfaces = vec![SurfaceBoundary { face: GeometryRef::face("B", "f").unwrap(), condition: SurfaceCondition::Wall {} }]),
            bad(&|d| {
                d.obstacles = vec!["B".into()];
                d.surfaces = vec![SurfaceBoundary { face: GeometryRef::body("B").unwrap(), condition: SurfaceCondition::Wall {} }];
            }),
            bad(&|d| d.tolerance = Some(0.0)),
            bad(&|d| d.acceleration = [0.0, f64::NAN, 0.0]),
        ] {
            assert!(m.add_fluid_domain("Bad", d.clone()).is_err(), "{d:?}");
        }
        let mut solid = water();
        solid.fluid = None;
        m.add_material("Steel", Material { elastic: None, density: Some(7.85e-9), ..Material::default() }).unwrap();
        assert!(m.add_fluid_domain("Dry", FluidDomain { material: "Steel".into(), ..channel() }).is_err());
        assert!(m.replace_material("Water", solid).is_err(), "a domain needs the fluid behaviour");
        assert!(matches!(m.remove_material("Water"), Err(Error::InUse(_))));

        // A step refers to the domain, which then cannot be deleted.
        let fluid = Procedure::Fluid {
            domain: "Flow".into(),
            end_time: Some(1.0),
            steps: None,
            time_step: None,
            frames: 4,
            magic: None,
            collision: None,
            turbulence: None,
            steady: None,
        };
        m.add_step(Step::new("Run", fluid.clone()), None).unwrap();
        assert!(matches!(m.remove_fluid_domain("Flow"), Err(Error::InUse(_))));
        // The collision operator: named in lower case, left out when not given, nothing else.
        let regularized: Procedure =
            serde_json::from_value(serde_json::json!({"type": "fluid", "domain": "Flow", "end_time": 1.0, "collision": "regularized"})).unwrap();
        assert!(matches!(regularized, Procedure::Fluid { collision: Some(crate::Collision::Regularized), .. }));
        assert!(!serde_json::to_string(&fluid).unwrap().contains("collision"));
        assert!(
            serde_json::from_value::<Procedure>(serde_json::json!({"type": "fluid", "domain": "Flow", "end_time": 1.0, "collision": "cumulant"}))
                .is_err()
        );
        assert_eq!(serde_json::from_str::<Procedure>(&serde_json::to_string(&regularized).unwrap()).unwrap(), regularized);
        // The JSON form reads back to the same model, through the same checks.
        let text = serde_json::to_string(&m).unwrap();
        assert_eq!(serde_json::from_str::<Model>(&text).unwrap(), m);
        assert!(serde_json::from_str::<Model>(&text.replace("\"spacing\":0.5", "\"spacing\":0.3")).is_err());
        // Fluid steps take nothing structural, and need a domain and one length.
        let mut listed = Step::new("Listed", fluid.clone());
        listed.output.insert(crate::OutputVariable::Stress);
        assert!(m.add_step(listed, None).is_err());
        for bad in [
            Procedure::Fluid {
                domain: "Nope".into(),
                end_time: Some(1.0),
                steps: None,
                time_step: None,
                frames: 4,
                magic: None,
                collision: None,
                turbulence: None,
                steady: None,
            },
            Procedure::Fluid {
                domain: "Flow".into(),
                end_time: Some(1.0),
                steps: Some(5),
                time_step: None,
                frames: 4,
                magic: None,
                collision: None,
                turbulence: None,
                steady: None,
            },
            Procedure::Fluid {
                domain: "Flow".into(),
                end_time: None,
                steps: None,
                time_step: None,
                frames: 4,
                magic: None,
                collision: None,
                turbulence: None,
                steady: None,
            },
            Procedure::Fluid {
                domain: "Flow".into(),
                end_time: Some(-1.0),
                steps: None,
                time_step: None,
                frames: 4,
                magic: None,
                collision: None,
                turbulence: None,
                steady: None,
            },
            Procedure::Fluid {
                domain: "Flow".into(),
                end_time: None,
                steps: Some(0),
                time_step: None,
                frames: 4,
                magic: None,
                collision: None,
                turbulence: None,
                steady: None,
            },
            Procedure::Fluid {
                domain: "Flow".into(),
                end_time: Some(1.0),
                steps: None,
                time_step: Some(0.0),
                frames: 4,
                magic: None,
                collision: None,
                turbulence: None,
                steady: None,
            },
            Procedure::Fluid {
                domain: "Flow".into(),
                end_time: Some(1.0),
                steps: None,
                time_step: None,
                frames: 0,
                magic: None,
                collision: None,
                turbulence: None,
                steady: None,
            },
            Procedure::Fluid {
                domain: "Flow".into(),
                end_time: Some(1.0),
                steps: None,
                time_step: None,
                frames: 4,
                magic: Some(0.0),
                collision: None,
                turbulence: None,
                steady: None,
            },
            Procedure::Fluid {
                domain: "Flow".into(),
                end_time: Some(1.0),
                steps: None,
                time_step: None,
                frames: 4,
                magic: Some(f64::NAN),
                collision: None,
                turbulence: None,
                steady: None,
            },
        ] {
            assert!(m.add_step(Step::new("Bad", bad.clone()), None).is_err(), "{bad:?}");
        }
        // The steady-state rule's numbers are checked too.
        let steady = |st: crate::Steady| Procedure::Fluid {
            domain: "Flow".into(),
            end_time: Some(1.0),
            steps: None,
            time_step: None,
            frames: 4,
            magic: None,
            collision: None,
            turbulence: None,
            steady: Some(st),
        };
        let rule = crate::Steady { monitors: vec![], tolerance: 1e-3, window: None, min_time: None, divergence: None, averaging: None };
        let monitored = crate::Steady {
            monitors: vec![
                crate::FlowMonitor::MassFlow { face: BoxFace::XMax },
                crate::FlowMonitor::Force {
                    faces: vec![FluidBoundary::Box(BoxFace::YMin), FluidBoundary::Surface(crate::GeometryRef::face("Ball", "F3").unwrap())],
                },
                crate::FlowMonitor::Heat { faces: vec![FluidBoundary::Box(BoxFace::XMin)] },
                crate::FlowMonitor::MeanTemperature {},
            ],
            tolerance: 1e-4,
            window: Some(0.1),
            min_time: Some(0.2),
            divergence: Some(0.5),
            averaging: None,
        };
        m.add_step(Step::new("Steady", steady(monitored)), None).unwrap();
        // The rule and its monitors read back from JSON to the same model.
        let text = serde_json::to_string(&m).unwrap();
        assert!(text.contains("\"mean_temperature\"") && text.contains("\"F3\""), "{text}");
        assert_eq!(serde_json::from_str::<Model>(&text).unwrap(), m);
        m.remove_step("Steady").unwrap();
        let force = |n: usize| crate::FlowMonitor::Force { faces: vec![FluidBoundary::Box(BoxFace::XMin); n] };
        for bad in [
            crate::Steady { tolerance: 0.0, ..rule.clone() },
            crate::Steady { tolerance: f64::NAN, ..rule.clone() },
            crate::Steady { tolerance: 1.0, ..rule.clone() },
            crate::Steady { window: Some(-1.0), ..rule.clone() },
            crate::Steady { window: Some(f64::INFINITY), ..rule.clone() },
            crate::Steady { min_time: Some(-1.0), ..rule.clone() },
            crate::Steady { min_time: Some(f64::NAN), ..rule.clone() },
            crate::Steady { divergence: Some(1e-4), ..rule.clone() },
            crate::Steady { divergence: Some(f64::INFINITY), ..rule.clone() },
            crate::Steady { monitors: vec![force(0)], ..rule.clone() },
            crate::Steady { monitors: vec![force(1001)], ..rule.clone() },
            crate::Steady { monitors: vec![crate::FlowMonitor::MeanTemperature {}; 101], ..rule.clone() },
        ] {
            assert!(m.add_step(Step::new("Bad", steady(bad.clone())), None).is_err(), "{bad:?}");
        }
        m.remove_step("Run").unwrap();
        m.remove_fluid_domain("Flow").unwrap();
        m.remove_material("Water").unwrap();
    }

    #[test]
    fn boundaries_are_set_one_at_a_time() {
        let mut m = Model::default();
        m.add_material("Water", water()).unwrap();
        let mut d = channel();
        d.obstacles = vec!["Ball".into()];
        m.add_fluid_domain("Flow", d).unwrap();
        m.set_fluid_boundary("Flow", FluidBoundary::Box(BoxFace::YMax), FaceCondition::Slip {}).unwrap();
        assert_eq!(m.fluid_domain("Flow").unwrap().faces.y_max, FaceCondition::Slip {});
        let face = GeometryRef::face("Ball", "F1").unwrap();
        let spin = FaceCondition::MovingWall { velocity: [0.0, 1.0, 0.0] };
        m.set_fluid_boundary("Flow", FluidBoundary::Surface(face.clone()), spin).unwrap();
        assert_eq!(m.fluid_domain("Flow").unwrap().surfaces.len(), 1);
        m.set_fluid_boundary("Flow", FluidBoundary::Surface(face.clone()), FaceCondition::Wall {}).unwrap();
        assert!(m.fluid_domain("Flow").unwrap().surfaces.is_empty(), "walls are the default");
        // Wrong kinds and broken pairs are refused and change nothing.
        let before = m.clone();
        assert!(m.set_fluid_boundary("Flow", FluidBoundary::Box(BoxFace::ZMin), FaceCondition::Wall {}).is_err());
        assert!(m.set_fluid_boundary("Flow", FluidBoundary::Surface(face), FaceCondition::Slip {}).is_err(), "obstacle faces are walls");
        assert!(m.set_fluid_boundary("Nope", FluidBoundary::Box(BoxFace::XMin), FaceCondition::Wall {}).is_err());
        assert_eq!(m, before);

        // The JSON forms.
        let at: FluidBoundary = serde_json::from_value(serde_json::json!("x_min")).unwrap();
        assert_eq!(at, FluidBoundary::Box(BoxFace::XMin));
        let at: FluidBoundary = serde_json::from_value(serde_json::json!({"body": "Ball", "face": "F1"})).unwrap();
        assert!(matches!(at, FluidBoundary::Surface(_)));
        assert!(serde_json::from_value::<FaceCondition>(serde_json::json!({"type": "wall", "velocity": [1, 0, 0]})).is_err());
        // A pressure face is uniform unless it says otherwise, and writes `outflow` only then.
        let uniform = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Uniform };
        assert_eq!(serde_json::from_value::<FaceCondition>(serde_json::json!({"type": "pressure", "pressure": 0})).unwrap(), uniform);
        assert_eq!(serde_json::to_value(uniform).unwrap(), serde_json::json!({"type": "pressure", "pressure": 0.0}));
        let convective = serde_json::json!({"type": "pressure", "pressure": 0.0, "outflow": "convective"});
        let c: FaceCondition = serde_json::from_value(convective.clone()).unwrap();
        assert_eq!(c, FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Convective });
        assert_eq!(serde_json::to_value(c).unwrap(), convective);
        assert!(serde_json::from_value::<FaceCondition>(serde_json::json!({"type": "pressure", "pressure": 0, "outflow": "sponge"})).is_err());
    }
}
