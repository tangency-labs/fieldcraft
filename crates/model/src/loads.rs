//! Loads. Which steps a load is active in is listed explicitly by each step ([`crate::steps`]).
//!
//! In a 2D model (plane stress, plane strain, axisymmetric sections) loads lie in the XY plane:
//! surface loads act on the element edges of the set, forces are totals over the thickness (or
//! over the full circumference for axisymmetric models), and a vector with a Z component is an
//! error.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

use crate::bcs::Target;
use crate::{Error, Model, Result, context, finite, insert_new, nonnegative, not_in_use, positive, remove_existing, replace_existing, require};

/// The Stefan–Boltzmann constant σ in W/(m²·K⁴) (CODATA 2018, exact since the SI redefinition
/// of 2019): the model is stored in SI, so it applies as it is.
pub const STEFAN_BOLTZMANN: f64 = 5.670_374_419e-8;

/// Radiation from a surface to its surroundings: a heat flux ε σ (T⁴ − T_amb⁴) out of it
/// ([`STEFAN_BOLTZMANN`]). The surface sees only the surroundings: radiation between
/// surfaces (view factors) is not modelled.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SurfaceRadiation {
    /// The surface's emissivity, above 0 and at most 1.
    pub emissivity: f64,
    /// The surroundings' temperature (kelvin, absolute).
    pub ambient: f64,
}

impl SurfaceRadiation {
    pub fn validate(&self) -> Result<()> {
        if !(self.emissivity > 0.0 && self.emissivity <= 1.0) {
            return Err(Error::Invalid(format!("an emissivity must be a number above 0 and at most 1, got {}", self.emissivity)));
        }
        nonnegative("ambient temperature (absolute)", self.ambient).map(|_| ())
    }
}

/// How a force or moment on a set is shared among the set's nodes.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Distribution {
    /// The vector is the total. On a set with element faces it is split by the faces' consistent
    /// nodal weights (area-weighted, [`Model::surface_weights`]); on a set of nodes only, equally.
    /// The total is the same on every mesh of the set.
    #[default]
    Total,
    /// The vector is applied at every node of the set, so the total grows with the number of
    /// nodes: it depends on the mesh.
    PerNode,
}

/// A load.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case", deny_unknown_fields)]
pub enum Load {
    /// A force on the target's nodes, shared as `distribution` says (a reference point takes it
    /// whole).
    Force {
        target: Target,
        vector: [f64; 3],
        #[serde(default)]
        distribution: Distribution,
    },
    /// A moment on the target's nodes, shared as `distribution` says. Solid elements have no
    /// rotations, so a moment on a solid goes through a coupled reference point.
    Moment {
        target: Target,
        vector: [f64; 3],
        #[serde(default)]
        distribution: Distribution,
    },
    /// A force per area on the set's element faces, as consistent nodal forces.
    Traction { set: String, vector: [f64; 3] },
    /// Uniform pressure on the set's element faces, positive pushing into the elements (against
    /// the outward normal).
    Pressure { set: String, magnitude: f64 },
    /// Gravity: an acceleration of the mass of the set's elements, or of the whole model when no
    /// set is given.
    Gravity {
        acceleration: [f64; 3],
        #[serde(default, skip_serializing_if = "Option::is_none")]
        set: Option<String>,
    },
    /// A force per volume on the set's elements.
    BodyForce { set: String, vector: [f64; 3] },
    /// The centrifugal load of a steady rotation at `angular_velocity` (radians per time unit)
    /// about the axis through `point` along `axis`: a force per volume ρ ω² r on the mass of the
    /// set's elements (or the whole model), r the distance from the axis, pointing away from it.
    Centrifugal {
        angular_velocity: f64,
        point: [f64; 3],
        axis: [f64; 3],
        #[serde(default, skip_serializing_if = "Option::is_none")]
        set: Option<String>,
    },
    /// Heat flux per area into the set's element faces (positive heats the body).
    SurfaceFlux { set: String, flux: f64 },
    /// Convection from the set's element faces to a surrounding medium at `sink` temperature:
    /// heat flux h (T_sink − T) into the body, with film coefficient `film` (h).
    Convection { set: String, film: f64, sink: f64 },
    /// Radiation from the set's element faces to surroundings at `ambient` (an absolute
    /// temperature): heat flux ε σ (T⁴ − T_amb⁴) out of the body, with surface emissivity
    /// `emissivity` (0 < ε ≤ 1) and σ [`STEFAN_BOLTZMANN`]. The faces see only the ambient:
    /// radiation between surfaces (view factors) is not modelled.
    Radiation { set: String, emissivity: f64, ambient: f64 },
    /// Heat generated per volume in the set's elements.
    HeatSource { set: String, power: f64 },
    /// A temperature field at the target's nodes in a stress step: it strains the material
    /// thermally (ε = α (T − T_ref)); nodes it does not reach stay at their material's reference
    /// temperature.
    Temperature { target: Target, value: f64 },
}

/// The centrifugal acceleration ω² r at `x` of a rotation at `angular_velocity` about the axis
/// through `point` along `axis` (any length, not zero): r is the part of `x − point` normal to
/// the axis, so the acceleration points away from the axis.
pub fn centrifugal_acceleration(angular_velocity: f64, point: [f64; 3], axis: [f64; 3], x: [f64; 3]) -> [f64; 3] {
    let len2 = axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2];
    let d = [x[0] - point[0], x[1] - point[1], x[2] - point[2]];
    let along = if len2 > 0.0 { (d[0] * axis[0] + d[1] * axis[1] + d[2] * axis[2]) / len2 } else { 0.0 };
    let w2 = angular_velocity * angular_velocity;
    [0, 1, 2].map(|k| w2 * (d[k] - along * axis[k]))
}

impl Load {
    /// `force`, `moment`, `traction`, `pressure`, `gravity` or `body_force`.
    pub fn kind(&self) -> &'static str {
        match self {
            Load::Force { .. } => "force",
            Load::Moment { .. } => "moment",
            Load::Traction { .. } => "traction",
            Load::Pressure { .. } => "pressure",
            Load::Gravity { .. } => "gravity",
            Load::BodyForce { .. } => "body_force",
            Load::Centrifugal { .. } => "centrifugal",
            Load::SurfaceFlux { .. } => "surface_flux",
            Load::Convection { .. } => "convection",
            Load::Radiation { .. } => "radiation",
            Load::HeatSource { .. } => "heat_source",
            Load::Temperature { .. } => "temperature",
        }
    }

    /// Whether it is a heat transfer load (flux, convection, radiation, heat source), active in
    /// heat transfer and coupled temperature–displacement steps.
    pub fn is_thermal(&self) -> bool {
        matches!(self, Load::SurfaceFlux { .. } | Load::Convection { .. } | Load::Radiation { .. } | Load::HeatSource { .. })
    }

    /// The set the load acts on, if any.
    pub fn set(&self) -> Option<&str> {
        match self {
            Load::Force { target, .. } | Load::Moment { target, .. } | Load::Temperature { target, .. } => target.set(),
            Load::Pressure { set, .. }
            | Load::BodyForce { set, .. }
            | Load::Traction { set, .. }
            | Load::SurfaceFlux { set, .. }
            | Load::Convection { set, .. }
            | Load::Radiation { set, .. }
            | Load::HeatSource { set, .. } => Some(set),
            Load::Gravity { set, .. } | Load::Centrifugal { set, .. } => set.as_deref(),
        }
    }

    /// The reference point the load acts on, if any.
    pub fn reference_point(&self) -> Option<&str> {
        match self {
            Load::Force { target, .. } | Load::Moment { target, .. } | Load::Temperature { target, .. } => target.reference_point(),
            _ => None,
        }
    }

    pub fn validate(&self) -> Result<()> {
        let vector = match self {
            Load::Force { vector, .. } | Load::Moment { vector, .. } | Load::BodyForce { vector, .. } | Load::Traction { vector, .. } => vector,
            Load::Gravity { acceleration, .. } => acceleration,
            Load::Pressure { magnitude, .. } => {
                finite("pressure", *magnitude)?;
                return Ok(());
            }
            Load::Centrifugal { angular_velocity, point, axis, .. } => {
                let w = finite("angular velocity", *angular_velocity)?;
                if w == 0.0 || !(w * w).is_finite() {
                    return Err(Error::Invalid(format!(
                        "a centrifugal load's angular velocity must be a non-zero number whose square is finite, got {w}"
                    )));
                }
                for c in point.iter().chain(axis) {
                    finite("centrifugal axis", *c)?;
                }
                if axis.iter().all(|c| *c == 0.0) {
                    return Err(Error::Invalid("a centrifugal load's axis direction must not be zero".into()));
                }
                return Ok(());
            }
            Load::SurfaceFlux { flux: v, .. } | Load::HeatSource { power: v, .. } => {
                if finite(self.kind(), *v)? == 0.0 {
                    return Err(Error::Invalid(format!("a {} must not be zero", self.kind())));
                }
                return Ok(());
            }
            Load::Convection { film, sink, .. } => {
                positive("film coefficient", *film)?;
                finite("sink temperature", *sink)?;
                return Ok(());
            }
            Load::Radiation { emissivity, ambient, .. } => return SurfaceRadiation { emissivity: *emissivity, ambient: *ambient }.validate(),
            Load::Temperature { target, value } => {
                finite("temperature", *value)?;
                if let Target::ReferencePoint(r) = target {
                    return Err(Error::Invalid(format!("reference points have no temperature (`{r}`); apply it to a set")));
                }
                return Ok(());
            }
        };
        for c in vector {
            finite(self.kind(), *c)?;
        }
        if vector.iter().all(|c| *c == 0.0) {
            return Err(Error::Invalid(format!("a {} vector must not be zero", self.kind())));
        }
        Ok(())
    }
}

impl Model {
    pub fn loads(&self) -> &BTreeMap<String, Load> {
        &self.loads
    }

    pub fn load(&self, name: &str) -> Option<&Load> {
        self.loads.get(name)
    }

    fn check_load(&self, name: &str, load: &Load) -> Result<()> {
        load.validate().map_err(|e| context("load", name, e))?;
        if let Some(s) = load.set() {
            require(&self.sets, "set", s).map_err(|e| context("load", name, e))?;
        }
        if let Some(r) = load.reference_point() {
            require(&self.reference_points, "reference point", r).map_err(|e| context("load", name, e))?;
        }
        Ok(())
    }

    pub fn add_load(&mut self, name: &str, load: Load) -> Result<()> {
        self.check_load(name, &load)?;
        insert_new(&mut self.loads, "load", name, load)
    }

    pub fn replace_load(&mut self, name: &str, load: Load) -> Result<()> {
        self.check_load(name, &load)?;
        replace_existing(&mut self.loads, "load", name, load)
    }

    /// The steps that list a load, as `step \`name\`` strings.
    pub fn load_users(&self, name: &str) -> Vec<String> {
        self.steps.iter().filter(|s| s.loads.contains(name)).map(|s| format!("step `{}`", s.name)).collect()
    }

    /// Remove a load that no step lists.
    pub fn remove_load(&mut self, name: &str) -> Result<()> {
        not_in_use("load", name, self.load_users(name))?;
        remove_existing(&mut self.loads, "load", name).map(|_| ())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn heat_loads_validate() {
        let set = || "s".to_string();
        assert!(Load::SurfaceFlux { set: set(), flux: 2.0 }.validate().is_ok());
        assert!(Load::SurfaceFlux { set: set(), flux: 0.0 }.validate().is_err());
        assert!(Load::HeatSource { set: set(), power: f64::NAN }.validate().is_err());
        assert!(Load::Convection { set: set(), film: 0.0, sink: 1.0 }.validate().is_err());
        assert!(Load::Convection { set: set(), film: 1.0, sink: f64::INFINITY }.validate().is_err());
        assert!(Load::Temperature { target: Target::ReferencePoint("RP".into()), value: 1.0 }.validate().is_err());
        assert!(Load::Temperature { target: Target::Set(set()), value: -40.0 }.validate().is_ok());
        let l = Load::Convection { set: set(), film: 10.0, sink: 20.0 };
        assert!(l.is_thermal() && l.set() == Some("s") && l.kind() == "convection");
        assert!(!Load::Temperature { target: Target::Set(set()), value: 1.0 }.is_thermal());
        let rad = |emissivity, ambient| Load::Radiation { set: set(), emissivity, ambient }.validate();
        assert!(rad(0.8, 293.15).is_ok() && rad(1.0, 0.0).is_ok());
        for (e, t) in [(0.0, 300.0), (1.5, 300.0), (f64::NAN, 300.0), (0.5, -1.0), (0.5, f64::INFINITY)] {
            assert!(rad(e, t).is_err(), "{e} {t}");
        }
        let l = Load::Radiation { set: set(), emissivity: 0.9, ambient: 300.0 };
        assert!(l.is_thermal() && l.set() == Some("s") && l.kind() == "radiation");
    }
}
